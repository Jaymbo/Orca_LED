import os
import re
import time
import streamlit as st
from pathlib import Path
import pipeline
from LED_extraction import LEDExtractor
from csv_to_viz import extract
from visualization import MoleculeVisualizer
import logging
from xlsx_to_sdf import SdfXyzMerger
from rdkit import Chem
import cProfile
import pstats
import numpy as np
import asyncio
from database import Database

BASE_PATH = Path(Path(__file__).resolve().parent.parent / "calculations")
if not BASE_PATH.exists():
    BASE_PATH.mkdir(parents=True, exist_ok=True)
open_topic = ""
state = 0
file_cache ={}

def profile(func):
    def wrapper(*args, **kwargs):
        profiler = cProfile.Profile()
        profiler.enable()
        result = func(*args, **kwargs)
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.print_stats(40)  # Print only the first 40 lines
        return result
    return wrapper

class Dashboard:
    class Visualizer:
        @staticmethod
        def visualize_xyz(molecule_path: Path):
            viewer = MoleculeVisualizer.render_xyz(molecule_path)
            st.session_state.previous_xyz = viewer  # Nur wenn benötigt
            st.components.v1.html(viewer._make_html(), height=600)

        @staticmethod
        def visualize_in_3Dmol(mols, viz, state: int):
            MoleculeVisualizer.show_led_analysis(
                mols=mols,
                viz=viz,
                state=state,
            )

    @staticmethod
    @st.fragment(run_every="600s")
    @profile
    def upload_file_and_start_calculation():

        def chose_topic():
            # TODO für wenn topic noch nicht ausgewählt wurde und dann nichts neues eingegeben wird muss neues topic gelöscht werden mit inhalt
            # besser noch einfach nur einmal hochladen und verschieben und nicht bei jedem mal neu hochladen
            neues_topic_text = "Neues Topic erstellen"
            topic_options = [f for f in os.listdir(BASE_PATH) if os.path.isdir(os.path.join(BASE_PATH, f))]
            topic_options.append(neues_topic_text)
            index = topic_options.index(open_topic) if open_topic in topic_options else len(topic_options) - 1
            topic = st.selectbox("Topic auswählen", topic_options, index=index)
            if topic == "Neues Topic erstellen":
                new_topic = st.text_input("Neuen Topic-Namen eingeben")
                if new_topic:
                    new_topic_path = BASE_PATH / new_topic
                    new_topic_path.mkdir(parents=True, exist_ok=True)
                    topic_options.append(new_topic)
                    index = topic_options.index(new_topic)
                    topic = new_topic
                    st.success("Neues Topic erfolgreich erstellt.")
                    new_topic_text_path = BASE_PATH / neues_topic_text    
                    if new_topic_text_path.exists():
                        for item in new_topic_text_path.iterdir():
                            if item.is_dir():
                                for sub_item in item.iterdir():
                                    sub_item.unlink()
                                item.rmdir()
                            else:
                                item.unlink()
                        new_topic_text_path.rmdir()
            return topic

        def upload_file(file_paths=[]):
            
            class FileHandler:
                @staticmethod
                def handle_file_upload(topic, uploaded_file) -> Path:
                    """Verarbeitet hochgeladene Dateien und gibt gespeicherten Pfad zurück"""
                    try:
                        clean_name = uploaded_file.name.translate(str.maketrans(" ", "_", "!$%&()=+,-/:;<=>?@[\]^`{|}~"))
                        save_path = BASE_PATH / f"{topic}/{Path(clean_name).stem}/{clean_name}"
                        save_path.parent.mkdir(parents=True, exist_ok=True)
                        save_path.write_bytes(uploaded_file.getbuffer())
                        logging.info(f"File uploaded: {uploaded_file.name} with topic: {topic}")
                        return save_path
                    except Exception as e:
                        logging.error(f"Upload failed: {str(e)}")
                        raise

                @staticmethod
                def sanitize_saving(topic, content: list[str], name: str):
                    clean_name = name.translate(str.maketrans(" ", "_", "!$%&()=+,-/:;<=>?@[\]^`{|}~"))
                    paths = [
                        BASE_PATH / f"{topic}/{Path(clean_name).stem}_{i+1}/{Path(clean_name).stem}_{i+1}.mol2"
                        for i in range(len(content))
                    ]
                    for path, text in zip(paths, content):
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(text)
                    return paths
            uploaded_files = st.file_uploader("Wählen Sie eine Datei aus", type=["xyz", "mol2", "sdf"], accept_multiple_files=True)
            
            if uploaded_files:
                for f in uploaded_files:
                    if f.name.endswith(".sdf"):
                        topic = f.name.split(".")[0]
                        if not os.path.exists(BASE_PATH / topic):
                            os.makedirs(BASE_PATH / topic)
                        path = FileHandler.handle_file_upload(topic, f)
                        file_paths.extend(FileHandler.sanitize_saving(topic, MoleculeVisualizer.convert_sdf_to_mol2(path), f.name))
                        #move path to topic folder
                        path.rename(BASE_PATH / topic / path.name)
                        #remove empty folder
                        path.parent.rmdir()

                    else:
                        topic = chose_topic()
                        file_paths.append(FileHandler.handle_file_upload(topic, f))
                Dashboard.Visualizer.visualize_xyz(file_paths)
            return file_paths
    
        st.title("Neue Berechnung starten")
        file_paths = upload_file()
        header_input = st.text_area("Zusätzliche Header-Einstellungen für die .inp-Datei (optional)", """! DLPNO-CCSD(T) def2-svp def2-svp/C DEF2/J RIJCOSX tightSCF normalPNO LED

    %maxcore 160000

    %mdci
        MaxIter 200
    end""")
        if st.button("Berechnung starten"):
            if file_paths:
                st.info("Die Berechnung wurde in Auftrag gegeben...")
                async def process_file(file_path):
                    pipeline.ORCAInputFileCreator(str(file_path), header_input).create_inp_files(file_cache)

                async def process_all_files():
                    tasks = [process_file(file_path) for file_path in file_paths]
                    await asyncio.gather(*tasks)

                asyncio.run(process_all_files())
                st.success("Die Berechnung wurde vorbereitet.")
            else:
                st.error("Bitte wählen Sie eine Datei aus.")

    @staticmethod
    @st.fragment(run_every="600s")
    @profile
    def check_progress_of_all_jobs():
        class JobHandler:
            @staticmethod
            def check_orca_termination(content):
                return "****ORCA TERMINATED NORMALLY****" in content

            @staticmethod
            def get_progress_of_job(output, path):
                runtime = None
                error_status = JobHandler.get_error_file(path)
                
                status, runtime = JobHandler.check_slurm_job_status_and_duration(path)
                if status != "Not Finished":
                    return status, runtime

                if JobHandler.check_orca_termination(output):
                    return "Progress: 100%", runtime
                
                if error_status is not None:
                    return error_status, None
                count = 0
                output = output.split("\n")
                if "INITIAL GUESS DONE" in output:
                    count += 1
                for line in output:
                    if "TIMINGS" in line:
                        count += 1
                return f"Progress: {int(count *50)}%", runtime

            @staticmethod
            def get_error_file(path):
                total_path = Path(f"{path}_err.err")
                if os.path.exists(total_path):
                    try:
                        content = total_path.read_text()
                        if "multiplicity" in content:
                            return "multiplicity"
                        if "OUT OF MEMORY ERROR!" in content:
                            return "memory"
                        if "Segmentation fault" in content:
                            return "Seg fault"
                        if "Wrong syntax in xyz coordinates" in content:
                            return "Syntax"
                        if "Tool-Scanner" in content:
                            return "Scanner"
                        if "CANCELLED AT" in content:
                            return "CANCELLED"
                        if "mpirun noticed that process" in content:
                            return "mpirun"
                        if "CalcSigma" in content:
                            return "CalcSigma"
                        if "out of memory" in content:
                            return "memory"
                        if "CANCELLED" in content:
                            return "CANCELLED"
                        if "aborting the run" in content:
                            return "aborted"
                    except Exception as e:
                        logging.error(f"Error reading error file: {str(e)}")
                        return f"Error: {str(e)}"
                return None

            @staticmethod
            def check_slurm_job_status_and_duration(base_path: Path) -> tuple:
                """Gibt Status und Laufzeit für SLURM-Job zurück"""
                slurm_output_path = base_path.with_name(base_path.name + "_out.out")
                
                if not slurm_output_path.exists():
                    return "Not Finished", None

                status = "Not Finished"
                runtime = None
                content = slurm_output_path.read_text()
                
                if "DUE TO TIME LIMIT" in content:
                    status = "Failed: Time Limit"
                elif "CANCELLED" in content:
                    status = "Cancelled"
                elif "CPU Utilized:" in content:
                    time_parts = content.split("CPU Utilized:")[1].split("\n")[0].split(":")[0:3]
                    runtime = ":".join(time_parts).strip()
                return status, runtime

            @staticmethod
            def get_color_and_progress(progress: str) -> tuple:
                """Gibt Farbe und Prozentwert für Progress-Balken zurück"""
                if progress == "Not Started":
                    return "grey", 100
                if "Progress:" not in progress:
                    return "red", 100
                if "100%" in progress:
                    return "green", 100
                
                progress_value = float(progress.replace("Progress:", "").replace("%", "").strip())
                return "orange", progress_value

            @staticmethod
            def energy_extraction(context):
                matches = list(re.finditer(r"FINAL SINGLE POINT ENERGY \s+([-+]?\d+\.\d+)", context))
                if not matches:
                    return None
                return float(matches[-1].group(1))

        def check_progress_of_single_topic(topic_path):
            def check_progress_of_single_file(path):
                if path.is_dir():
                    output_file = path / f"{path.stem}.out"
                    if output_file.exists():
                        lines = output_file.read_text()

                        progress, runtime = JobHandler.get_progress_of_job(lines, path / path.stem)
                        jobs_progress = (progress, runtime, path.stem, JobHandler.energy_extraction(lines))
                    else:
                        jobs_progress = ("Not Started", None, path.stem, None)
                return jobs_progress

            subtopics = {}
            for subfolder_name in os.listdir(topic_path):
                subtopics[subfolder_name] = {}
                subfolder_path = Path(topic_path) / subfolder_name
                if os.path.isdir(subfolder_path):
                    subfolders = sorted(
                        [f for f in subfolder_path.iterdir() if f.is_dir()],
                        key=lambda x: ("subsys_" in x.name, x.name)
                    ) or sorted(
                        [f for f in subfolder_path.iterdir() if f.is_dir()],
                        key=lambda x: ("fragment_" in x.name, x.name)
                    )
                    check = any(folder.name.startswith(("fragment_", "subsys_")) for folder in subfolders)

                    if check:
                        jobs_progress = {}
                        for subsubfolder_name in subfolders:
                            job_path = subfolder_path / subsubfolder_name
                            jobs_info = check_progress_of_single_file(job_path)
                            jobs_progress[subsubfolder_name] = jobs_info
                        if jobs_progress:
                            subtopics[subfolder_name] = (jobs_progress, subfolder_path)

                if subtopics[subfolder_name] == {}:
                    del subtopics[subfolder_name]
            return subtopics

        def update_dashboard(topics_progress):
            global open_topic, state

            def update_single_topic(jobs, folder, mols, viz, check_sdf=True):
                total_jobs = len(jobs)
                completed_jobs = sum(1 for progress, _, _, _ in jobs.values() if "Progress: 100%" == progress)
                # not_started_jobs = sum(1 for progress, _, _, _ in jobs.values() if progress == "Not Started")

                if completed_jobs == total_jobs:
                    # st.markdown(f"<div style='background-color: green; height: 10px; width: 100%'></div>", unsafe_allow_html=True)
                    energy = 2 * jobs[list(jobs.keys())[0]][3]
                    for i in range(completed_jobs):
                        energy -= jobs[list(jobs.keys())[i]][3]
                    # st.text(f"{energy * 627.509474:.6f} kcal/mol")
                    
                    LEDExtractor(folder).extract_LED_energy()
                    xyz_file = folder / f"{folder.name}.xyz"
                    extract(folder)
                    if check_sdf:
                        mols, viz = SdfXyzMerger(xyz_file, folder, mols, viz).run()
                    # nur wenn drauf geklickt wird
                    # if st.button(f"Visualisierung für {folder.name} erstellen"):
                    #     Dashboard.Visualizer.visualize_in_3Dmol(Path(folder), Path(folder) / "viz.py")

                    # st.text(f"Alle Berechnungen abgeschlossen!")
                    # logging.info(f"All calculations completed for topic: {topic}")
                # elif not_started_jobs == total_jobs:
                    # st.markdown(f"<div style='background-color: grey; height: 10px; width: 100%'></div>", unsafe_allow_html=True)
                    # st.text(f"Alle Berechnungen noch nicht gestartet.")
                    # logging.info(f"No calculations started for topic: {topic}")
                # else:
                #     num_columns = 7
                #     cols = st.columns(num_columns)

                #     for i, (subfolder_name, (progress, runtime, path, energy)) in enumerate(jobs.items()):
                #         subfolder_name = subfolder_name.replace("fragment_", "")
                #         color, progress_value = JobHandler.get_color_and_progress(progress)
                #         if len(subfolder_name) > 6:
                #             subfolder_name = subfolder_name[:6] + "..."

                #         col = cols[i % num_columns]
                    #     with col:
                    #         st.text(f"{subfolder_name}:\n{progress_value}%\n{runtime if runtime else 'No Time'}")
                    #         st.markdown(f"<div style='background-color: {color}; height: 10px; width: {progress_value}%'></div>", unsafe_allow_html=True)

                    # st.text(f"Fortschritt des Topics: {completed_jobs}/{total_jobs} abgeschlossen")
                    # logging.info(f"Progress for topic {topic}: {completed_jobs}/{total_jobs} completed.")
            
                return mols, viz
            def download_sdf_file():
                viz_file = Path(str(sdf_file).split('.')[0] + f".py")
                print(viz_file)
                if sdf_file.exists():
                    with open(sdf_file, "rb") as file:
                        st.download_button(label="Download SDF", data=file, file_name=sdf_file.name)
                if viz_file.exists():
                    with open(viz_file, "rb") as file:
                        st.download_button(label="Download viz.py", data=file, file_name=viz_file.name)

            def load_sdf():
                if not sdf_file.exists() or sdf_file.stat().st_size == 0:
                    mols = []
                else:
                    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
                    mols = [mol for mol in suppl if mol is not None]
                return mols

            def load_mols():
                if not sdf_file.exists():
                    sdf_file.touch()
                    mols = []
                else:
                    mols = load_sdf()
                return mols
    
            def load_viz():
                if not viz_file.exists():
                    viz_file.touch()
                    lines = ""
                else:
                    with open(viz_file, "r") as file:
                        lines = file.read()
                return lines
                
            def save_mols():
                writer = Chem.SDWriter(sdf_file)
                for mol in mols:
                    writer.write(mol)
                writer.close()
            
            def save_viz():
                with open(viz_file, "w") as file:
                    file.writelines(viz)

            st.title('ORCA-Status')
            for topic_name, topic in topics_progress.items():
                st.subheader(f"Topic: {topic_name}")
                if st.button(f"Fortschritt für {topic_name} anzeigen") or open_topic == topic_name:
                    open_topic = topic_name
                    sdf_file = BASE_PATH / f"{topic_name}/{topic_name}.sdf"
                    viz_file = BASE_PATH / f"{topic_name}/{topic_name}.py"
                    download_sdf_file()
                    mols = load_mols()
                    viz = load_viz()
                    sdf_time = sdf_file.stat().st_mtime
                    # wenn sdf time älter als 1 min ist True
                    if sdf_time < (time.time() - 60):
                        check_sdf = True
                    else:
                        check_sdf = False
                    sdf_file.touch()
                    for subtopic_name, subtopic_t in topic.items():
                        subtopic, subtopic_path = subtopic_t
                        # st.subheader(f"Subtopic: {subtopic_name}")
                        mols, viz = update_single_topic(subtopic, subtopic_path, mols, viz, check_sdf=check_sdf)
                    
                    save_mols()
                    save_viz()

                    Dashboard.Visualizer.visualize_in_3Dmol(mols, viz, state)
                    if st.button(f"Fortschritt für {topic_name} einklappen"):
                        open_topic = ""
                        st.empty()

        logging.info("Checking progress of all jobs.")
        topics_dir = BASE_PATH
        topics = {}
        async def process_topic(topic_name):
            topic_path = Path(topics_dir) / topic_name
            topic_progress = await asyncio.to_thread(check_progress_of_single_topic, topic_path)
            if topic_progress:
                topics[topic_name] = topic_progress

        async def process_all_topics():
            tasks = [process_topic(topic_name) for topic_name in os.listdir(topics_dir)]
            await asyncio.gather(*tasks)

        asyncio.run(process_all_topics())
        update_dashboard(topics)
        logging.info("Finished checking progress of all jobs.")

Dashboard.check_progress_of_all_jobs()
Dashboard.upload_file_and_start_calculation()
