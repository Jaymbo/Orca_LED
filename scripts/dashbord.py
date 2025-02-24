import time

start_time = time.time()
import os
print(f"os imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
import re
print(f"re imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
import streamlit as st
print(f"streamlit imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from pathlib import Path
print(f"pathlib imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
import pipeline
print(f"pipeline imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from LED_extraction import extract_LED_energy
print(f"LED_extraction imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from csv_to_viz import extrakt
print(f"csv_to_viz imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
import subprocess
print(f"subprocess imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from visualization import MoleculeVisualizer
print(f"visualization imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from config_manager import FragmentConfig
print(f"config_manager imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
import logging
print(f"logging imported in {time.time() - start_time:.4f} seconds")

start_time = time.time()
from logging.handlers import RotatingFileHandler
print(f"RotatingFileHandler imported in {time.time() - start_time:.4f} seconds")

from database import Database

topics_progress_b, total_files_b, completed_files_b, pending_jobs_b = {}, 0, 0, []

def setup_logging():
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    log_file = Path(__file__).parent / 'logs' / 'app.log'
    
    log_file.parent.mkdir(exist_ok=True)
    
    logging.basicConfig(
        handlers=[
            RotatingFileHandler(
                log_file,
                maxBytes=1e6,  # 1MB
                backupCount=3
            )
        ],
        level=logging.INFO,
        format=log_format,
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def visualize_xyz(molecule_path: Path):
    viewer = MoleculeVisualizer.render_xyz(molecule_path)
    st.session_state.previous_xyz = viewer  # Nur wenn benötigt
    st.components.v1.html(viewer._make_html(), height=600)

def visualize_in_3Dmol(molecule_path: Path, lines_path: Path):
    MoleculeVisualizer.show_led_analysis(
        molecule_dir=molecule_path,
        lines_path=lines_path
    )

def split_sdf_file_to_mol2(file_path: Path):
    return MoleculeVisualizer.convert_sdf_to_mol2(file_path)

def handle_file_upload(uploaded_file) -> Path:
    """Verarbeitet hochgeladene Dateien und gibt gespeicherten Pfad zurück"""
    try:
        SANITIZE_TABLE = str.maketrans(" ", "_", "!$%&()=+,-/:;<=>?@[\]^`{|}~")
        clean_name = uploaded_file.name.translate(SANITIZE_TABLE)
        name = Path(clean_name).stem
        save_path = Path(f"/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}/{clean_name}")
        
        save_path.parent.mkdir(parents=True, exist_ok=True)
        save_path.write_bytes(uploaded_file.getbuffer())
        logging.info(f"File uploaded: {uploaded_file.name}")
        return save_path
    except Exception as e:
        logging.error(f"Upload failed: {str(e)}")
        raise

def sanitize_saving(content: list[str], name: str):
    paths = []
    SANITIZE_TABLE = str.maketrans(" ", "_", "!$%&()=+,-/:;<=>?@[\]^`{|}~")
    clean_name = name.translate(SANITIZE_TABLE)
    name = Path(clean_name).stem
    for i in range(len(content)):
        save_path = Path(f"/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}_{i+1}/{name}_{i+1}.mol2")
        save_path.parent.mkdir(parents=True, exist_ok=True)
        save_path.write_text(content[i])
        paths.append(save_path)
    return paths

def upload_file_and_start_calculation():
    config = FragmentConfig()
    default_fragment = ""
    file_paths = []
    
    st.title("Neue Berechnung starten")
    uploaded_files = st.file_uploader("Wählen Sie eine Datei aus", type=["xyz", "mol2", "sdf"], accept_multiple_files=True)
    
    if uploaded_files:
        default_fragment = config.get_fragment(uploaded_files[0].name)
        for f in uploaded_files:
            if f.name.endswith(".sdf"):
                path = handle_file_upload(f)
                file_paths.extend(sanitize_saving(split_sdf_file_to_mol2(path), f.name))
                path.unlink()
                path.parent.rmdir()
            else:
                file_paths.append(handle_file_upload(f))
            
            print(file_paths)
        visualize_xyz(file_paths)

    fragment_input = st.text_input("Fragmentierungsangabe", value=default_fragment)
    header_input = st.text_area("Zusätzliche Header-Einstellungen für die .inp-Datei (optional)", """! DLPNO-CCSD(T) def2-svp def2-svp/C DEF2/J RIJCOSX tightSCF normalPNO LED

%mdci DoDIDplot true end

%maxcore 160000

%mdci
  MaxIter 200
end""")
    if st.button("Berechnung starten"):
        if file_paths:
            for i in range(len(file_paths)):
                file = file_paths[i]
                config.save_fragment(file.name, fragment_input)
                st.info("Die Berechnung wurde in auftrag gegeben...")
                paths, base, frag_len, index_compound = pipeline.ORCAInputFileCreator(str(file_paths[i]), fragment_input, header_input).create_inp_files()
                if fragment_input == "-1":
                    ret = pipeline.ShellScriptCreator.sh_script_erstellen(paths, base, frag_len, index_compound, time="2:00:00", mem=20, main=False)
                else:
                    ret = pipeline.ShellScriptCreator.sh_script_erstellen(paths, base, frag_len, index_compound)
                st.info("Die Berechnung wurde vorbereitet.")
                for r in ret:
                    Database.process_candidate(Path(r).parent)
        else:
            st.error("Bitte wählen Sie eine Datei aus.")

def check_orca_termination(content):
    return "****ORCA TERMINATED NORMALLY****" in content

def get_progress_of_job(output, path):
    runtime = None
    error_status = get_error_file(path)
    if error_status is not None:
        return error_status, None
    
    status, runtime = check_slurm_job_status_and_duration(path)
    if status != "Not Finished":
        return status, runtime

    if check_orca_termination(output):
        return "Progress: 100%", runtime
    count = 0
    output = output.split("\n")
    if "INITIAL GUESS DONE" in output:
        count += 1
    for line in output:
        if "timings" in line.lower():
            count += 1
    if "LED" in output:
        return f"Progress: {int(count / 7 * 100)}%", runtime
    else:
        return f"Progress: {int(count / 4 * 100)}%", runtime

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
            return f"Error: {str(e)}"
    return None

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

@st.fragment(run_every="60s")
def check_progress_of_all_jobs():
    topics_dir = "/lustre/work/ws/ws1/tu_zxofv28-my_workspace/"
    global total_files_b, completed_files_b, pending_jobs_b, topics_progress_b
    topics_progress = {}
    total_files = 0
    completed_files = 0
    pending_jobs = []

    for topic_name in os.listdir(topics_dir):
        topic_path = Path(topics_dir) / topic_name
        if os.path.isdir(topic_path):
            subfolders = sorted(
                [f for f in topic_path.iterdir() if f.is_dir()],
                key=lambda x: ("fragment_" in x.name, x.name)
            )

            check = False
            for folder in subfolders:
                if folder.name.startswith("fragment_"):
                    check = True
                    break

            if check:
                jobs_progress = {}
                for subfolder_name in subfolders:
                    job_path = Path(os.path.join(topic_path, subfolder_name))
                    if os.path.isdir(job_path):
                        output_file = os.path.join(subfolder_name, f"{subfolder_name.stem}.out")
                        if os.path.exists(output_file):
                            total_files += 1
                            with open(output_file, 'r') as f:
                                lines = f.read()

                            progress, runtime = get_progress_of_job(lines, Path(os.path.join(job_path, f"{subfolder_name.stem}")))
                            jobs_progress[subfolder_name.stem] = (progress, runtime, job_path.stem, energy_extraction(lines))

                            if progress != "Progress: 100%":
                                pending_jobs.append(f"{topic_name} - {subfolder_name.stem}")
                            else:
                                completed_files += 1
                        else:
                            jobs_progress[subfolder_name.stem] = ("Not Started", None, job_path.stem, None)

                if jobs_progress:
                    topics_progress[topic_name] = (jobs_progress, topic_path)
    topics_progress = dict(sorted(topics_progress.items(), key=lambda x: x[0]))
    if topics_progress_b != topics_progress or pending_jobs_b != pending_jobs or total_files_b != total_files or completed_files_b != completed_files:
        topics_progress_b = topics_progress
        pending_jobs_b = pending_jobs
        total_files_b = total_files
        completed_files_b = completed_files
    update_dashboard(topics_progress, total_files, completed_files, pending_jobs)

def energy_extraction(context):
    matches = list(re.finditer(r"FINAL SINGLE POINT ENERGY \s+([-+]?\d+\.\d+)", context))
    if not matches:
        return None
    return float(matches[-1].group(1))

def update_dashboard(topics_progress, total_files, completed_files, pending_jobs):

    st.title('ORCA-Status')

    st.write(f"Gesamtfortschritt: {completed_files} von {total_files} abgeschlossen.")
    overall_progress = completed_files / total_files if total_files > 0 else 0
    st.progress(overall_progress)

    if pending_jobs:
        st.subheader("Noch nicht abgeschlossene Jobs:")
        for job in pending_jobs:
            st.text(job)
    else:
        st.text("Alle Jobs sind abgeschlossen.")

    for topic, (jobs, folder) in topics_progress.items():
        st.subheader(f"Topic: {topic}")

        total_jobs = len(jobs)
        completed_jobs = sum(1 for progress, _, _, _ in jobs.values() if "Progress: 100%" == progress)
        not_started_jobs = sum(1 for progress, _, _, _ in jobs.values() if progress == "Not Started")

        if completed_jobs == total_jobs:
            st.markdown(f"<div style='background-color: green; height: 10px; width: 100%'></div>", unsafe_allow_html=True)
            energy = 2 * jobs[list(jobs.keys())[0]][3]
            for i in range(completed_jobs):
                energy -= jobs[list(jobs.keys())[i]][3]
            st.text(f"{energy * 627.509474:.2f} kcal/mol")
            extract_LED_energy(folder)
            extrakt(folder)
            visualize_in_3Dmol(Path(folder), Path(folder) / "viz.py")

            st.text(f"Alle Berechnungen abgeschlossen!")
        elif not_started_jobs == total_jobs:
            st.markdown(f"<div style='background-color: grey; height: 10px; width: 100%'></div>", unsafe_allow_html=True)
            st.text(f"Alle Berechnungen noch nicht gestartet.")
        else:
            num_columns = 7
            cols = st.columns(num_columns)

            for i, (subfolder_name, (progress, runtime, path, energy)) in enumerate(jobs.items()):
                subfolder_name = subfolder_name.replace("fragment_", "")
                color, progress_value = get_color_and_progress(progress)
                if len(subfolder_name) > 6:
                    subfolder_name = subfolder_name[:6] + "..."

                col = cols[i % num_columns]
                with col:
                    if color == "red" or color == "grey":
                        st.text(f"{subfolder_name}:\n{progress}\nNo Time")
                    else:
                        if runtime:
                            st.text(f"{subfolder_name}:\n{progress_value}%\n{runtime}")
                        else:
                            st.text(f"{subfolder_name}:\n{progress_value}%\nNo Time")
                    st.markdown(f"<div style='background-color: {color}; height: 10px; width: {progress_value}%'></div>", unsafe_allow_html=True)

            st.text(f"Fortschritt des Topics: {completed_jobs}/{total_jobs} abgeschlossen")
    
setup_logging()
check_progress_of_all_jobs()
upload_file_and_start_calculation()
