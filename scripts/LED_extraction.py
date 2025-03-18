import os
import sys
import concurrent.futures
import logging
import time
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from LEDAW.ledaw_package.nbody_engine import engine_LED_N_body

def check_orca_termination(content):
    return "****ORCA TERMINATED NORMALLY****" in content


def read_file(filepath):
    try:
        with open(filepath, "r") as file:
            content = file.read()
        print(f"Successfully read file: {filepath}")
        return filepath, content
    except FileNotFoundError:
        print(f"Datei {filepath} nicht gefunden.")
        return filepath, None

def extract_LED_energy(base):
    xlsx_file = f"{base}/Summary_fp-LED_matrices.xlsx"
    xyz_file = f"{base}/{base.stem}.xyz"
    if not os.path.exists(xyz_file):
        print(f"File {xyz_file} does not exist.")
        return 
    
    if os.path.exists(xlsx_file):
        print(f"File {xlsx_file} exists.")
        return
    
    print(f"Extracting LED energies from {base}")
    base_path = Path(base)
    dateinamen = [folder for folder in base_path.iterdir() if folder.is_dir()]
    dateinamen.sort(key=lambda x: (x.name.startswith('fragment_'), x.name))
    dateinamen.sort(key=lambda x: (x.name.startswith('subsys_'), x.name))
    
    filepaths = [str(folder / f"{folder.name}.out") for folder in dateinamen]
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(read_file, filepaths))
    
    name, content = zip(*results)
    content = [c for c in content if c is not None]
    print(f"len(dateinamen): {len(dateinamen)}"	)
    if not content or " LED " not in content[0] or len(dateinamen) <= 2:
        print("LED keyword not found in the first file content.")
        return
    
    LEDAW_output_path = f"{base_path}"
    print(f"LEDAW output path: {LEDAW_output_path}")

    terminated = [check_orca_termination(file_content) for file_content in content]

    if not all(terminated):
        print("Not all files terminated normally.")
        return
    
    method = "DLPNO-CCSD(T)"
    conversion_factor = 627.5096080305927  # for kj/mol, use: 2625.5
    alternative_filenames = ['' for _ in range(len(dateinamen))]
    relabel_mapping = [i+1 for i in range(len(dateinamen)-1)]
    use_ref_as_rhf_in_hfld = None
    engine_LED_N_body(main_filenames=filepaths, 
                        alternative_filenames=alternative_filenames, 
                        conversion_factor=conversion_factor,
                        method=method,
                        LEDAW_output_path=LEDAW_output_path)
    print("LED engine execution completed successfully.")