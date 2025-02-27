import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from LEDAW.ledaw_package.nbody_engine import engine_LED_N_body

def check_orca_termination(content):
    """
    Überprüft, ob ein ORCA-Lauf korrekt terminiert wurde.
    :param file_path: Pfad zur ORCA-Ausgabedatei
    :return: True, wenn ORCA normal terminiert wurde, sonst False
    """
    return "****ORCA TERMINATED NORMALLY****" in content

def extract_LED_energy(base):
    print(f"Extracting LED energies from {base}")
    dateinamen = [folder for folder in os.listdir(base) if os.path.isdir(os.path.join(base, folder))]
    dateinamen.sort(key=lambda x: (x.startswith('fragment_'), x))
    name = []
    content = []
    filepaths = []
    for folder in dateinamen:
        try:
            filepath = f"{base}/{folder}/{folder}.out"
            name.append(filepath)
            filepaths.append(filepath)
            with open(filepath, "r") as file:
                content.append(file.read())
            print(f"Successfully read file: {filepath}")
        except FileNotFoundError:
            print(f"Datei {filepath} nicht gefunden.")
            continue

    if " LED " not in content[0]:
        print("LED keyword not found in the first file content.")
        return
    
    LEDAW_output_path = f"./{base.stem}"
    print(f"LEDAW_output_path: {LEDAW_output_path}")

    terminated = []
    for file_content in content:
        terminated.append(check_orca_termination(file_content))

    if not all(terminated):
        print("Not all files terminated normally.")
        return
    method = "DLPNO-CCSD(T)"
    ### Convert Hartree to kcal/mol
    conversion_factor = 627.5095  # for kj/mol, use: 2625.5
    alternative_filenames = ['' for i in range(len(dateinamen))]
    try:
        engine_LED_N_body(main_filenames=filepaths, 
                    alternative_filenames=alternative_filenames, 
                    conversion_factor=conversion_factor, 
                    method=method,
                    LEDAW_output_path=LEDAW_output_path)
        print("LED engine execution completed successfully.")
    except Exception as e:
        print(f"Error in LED engine: {e}")