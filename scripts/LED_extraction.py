import os
import sys
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from LEDAW.ledaw_package.nbody_engine import engine_LED_N_body

class LEDExtractor:
    def __init__(self, base):
        self.base = base
        self.base_path = Path(base)
        self.xlsx_file = f"{base}/Summary_fp-LED_matrices.xlsx"
        self.xyz_file = f"{base}/{self.base_path.stem}.xyz"
        self.filepaths = []
        self.content = []

    def check_orca_termination(self, content):
        return "****ORCA TERMINATED NORMALLY****" in content

    def read_file(self, filepath):
        try:
            with open(filepath, "r") as file:
                content = file.read()
            print(f"Successfully read file: {filepath}")
            return filepath, content
        except FileNotFoundError:
            print(f"Datei {filepath} nicht gefunden.")
            return filepath, None

    def extract_LED_energy(self):
        if not os.path.exists(self.xyz_file):
            return 
        
        if os.path.exists(self.xlsx_file):
            return
        
        print(f"Extracting LED energies from {self.base}")
        dateinamen = [folder for folder in self.base_path.iterdir() if folder.is_dir()]
        dateinamen.sort(key=lambda x: (x.name.startswith('fragment_'), x.name))
        dateinamen.sort(key=lambda x: (x.name.startswith('subsys_'), x.name))
        
        self.filepaths = [str(folder / f"{folder.name}.out") for folder in dateinamen]
        
        for filepath in self.filepaths:
            _, content = self.read_file(filepath)
            if content is not None:
                self.content.append(content)
        
        print(f"len(dateinamen): {len(dateinamen)}")
        if not self.content or " LED " not in self.content[0] or len(dateinamen) <= 2:
            print("LED keyword not found in the first file content.")
            return
        
        LEDAW_output_path = f"{self.base_path}"
        print(f"LEDAW output path: {LEDAW_output_path}")

        terminated = [self.check_orca_termination(file_content) for file_content in self.content]

        if not all(terminated):
            print("Not all files terminated normally.")
            return
        
        method = "DLPNO-CCSD(T)"
        conversion_factor = 627.5096080305927  # for kj/mol, use: 2625.5
        alternative_filenames = ['' for _ in range(len(dateinamen))]
        engine_LED_N_body(main_filenames=self.filepaths, 
                          alternative_filenames=alternative_filenames, 
                          conversion_factor=conversion_factor,
                          method=method,
                          LEDAW_output_path=LEDAW_output_path)
        print("LED engine execution completed successfully.")