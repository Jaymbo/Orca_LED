import numpy as np
import openpyxl
from pathlib import Path
from xbpy.rdutil.io import read_molecules

class SdfXyzMerger:
    def __init__(self, xyz_file, folder, mols=[], viz = ""):
        self.xlsx_file = folder / "Summary_fp-LED_matrices.xlsx"
        self.xlsx_file2 = folder / "Summary_Standard_LED_matrices.xlsx"
        self.mols = mols
        self.viz = viz
        self.xyz_file = xyz_file
        self.viz_file = folder / "viz.py"
        self.coordinates = []
        self.xlsx_data = None
    
    def load_xyz(self):
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
            atom_count = int(lines[0].strip())
            self.coordinates = np.array([list(map(float, line.split()[1:])) for line in lines[2:2+atom_count]])
    
    def load_xlsx(self, xlsx_file=None):
        if not xlsx_file or not xlsx_file.exists():
            return
        workbook = openpyxl.load_workbook(xlsx_file, data_only=True)
        sheet = workbook.active
        self.xlsx_data = [[cell.value for cell in row] for row in sheet.iter_rows()]
    
    def match_and_update(self):
        found = False
        mol = None
        for i, mol in enumerate(self.mols):
            conf = mol.GetConformer()
            sdf_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            if self.compare_coordinates(sdf_coords, self.coordinates):
                print(f"Found matching coordinates in SDF file for molecule {i+1}.")
                found = True
                break
        if not found:
            worked = self.append_xyz_to_sdf()
            if worked:
                mol = self.mols[-1]
            else:
                return False
        if mol is None:
            print("No matching molecule found.")
            return False
        elif not self.properties_exist(mol):
            self.write_properties(mol)
            self.xlsx_file = self.xlsx_file2
            self.load_xlsx(self.xlsx_file)
            self.write_properties(mol)
            return True
        else:
            return False
    def properties_exist(self, mol):
        try:
            prop = mol.GetProp(Path(self.xlsx_file).name)
            if prop:
                return True
        except Exception:
            return False
        
    def write_properties(self, mol):
        energies = {}
        for i, row in enumerate(self.xlsx_data[1:], start=1):
            for j, value in enumerate(row[i:], start=i):
                if value is not None:
                    energies[f"{i}-{j}"] = value
        mol.SetProp(Path(self.xlsx_file).name, str(energies))
        
    def compare_coordinates(self, sdf_coords, xyz_coords, tolerance=0.1):
        if sdf_coords.shape != xyz_coords.shape:
            return False
        return np.allclose(sdf_coords, xyz_coords, atol=tolerance)
        
    def append_xyz_to_sdf(self):
        try:
            mol = next(read_molecules(str(self.xyz_file)))
            self.mols.append(mol)
            return True
        except Exception:
            return False

    def add_viz(self):
        with open(self.viz_file, 'r') as f:
            lines = f.readlines()
        if self.viz == "":
            self.viz = ''.join(lines)
        else:
            zs_viz = []
            counter = 1
            for line in self.viz.split("\n"):
                if line.startswith("Line"):
                    counter += 1
                if line.startswith("for"):
                    break
                zs_viz.append(line + "\n")
            start = False
            for line in lines:
                if not start and line.startswith("Line"):
                    start = True
                    zs_viz.append(line)
                    continue
                elif start and line.startswith("cmd"):
                    zs_viz.append(f"""cmd.load_cgo(Lines, "Lines", state={counter})""")
                elif start:
                    zs_viz.append(line)
            self.viz = ''.join(zs_viz)

    def run(self):
        if not self.xlsx_file.exists() or not self.xlsx_file2.exists():
            return self.mols, self.viz
        self.load_xyz()
        self.load_xlsx(self.xlsx_file)
        if self.xlsx_data is None:
            return self.mols, self.viz
        worked = self.match_and_update()
        if worked:
            self.add_viz()
        return self.mols, self.viz
