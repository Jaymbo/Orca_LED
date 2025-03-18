from rdkit import Chem
import numpy as np
import openpyxl
from pathlib import Path

class SdfXyzMerger:
    def __init__(self, sdf_file, xyz_file, folder):
        self.sdf_file = sdf_file
        self.xyz_file = xyz_file
        self.xlsx_file = folder + "/Summary_fp-LED_matrices.xlsx"
        self.xlsx_file2 = folder + "/Summary_Standard_LED_matrices.xlsx"
        self.sdf_molecules = []
        self.coordinates = []
        self.xlsx_data = None
    
    def load_sdf(self):
        suppl = Chem.SDMolSupplier(self.sdf_file, removeHs=False)
        self.sdf_molecules = [mol for mol in suppl if mol is not None]
    
    def load_xyz(self):
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
            atom_count = int(lines[0].strip())
            self.coordinates = np.array([list(map(float, line.split()[1:])) for line in lines[2:2+atom_count]])
            self.coordinates = self.coordinates[self.coordinates[:, 0].argsort()]
    
    def load_xlsx(self):
        workbook = openpyxl.load_workbook(self.xlsx_file, data_only=True)
        sheet = workbook.active
        self.xlsx_data = [[cell.value for cell in row] for row in sheet.iter_rows()]
    
    def match_and_update(self):
        for mol in self.sdf_molecules:
            conf = mol.GetConformer()
            sdf_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            sdf_coords = sdf_coords[sdf_coords[:, 0].argsort()]
            energies = {}
            if self.compare_coordinates(sdf_coords, self.coordinates):
                for i in range(1,len(self.xlsx_data)):
                    row = self.xlsx_data[i]
                    for j in range(i,len(row)):
                        value = row[j]
                        if value is not None:
                            energies[f"{i}-{j}"] = value
                mol.SetProp(Path(self.xlsx_file).name, str(energies))
                break

    def compare_coordinates(self, sdf_coords, xyz_coords, tolerance=0.1):
        if sdf_coords.shape != xyz_coords.shape:
            print("Shapes do not match.")
            return False
        return np.allclose(sdf_coords, xyz_coords, atol=tolerance)
    
    def save_updated_sdf(self):
        writer = Chem.SDWriter(self.sdf_file)
        for mol in self.sdf_molecules:
            writer.write(mol)
        writer.close()
    
    def run(self):
        self.load_sdf()
        self.load_xyz()
        self.load_xlsx()
        self.match_and_update()
        self.xlsx_file = self.xlsx_file2
        self.load_xlsx()
        self.match_and_update()
        self.save_updated_sdf()
