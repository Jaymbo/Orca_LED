from rdkit import Chem
import numpy as np
import openpyxl
from pathlib import Path
from XBPy.module.xbpy.rdutil.io import read_molecules
import time
import os

class SdfXyzMerger:
    def __init__(self, sdf_file, xyz_file, folder):
        self.sdf_file = sdf_file
        # wenns vor mehr als 3 min bearbeitet wurde
        if not self.sdf_file.exists():
            self.make_sdf()
        self.xlsx_file = folder / "Summary_fp-LED_matrices.xlsx"
        self.xlsx_file2 = folder / "Summary_Standard_LED_matrices.xlsx"
        if self.xlsx_file == None or self.xlsx_file2 == None or not self.xlsx_file.exists() or not self.xlsx_file2.exists():
            print("Excel file not found.")
            self.xlsx_file = None
            self.xlsx_file2 = None
        else:
            self.xyz_file = xyz_file
            if self.xyz_file.stat().st_mtime < time.time() - 3 * 60:
                return
            self.sdf_molecules = []
            self.coordinates = []
            self.xlsx_data = None
            self.run()
    
    def load_sdf(self):
        if not self.sdf_file.exists() or self.sdf_file.stat().st_size == 0:
            self.sdf_molecules = []
        else:
            suppl = Chem.SDMolSupplier(self.sdf_file, removeHs=False)
            self.sdf_molecules = [mol for mol in suppl if mol is not None]
    
    def load_xyz(self):
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()
            atom_count = int(lines[0].strip())
            self.coordinates = np.array([list(map(float, line.split()[1:])) for line in lines[2:2+atom_count]])
            self.coordinates = self.coordinates[self.coordinates[:, 0].argsort()]
    
    def load_xlsx(self, xlsx_file=None):
        if self.xlsx_file == None or not xlsx_file.exists():
            return
        workbook = openpyxl.load_workbook(xlsx_file, data_only=True)
        sheet = workbook.active
        self.xlsx_data = [[cell.value for cell in row] for row in sheet.iter_rows()]
    
    def match_and_update(self):
        found = False
        mol = None
        for i in range(len(self.sdf_molecules)):
            mol = self.sdf_molecules[i]
            conf = mol.GetConformer()
            sdf_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            sdf_coords = sdf_coords[sdf_coords[:, 0].argsort()]
            if self.compare_coordinates(sdf_coords, self.coordinates):
                print(f"Found matching coordinates in SDF file for molecule {i+1}.")
                found = True
                break
        if not found:
            if self.append_xyz_to_sdf():
                mol =  self.sdf_molecules[-1]
        if mol is None:
            print("No matching molecule found.")
        else:
            self.write_properties(mol)
            self.xlsx_file = self.xlsx_file2
            self.load_xlsx(self.xlsx_file)
            self.write_properties(mol)
        
    def write_properties(self, mol):
        energies = {}
        for i in range(1,len(self.xlsx_data)):
            row = self.xlsx_data[i]
            for j in range(i,len(row)):
                value = row[j]
                if value is not None:
                    energies[f"{i}-{j}"] = value
        mol.SetProp(Path(self.xlsx_file).name, str(energies))

        
    def compare_coordinates(self, sdf_coords, xyz_coords, tolerance=0.1):
        if sdf_coords.shape != xyz_coords.shape:
            return False
        return np.allclose(sdf_coords, xyz_coords, atol=tolerance)
    
    def save_updated_sdf(self):
        writer = Chem.SDWriter(self.sdf_file)
        for mol in self.sdf_molecules:
            writer.write(mol)
        writer.close()
    
    def make_sdf(self):
        """
        Create an empty sdf file with the name of the sdf file
        """
        if not self.sdf_file.exists():
            self.sdf_file.touch()
        
    def append_xyz_to_sdf(self):
        """
        Append the coordinates from the xyz file to the sdf file
        """
        try:
            mol = read_molecules(str(self.xyz_file))
            mol = next(mol)
            if self.sdf_file.exists():
                with open(self.sdf_file, 'r') as f:
                    lines = f.readlines()
                if len(lines) == 0:
                    supplier = []
                else:
                    supplier = Chem.SDMolSupplier(self.sdf_file, removeHs=False)
            else:
                supplier = []
            mols = [mol for mol in supplier]
            mols.append(mol)
            writer = Chem.SDWriter(self.sdf_file)
            for mol in mols:
                writer.write(mol)
            writer.close()
            os.utime(self.xyz_file, (time.time(), time.time()))
            self.load_sdf()
            print(f"Appended XYZ coordinates to SDF file: {self.sdf_file}")
            return True
        except Exception:
            return False


    def run(self):
        self.load_sdf()
        self.load_xyz()
        self.load_xlsx(self.xlsx_file)
        if self.xlsx_data is None:
            return
        self.match_and_update()
        self.save_updated_sdf()
