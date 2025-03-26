import logging
from pathlib import Path
import datetime
from shutil import copy2
import numpy as np
from scipy.optimize import linear_sum_assignment
import time
import os

BASE_PATH = Path(__file__).resolve().parent.parent

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Database:
    def __init__(self, dir: Path, file_cache) -> None:
        self.base = BASE_PATH / "database/"
        self.base.mkdir(parents=True, exist_ok=True)
        self.dir: Path = dir
        self.get_file_paths()
        self.rmsd_value = 0.001
        self.ends = ["_out.out", ".out", ".densities", "_err.err", ".property.txt", ".bibtex", ".cube", ".densitiesinfo", ".sh", ".inp"]
        self.file_cache = file_cache

    def get_file_paths(self):
        self.filename = self.dir.parent / f"{self.dir.stem}.xyz"
        self.header_filename = self.dir / f"{self.dir.stem}.inp"

    def get_filecontent(self) -> str:
        return self.file_cache.get(self.filename)

    def atoms_from_filecontent(self, content: str) -> list:
        lines = [line for line in content.splitlines() if line.strip()]
        return [line.split()[0] for line in lines[2:]]

    def header_from_file(self) -> str:
        content = self.file_cache.get(self.header_filename)
        return "".join(line.strip() for line in content.splitlines())

    def atoms_to_str(self, atoms: list) -> str:
        sorted_atoms = sorted(atom.upper() for atom in atoms)
        unique, counts = np.unique(sorted_atoms, return_counts=True)
        return "".join(f"{a}{c}" for a, c in zip(unique, counts))

    def header_to_str(self, header: str) -> str:
        header_str = "".join(char for char in header.upper() if char.isalnum())
        return header_str[:55]

    def date_to_str(self) -> str:
        return datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")[:-4]

    def create_filename(self, atoms: list) -> tuple[Path, Path]:
        atoms_str = self.atoms_to_str(atoms)
        date_str = self.date_to_str()
        name = f"{atoms_str}_{date_str}"
        dirpath = self.base / name
        filepath = dirpath / name
        return dirpath, filepath

    def get_database_names(self) -> list:
        return [f.name for f in self.base.iterdir() if f.is_dir()]

    def find_matches(self, current_name: str, database_names: list) -> list:
        name_parts = current_name.split("_")[:1]
        return [name for name in database_names if name.split("_")[:1] == name_parts]

    def compute_distance_matrix(self, coords: np.ndarray) -> np.ndarray:
        diff = coords[:, None, :] - coords[None, :, :]
        return np.sqrt(np.sum(diff**2, axis=-1))

    def match_distance_matrices(self, D1: np.ndarray, D2: np.ndarray):
        cost_matrix = np.sum((D1[:, None, :] - D2[None, :, :])**2, axis=-1)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        return row_ind, col_ind, cost_matrix

    def rmsd_of_distance_matrices(self, D1: np.ndarray, D2: np.ndarray, row_ind, col_ind) -> float:
        D2_perm = D2[np.ix_(row_ind, col_ind)]
        diff = D1 - D2_perm
        return np.sqrt(np.sum(diff**2) / D1.size)

    def rmsd(self, xyz1: str, xyz2: str) -> float:
        lines1 = [line.strip() for line in xyz1.splitlines() if line.strip()]
        data1 = lines1[2:]
        atoms1 = [line.split()[0] for line in data1]
        coords1 = np.array([[float(x) for x in line.split()[1:4]] for line in data1])

        lines2 = [line.strip() for line in xyz2.splitlines() if line.strip()]
        data2 = lines2[2:]
        atoms2 = [line.split()[0] for line in data2]
        coords2 = np.array([[float(x) for x in line.split()[1:4]] for line in data2])

        if sorted(atoms1) != sorted(atoms2):
            return 100.0, None

        D1 = self.compute_distance_matrix(coords1)
        D2 = self.compute_distance_matrix(coords2)
        row_ind, col_ind, _ = self.match_distance_matrices(D1, D2)
        return self.rmsd_of_distance_matrices(D1, D2, row_ind, col_ind), col_ind

    def right_fragmentation(self, header_str, col_ind):
        fragmentation = self.get_fragmentation(header_str)
        if col_ind is not None:
            fragmentation = [np.sort(np.array([col_ind[int(f)] for f in frag], dtype=int)) for frag in fragmentation]
        return sorted(fragmentation, key=lambda x: x.tolist())

    def get_fragmentation(self, header_str):
        fragments = []
        lines = header_str.splitlines()
        frag_start = np.where(np.char.find(lines, "Fragments") >= 0)[0]
        frag_end = np.where(np.char.find(lines, "end") >= 0)[0]

        if frag_start.size > 0 and frag_end.size > 0:
            frag_lines = lines[frag_start[0] + 1:frag_end[0]]
            fragments = [np.sort(line.split("{")[1].split("}")[0].split(), axis=0) for line in frag_lines]
        return sorted(fragments, key=lambda x: x.tolist())

    def molecule_exists(self, candidate_xyz: str, matched: list, header_str: str) -> tuple[list, bool]:
        """Check if the molecule already exists in the database"""
        if not matched:
            return matched, False
        fragments = self.right_fragmentation(header_str, None)
        header_header = "".join(char for char in header_str.split("*XYZfile")[0].strip() if char.isalnum())
        def check_header(folder):
            xyz_path = self.base / folder / f"{folder}.inp"
            content = "".join(line.strip() for line in self.file_cache.get(xyz_path).splitlines())
            output = "".join(char for char in content.split("*XYZfile")[0].strip() if char.isalnum())
            return output == header_header

        def compute_rmsd_and_col(folder):
            xyz_path = self.base / folder / f"{folder}.xyz"
            xyz_content = self.file_cache.get(xyz_path)
            return self.rmsd(candidate_xyz, xyz_content)
        matched = [folder for folder in matched if check_header(folder)]
        if not matched:
            return matched, False
        results = [compute_rmsd_and_col(folder) for folder in matched]
        rmsd_list = np.array([result[0] for result in results], dtype=float)
        col_ind = [result[1] for result in results]

        def check_fragmentation(folder, col):
            inp_path = self.base / folder / f"{folder}.inp"
            inp_content = self.file_cache.get(inp_path)
            return self.right_fragmentation(inp_content, col) == fragments

        fragmentation_matches = np.array([check_fragmentation(folder, col) for folder, col in zip(matched, col_ind)])

        new_matched = np.array(matched)[fragmentation_matches]
        new_rmsd_list = rmsd_list[fragmentation_matches]

        sorted_indices = np.argsort(new_rmsd_list)
        matched = new_matched[sorted_indices].tolist()
        rmsd_list = new_rmsd_list[sorted_indices].tolist()

        exists = bool(rmsd_list) and min(rmsd_list) < self.rmsd_value
        return matched, exists

    def insert(self, dirpath: Path, filepath: Path) -> None:
        dirpath.mkdir(parents=True, exist_ok=True)
        filepath = Path(str(filepath.parent) + "/" + str(filepath.stem))
        xyz_path = filepath.with_suffix(".xyz")
        for end in self.ends:
            out_path = Path(str(filepath) + end)
            self.file_cache.set(out_path, "test")
        self.file_cache.set(xyz_path, self.file_cache.get(self.filename))
        self.create_symlink(xyz_path)

    def create_symlink(self, out_path: Path) -> None:
        symlink_base = self.header_filename.parent / self.header_filename.stem
        for end in self.ends:
            symlink = Path(f"{symlink_base}{end}")
            self.file_cache.link(out_path, symlink)

    def file_end(self, filename: Path) -> str:
        for match in self.ends:
            if filename.name.endswith(match):
                return match
        return

    @classmethod
    def process_candidate(cls, dir: Path, file_cache) -> Path:
        db = cls(dir, file_cache)
        candidate_xyz = db.get_filecontent()
        atoms = db.atoms_from_filecontent(candidate_xyz)
        header = db.header_from_file()
        new_dir, new_filepath = db.create_filename(atoms)
        
        database_names = db.get_database_names()
        matched = db.find_matches(new_filepath.name, database_names)
        matched, exists = db.molecule_exists(candidate_xyz, matched, header)    

        if not exists:
            db.insert(new_dir, new_filepath)
            return str(new_filepath) + ".xyz"
        else:
            existing_folder = db.base / matched[0]
            out_path = existing_folder / (matched[0] + ".out")
            db.create_symlink(out_path)
            return None

    def add_calculation(self):
        candidate_xyz = self.get_filecontent()
        atoms = self.atoms_from_filecontent(candidate_xyz)
        new_dir, new_filepath = self.create_filename(atoms)
        new_dir.mkdir(parents=True, exist_ok=True)

        for file in self.dir.iterdir():
            end = self.file_end(file)
            if end is not None:
                destination = new_dir / new_filepath.stem
                destination = Path(f"{destination}{end}")
                copy2(file, destination)

        destination = new_dir / new_filepath.stem
        destination = destination.with_suffix(".xyz")
        copy2(self.filename, destination)

        destination = new_dir / new_filepath.stem
        destination = destination.with_suffix(".inp")
        copy2(self.header_filename, destination)

    def cleanup(self):
        database_names = self.get_database_names()
        for folder in self.base.iterdir():
            try:
                matched = self.find_matches(folder.name, database_names)
                if len(matched) > 1:
                    xyz_files = [self.base / m / f"{m}.xyz" for m in matched if (self.base / m / f"{m}.xyz").exists()]
                    if not (xyz_files[0].exists() and xyz_files[1].exists()):
                        continue
                    content0 = self.file_cache.get(xyz_files[0])
                    content1 = self.file_cache.get(xyz_files[1])
                    current_rmsd = self.rmsd(content0, content1)
                    if current_rmsd < self.rmsd_value:
                        self.syslink_merge(xyz_files[0], xyz_files[1])
                        rm = self.base / matched[0]
                        for file in rm.iterdir():
                            file.unlink()
                        rm.rmdir()
            except Exception as e:
                logging.error("Error in folder %s: %s", folder.name, e)

    def syslink_merge(self, original: Path, merged: Path):
        for folder in BASE_PATH.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        for file in subfolder.iterdir():
                            if file.is_symlink() and file.resolve() == original:
                                file.unlink()
                                file.symlink_to(merged)

    def copy_to_db(self):
        for folder in BASE_PATH.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        Database(subfolder).add_calculation()
        self.cleanup()
