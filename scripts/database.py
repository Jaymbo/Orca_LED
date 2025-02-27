import logging
from pathlib import Path
import datetime
from shutil import copy2
import numpy as np
from scipy.optimize import linear_sum_assignment

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Database:
    def __init__(self, dir: Path) -> None:
        self.base = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/database/")
        self.base.mkdir(parents=True, exist_ok=True)
        self.dir: Path = dir
        self.get_file_paths()
        self.rmsd_value = 0.001
        self.ends = ["_out.out", ".out", ".densities", "_err.err", ".property.txt", ".bibtex", ".cube", ".densitiesinfo", ".sh", ".inp"]

    def get_file_paths(self):
        self.filename = self.dir.parent / f"{self.dir.stem}.xyz"
        self.header_filename = self.dir / f"{self.dir.stem}.inp"

    def get_filecontent(self) -> str:
        with self.filename.open("r") as f:
            return f.read()

    def atoms_from_filecontent(self, content: str) -> list:
        lines = [line for line in content.splitlines() if line.strip()]
        return [line.split()[0] for line in lines[2:]]

    def header_from_file(self) -> str:
        with self.header_filename.open("r") as f:
            return "".join(line.strip() for line in f)

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
        # logging.info("Row indices: %s, Column indices: %s", row_ind, col_ind)
        return self.rmsd_of_distance_matrices(D1, D2, row_ind, col_ind), col_ind

    def right_fragmentation(self, header_str, col_ind):
        fragmentation = self.get_fragmentation(header_str)
        for frag in fragmentation:
            for i, f in enumerate(frag):
                if col_ind is not None:
                    frag[i] = int(f)
                else:
                    frag[i] = col_ind[int(f)]
            frag.sort()
        fragmentation.sort()
        return fragmentation

    def get_fragmentation(self, header_str):
        header_lines = header_str.splitlines()
        frag = False
        fragments = []
        for line in header_lines:
            if "end" in line:
                frag = False
            if frag:
                print(line)
                fragments.append(line.split("{")[1].split("}")[0].split(" ").sort())
            if "Fragments" in line:
                frag = True
        fragments.sort()
        return fragments

    def molecule_exists(self, candidate_xyz: str, matched: list, header_str: str) -> tuple[list, bool]:
        """Check if the molecule already exists in the database"""
        rmsd_list = []
        fragments = self.right_fragmentation(header_str, None)
        col_ind = []

        for folder in matched:
            xyz_path = self.base / folder / f"{folder}.xyz"
            with xyz_path.open("r") as f:
                content = f.read()
            erg, col = self.rmsd(candidate_xyz, content)
            col_ind.append(col)
            rmsd_list.append(erg)

        new_matched = []
        new_erg = []
        for i, folder in enumerate(matched):
            inp_path = self.base / folder / f"{folder}.inp"
            with inp_path.open("r") as f:
                content = f.read()
            if self.right_fragmentation(content, col_ind[i]) == fragments:
                new_matched.append(folder)
                new_erg.append(rmsd_list[i])

        rmsd_list = new_erg
        matched = new_matched
        logging.info("Calculated RMSD values: %s", rmsd_list)
        matched = [x for _, x in sorted(zip(rmsd_list, matched))]
        exists = bool(rmsd_list) and min(rmsd_list) < self.rmsd_value
        return matched, exists

    def insert(self, dirpath: Path, filepath: Path) -> None:
        dirpath.mkdir(parents=True, exist_ok=True)
        filepath = Path(str(filepath.parent) + "/" + str(filepath.stem))
        xyz_path = filepath.with_suffix(".xyz")
        for end in self.ends:
            out_path = Path(str(filepath) + end)
            with out_path.open("w") as f:
                f.write("test")
        self.header_filename.unlink()
        copy2(self.filename, xyz_path)
        self.create_symlink(xyz_path)

    def create_symlink(self, out_path: Path) -> None:
        for file in out_path.parent.iterdir():
            end = self.file_end(file)
            if end is not None:
                symlink = Path(f"{self.header_filename.parent}/{self.header_filename.stem}{end}")
                if symlink.exists():
                    symlink.unlink()
                symlink.symlink_to(file)

    def file_end(self, filename: Path) -> str:
        for match in self.ends:
            if filename.name.endswith(match):
                return match
        return

    @classmethod
    def process_candidate(cls, dir: Path) -> Path:
        db = cls(dir)
        candidate_xyz = db.get_filecontent()
        atoms = db.atoms_from_filecontent(candidate_xyz)
        header = db.header_from_file()
        new_dir, new_filepath = db.create_filename(atoms)
        database_names = db.get_database_names()
        matched = db.find_matches(new_filepath.name, database_names)
        logging.info("Matched folders: %s", matched)
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
        header = self.header_from_file()
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

        logging.info("Calculation files copied to database folder: %s", new_dir)

    def cleanup(self):
        database_names = self.get_database_names()
        for folder in self.base.iterdir():
            try:
                matched = self.find_matches(folder.name, database_names)
                if len(matched) > 1:
                    xyz_files = [self.base / m / f"{m}.xyz" for m in matched if (self.base / m / f"{m}.xyz").exists()]
                    if not (xyz_files[0].exists() and xyz_files[1].exists()):
                        continue
                    with xyz_files[0].open("r") as f:
                        content0 = f.read()
                    with xyz_files[1].open("r") as f:
                        content1 = f.read()
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
        dir = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/")
        for folder in dir.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        for file in subfolder.iterdir():
                            if file.is_symlink() and file.resolve() == original:
                                file.unlink()
                                file.symlink_to(merged)

    def copy_to_db(self):
        dir = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/")
        for folder in dir.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        Database(subfolder).add_calculation()
        self.cleanup()
