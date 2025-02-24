from pathlib import Path
import datetime
from shutil import copy2
import numpy as np  # pylint: disable=import-error
from scipy.optimize import linear_sum_assignment  # pylint: disable=import-error
import subprocess


class Database:
    """
    The Database class manages xyz files and associated output files.
    Each dataset is stored in a folder whose name is composed of the atoms,
    the header, and the timestamp.
    """

    def __init__(self, dir: Path) -> None:
        self.base = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/database/")
        self.dir: Path = dir
        self.filename: Path
        self.header_filename: Path
        self.get_file_paths()
        self.rmsd_value = 0.001
        self.ends = ["_out.out", ".out", ".densities", "_err.err", ".gbw", ".property.txt", ".bibtex", ".cube", ".densitiesinfo", ".sh"]

    def get_file_paths(self):
        """
        xyz file is in parent file with the name dir.stem and the header file is in dir and .inp
        """
        self.filename = self.dir.parent / f"{self.dir.stem}.xyz"
        self.header_filename = self.dir / f"{self.dir.stem}.inp"

    def get_filecontent(self) -> str:
        with self.filename.open("r") as f:
            return f.read()

    def atoms_from_filecontent(self, content: str) -> list:
        lines = [line for line in content.splitlines() if line.strip()]
        # Assume atoms listed from the 3rd nonempty line onward.
        return [line.split()[0] for line in lines[2:]]

    def header_from_file(self) -> str:
        header_lines = []
        with self.header_filename.open("r") as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith("%maxcore"):
                    continue
                if stripped.startswith("%pal"):
                    break
                header_lines.append(stripped)
        return "".join(header_lines)

    def atoms_to_str(self, atoms: list) -> str:
        sorted_atoms = sorted(atom.upper() for atom in atoms)
        unique, counts = np.unique(sorted_atoms, return_counts=True)
        return "".join(f"{a}{c}" for a, c in zip(unique, counts))

    def header_to_str(self, header: str) -> str:
        return "".join(char for char in header.upper() if char.isalnum())

    def date_to_str(self) -> str:
        return datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")

    def create_filename(self, atoms: list, header: str) -> tuple[Path, Path]:
        atoms_str = self.atoms_to_str(atoms)
        header_str = self.header_to_str(header)
        date_str = self.date_to_str()
        name = f"{atoms_str}_{header_str}_{date_str}"
        dirpath = self.base / name
        filepath = dirpath / name
        return dirpath, filepath

    def get_database_names(self) -> list:
        return [f.name for f in self.base.iterdir() if f.is_dir()]

    def find_matches(self, current_name: str, database_names: list) -> list:
        name_parts = current_name.split("_")[:2]
        return [
            name for name in database_names
            if name.split("_")[:2] == name_parts
        ]

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
        # Process first xyz
        lines1 = [line.strip() for line in xyz1.splitlines() if line.strip()]
        data1 = lines1[2:]
        atoms1 = [line.split()[0] for line in data1]
        coords1 = np.array([[float(x) for x in line.split()[1:4]] for line in data1])

        # Process second xyz
        lines2 = [line.strip() for line in xyz2.splitlines() if line.strip()]
        data2 = lines2[2:]
        atoms2 = [line.split()[0] for line in data2]
        coords2 = np.array([[float(x) for x in line.split()[1:4]] for line in data2])

        if sorted(atoms1) != sorted(atoms2):
            return 100.0

        D1 = self.compute_distance_matrix(coords1)
        D2 = self.compute_distance_matrix(coords2)
        row_ind, col_ind, _ = self.match_distance_matrices(D1, D2)
        return self.rmsd_of_distance_matrices(D1, D2, row_ind, col_ind)

    def molecule_exists(self, candidate_xyz: str, matched: list) -> bool:
        xyz_files = [self.base / folder / f"{folder}.xyz" for folder in matched]
        rmsd_list = []
        for file_path in xyz_files:
            with file_path.open("r") as f:
                content = f.read()
            rmsd_list.append(self.rmsd(candidate_xyz, content))
        print("Calculated RMSD values:", rmsd_list)
        return bool(rmsd_list) and min(rmsd_list) < self.rmsd_value

    def insert(self, dirpath: Path, filepath: Path) -> None:
        dirpath.mkdir(parents=True, exist_ok=True) 
        filepath = Path(str(filepath.parent) + "/" + str(filepath.stem))
        xyz_path = filepath.with_suffix(".xyz")
        for end in self.ends:
            out_path = Path(str(filepath) + end)
            with out_path.open("w") as f:
                f.write("test")
        
        copy2(self.filename, xyz_path)
        # subprocess.run(["sbatch", sh_path])
        self.create_symlink(xyz_path)


    def create_symlink(self, out_path: Path) -> None:
        """
        für jede datei in out_path soll ein link zu dieser datei oder folder erstellt werden in dem ordner in dem filename liegt
        """
        for file in out_path.parent.iterdir():
            end = self.file_end(file)
            if end is not None:
                symlink = Path(f"{self.header_filename.parent}/{self.header_filename.stem}{end}")
                if not symlink.exists():
                    symlink.symlink_to(file)
                    symlink.chmod(0o444)
            
    def file_end(self, filename: Path) -> str:
        for match in self.ends:
            if filename.name.endswith(match):
                return match
        return

    @classmethod
    def process_candidate(cls, dir: Path) -> bool:
        """
        Executes all necessary steps:
        - Reads the xyz and header files,
        - Extracts the atoms and prepares the naming string,
        - Checks if the molecule is already in the database,
        - Inserts a new molecule or, if it already exists, still creates a symlink
          to the .out file.
        """
        db = cls(dir)
        candidate_xyz = db.get_filecontent()
        atoms = db.atoms_from_filecontent(candidate_xyz)
        header = db.header_from_file()
        new_dir, new_filepath = db.create_filename(atoms, header)
        database_names = db.get_database_names()
        matched = db.find_matches(new_filepath.name, database_names)
        print("Matched folders:", matched)

        if not db.molecule_exists(candidate_xyz, matched):
            db.insert(new_dir, new_filepath)
            return True
        else:
            # If the molecule exists, use the first matching database entry.
            existing_folder = db.base / matched[0]
            out_path = existing_folder / (matched[0] + ".out")
            db.create_symlink(out_path)
            return False
        
    def add_calculation(self):
        """
        Adds a calculation to the database.
        Imports a folder path and an xyz file from that folder.
        The header file is the first .inp file found in the folder.
        Copies every file from the folder to the new database folder.
        """

        # Create a Database instance using the found xyz and header files.
        candidate_xyz = self.get_filecontent()
        atoms = self.atoms_from_filecontent(candidate_xyz)
        header = self.header_from_file()
        new_dir, new_filepath = self.create_filename(atoms, header)

        # Create the new database folder.
        new_dir.mkdir(parents=True, exist_ok=True)

        # Copy every file from the original directory to the new database folder.
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
        
        print(f"Calculation files copied to database folder: {new_dir}")

    def cleanup(self):
        """
        merge files with rmsd < 0.0001
        """
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
            except:
                print(f"Error in folder {folder.name}")
    
    def syslink_merge(self, original: Path, merged: Path):
        """
        change all syslinks from the original file to the merged file
        """
        dir = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/")
        for folder in dir.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        for file in subfolder.iterdir():
                            if file.is_symlink() and file.resolve() == original:
                                file.unlink()
                                file.symlink_to(merged)
                                file.chmod(0o444)

    def copy_to_db(self):
        """
        copy every file in the folder to the database
        """
        dir = Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/")
        for folder in dir.iterdir():
            if folder.is_dir():
                for subfolder in folder.iterdir():
                    if subfolder.is_dir() and any(file.suffix == ".xyz" for file in folder.iterdir()):
                        Database(subfolder).add_calculation()


# Example usage:
if __name__ == "__main__":

    # simple test
    # Database.process_candidate(Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/MRD_A_docked_4/MRD_A_docked_4"))

    # copy everything to the database
    # Database(Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/MRD_A_docked_4/MRD_A_docked_4")).copy_to_db()
    
    # nach db merge cleanup
    # Database(Path("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/MRD_A_docked_4/MRD_A_docked_4")).cleanup()
    print("done")