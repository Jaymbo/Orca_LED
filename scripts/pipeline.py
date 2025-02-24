import os
import glob
import logging
from pathlib import Path
import subprocess
from openbabel import openbabel
from pymol import cmd as pycmd
from database import Database

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class XYZFileHandler:
    def __init__(self, input_xyz):
        self.input_xyz = input_xyz
        self.atom_count, self.atom_data = self.read_xyz_file()
        logging.info(f"Initialized XYZFileHandler with file: {input_xyz}")

    def read_xyz_file(self):
        logging.info(f"Reading XYZ file: {self.input_xyz}")
        with open(self.input_xyz, "r") as xyz_file:
            lines = xyz_file.readlines()
        atom_count = int(lines[0].strip())
        atom_data = [line.replace("\t", " ") for line in lines[2:2 + atom_count]]
        return atom_count, atom_data

    def write_fragment_xyz(self, fragment_atoms, output_filename):
        logging.info(f"Writing fragment XYZ file: {output_filename}")
        with open(output_filename, "w") as output_file:
            output_file.write(f"{len(fragment_atoms)}\n")
            output_file.write("Fragment generated by split_xyz\n")
            output_file.writelines(fragment_atoms)

    def split_xyz(self, fragments_list, output_prefix):
        logging.info(f"Splitting XYZ file into fragments with prefix: {output_prefix}")
        for i, fragment_indices in enumerate(fragments_list):
            fragment_atoms = [self.atom_data[idx] for idx in fragment_indices]
            output_filename = f"{output_prefix}/fragment_{i + 1:03}.xyz"
            self.write_fragment_xyz(fragment_atoms, output_filename)


class Mol2FileHandler:
    def __init__(self, input_mol2):
        self.input_mol2 = input_mol2
        logging.info(f"Initialized Mol2FileHandler with file: {input_mol2}")

    def convert_mol2_to_xyz(self, output_xyz):
        logging.info(f"Converting MOL2 file to XYZ: {self.input_mol2} -> {output_xyz}")
        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats("mol2", "xyz")
        mol = openbabel.OBMol()
        if not ob_conversion.ReadFile(mol, self.input_mol2):
            logging.error(f"Error reading file: {self.input_mol2}")
            raise ValueError(f"Error reading file: {self.input_mol2}")
        if not ob_conversion.WriteFile(mol, output_xyz):
            logging.error(f"Error writing file: {output_xyz}")
            raise ValueError(f"Error writing file: {output_xyz}")

    def extract_fragments_from_mol2(self):
        logging.info(f"Extracting fragments from MOL2 file: {self.input_mol2}")
        fragments = []
        with open(self.input_mol2, "r") as mol2_file:
            inside_atoms = False
            for line in mol2_file:
                if line.startswith("@<TRIPOS>ATOM"):
                    inside_atoms = True
                elif line.startswith("@<TRIPOS>BOND"):
                    inside_atoms = False
                elif inside_atoms:
                    parts = line.split()
                    if len(parts) > 7:
                        fragments.append(parts[7])
        fragment_assignments = {}
        for idx, group in enumerate(fragments):
            fragment_assignments.setdefault(group, []).append(idx)
        return fragment_assignments


class ORCAInputFileCreator:
    def __init__(self, file, fragments=None, header_in=None):
        self.file = file
        self.fragments = fragments
        self.header = header_in or """! DLPNO-CCSD(T) def2-svp def2-svp/C DEF2/J RIJCOSX veryTIGHTSCF TIGHTPNO LED

%mdci DoDIDplot true end

%maxcore 160000

%mdci
  MaxIter 200
end"""
        self.header += """\n%pal \n  nprocs """
        self.xyz_file = self.file.replace(".mol2", ".xyz") if self.file.endswith(".mol2") else self.file
        self.xyz_folder = os.path.dirname(self.xyz_file)
        os.makedirs(self.xyz_folder, exist_ok=True)
        self.frag_len = []
        logging.info(f"Initialized ORCAInputFileCreator with file: {file}")

    def create_inp_file_content(self, charge, npros, xyz_file, fragment_lines=None):
        inp_content = f"{self.header}{npros}\nend\n*XYZfile {charge} 1 {xyz_file}\n\n"
        if fragment_lines:
            inp_content += "".join(fragment_lines)
        return inp_content

    def create_inp_files(self):
        logging.info("Creating ORCA input files")
        if self.file.endswith(".mol2"):
            mol2_handler = Mol2FileHandler(self.file)
            if self.fragments == "":
                self.fragments = mol2_handler.extract_fragments_from_mol2()
            mol2_handler.convert_mol2_to_xyz(self.xyz_file)

        if self.fragments != "-1":
            fragment_lines = self.handle_fragments()
        else:
            fragment_lines = None

        xyz_files = sorted(glob.glob(os.path.join(self.xyz_folder, "*.xyz")))
        if not xyz_files:
            logging.warning("No .xyz files found in the specified folder.")
            return

        ind = 0
        for i, xyz_file_i in enumerate(xyz_files):
            sum_atoms = min(len(open(self.file).read().split("\n")) - 2, 48) if os.path.basename(self.xyz_file).split(".")[0] == os.path.basename(xyz_file_i).split(".")[0] else None
            if sum_atoms:
                self.frag_len = self.frag_len[:i] + [sum_atoms] + self.frag_len[i:]
                ind = i
            self.create_single_inp_file(xyz_file_i, Path(xyz_file_i).parent, self.frag_len, i, fragment_lines)
            path = Database.process_candidate(Path(xyz_file_i.split(".")[0]))
            if path is not None:
                self.create_single_inp_file(path, Path(path).parents[1], self.frag_len, i, fragment_lines)
                if self.fragments != "-1":
                    sh_path = ShellScriptCreator.single_sh_script_erstellen(path, Path(path).parents[1], i, self.frag_len, ind, time="2:00:00", mem=20, main=False)
                else:
                    sh_path = ShellScriptCreator.single_sh_script_erstellen(path, Path(path).parents[1], i, self.frag_len, ind)
                subprocess.run(["sbatch", sh_path])

    def handle_fragments(self):
        logging.info("Handling fragments")
        if self.fragments is None or self.fragments == "":
            self.fragments = FragmentExtractor.extract_fragments(self.file)
        fragment_groups = self.parse_fragments(self.fragments)

        self.calculate_frag_len(fragment_groups)

        xyz_handler = XYZFileHandler(self.xyz_file)
        xyz_handler.split_xyz(fragment_groups, self.xyz_folder)

        fragment_lines = self.create_fragment_lines(fragment_groups)
        return fragment_lines

    def parse_fragments(self, fragments):
        logging.info("Parsing fragments")
        fragment_groups = []
        for fragment in fragments.split(","):
            fragment_indices = []
            for part in fragment.split():
                if "-" in part:
                    start, end = map(int, part.split("-"))
                    fragment_indices.extend(range(start-1, end))
                else:
                    fragment_indices.append(int(part)-1)
            fragment_groups.append(fragment_indices)
        return fragment_groups

    def create_fragment_lines(self, fragment_groups):
        logging.info("Creating fragment lines")
        fragment_lines = ["%geom\n Fragments\n"]
        for i, group in enumerate(fragment_groups, start=1):
            fragment_atoms = " ".join(map(str, group))
            fragment_lines.append(f"  {i} {{{fragment_atoms}}} end\n")
        fragment_lines.append(" end\nend\n")
        return fragment_lines

    def create_single_inp_file(self, xyz_file_i, base, frag_len, i, fragment_lines):
        logging.info(f"Creating single ORCA input file for: {xyz_file_i}")
        base_name = os.path.splitext(os.path.basename(xyz_file_i))[0]
        charge = ChargeCalculator.calculate_charge(
            self.file if os.path.basename(self.xyz_file).split(".")[0] == os.path.basename(xyz_file_i).split(".")[0] else xyz_file_i
        )
        inp_content = self.create_inp_file_content(charge, frag_len[i], xyz_file_i, fragment_lines if self.fragments != "-1" else None)

        os.makedirs(base / base_name, exist_ok=True)
        inp_path = base / f"{base_name}/{base_name}.inp"
        with open(inp_path, "w") as inp_file:
            inp_file.write(inp_content)
        return inp_path

    def calculate_frag_len(self, fragment_groups):
        logging.info("Calculating fragment lengths")
        self.frag_len = [min(len(group), 48) for group in fragment_groups]


class ShellScriptCreator:
    def __init__(self, mem, nprocs, time, path, name, base, scrap=None, main=False):
        self.mem = mem
        self.nprocs = nprocs
        self.time = time
        self.path = path.split(".")[0]
        self.name = name
        self.base = base
        self.scrap = scrap
        self.main = main
        logging.info(f"Initialized ShellScriptCreator for: {name}")

    def create_sh_script_content(self):
        logging.info("Creating shell script content")
        if self.main:
            script_content = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem={self.mem}gb
#SBATCH --ntasks-per-node={self.nprocs}
#SBATCH --gres=scratch:{self.scrap}
#SBATCH --output={self.path}_out.out
#SBATCH --error={self.path}_err.err
#SBATCH --time={self.time}

name={self.name}

workspace_directory={self.base}
orca=/opt/bwhpc/common/chem/orca/6.0.1_shared_openmpi-4.1.6_avx2/orca

echo $name
module load chem/orca/6.0.1
module load mpi/openmpi/4.1
module list

echo "kopieren der .inp"
cp $workspace_directory/$name/$name.inp $SCRATCH

echo "ins scatch"
cd $SCRATCH

echo "ausführen"
$orca $name.inp > $name.out

echo "kopieren ohne tmp"
find $SCRATCH -type f ! -name '*tmp*' -exec cp {{}} $workspace_directory/$name/ \;
"""
        else:
            script_content = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem={self.mem}gb
#SBATCH --ntasks-per-node={self.nprocs}
#SBATCH --time={self.time}
#SBATCH --output={self.path}_out.out
#SBATCH --error={self.path}_err.err

name={self.name}

workspace_directory={self.base}
orca=/opt/bwhpc/common/chem/orca/6.0.1_shared_openmpi-4.1.6_avx2/orca

echo $name
module load chem/orca/6.0.1
module load mpi/openmpi/4.1
module list

echo "ausführen"
$orca $workspace_directory/$name/$name.inp > $workspace_directory/$name/$name.out
"""
        return script_content

    @staticmethod
    def single_sh_script_erstellen(path, base, i, frag_len, index_compound, time="40:00:00", time2="100:00:00", mem=360, scrap=700, main=True):
        name = Path(path).stem
        total_path = base / f"{name}/{name}.sh"
        if i == index_compound and main:
            script_creator = ShellScriptCreator(
                int(mem * 4 * frag_len[i] / 48), frag_len[i], time2, path, name, base, scrap, main
            )
        else:
            script_creator = ShellScriptCreator(
                int(mem * 2 * frag_len[i] / 48), frag_len[i], time, path, name, base
            )
        script_content = script_creator.create_sh_script_content()
        with open(total_path, "w") as file:
            file.write(script_content)
        return total_path


class ChargeCalculator:
    @staticmethod
    def get_valence_electrons(atom):
        valence_dict = {
            "H": 1, "C": 4, "N": 5, "O": 6, "CL": 7, "S": 6, "BR": 7, "P": 5, "SE": 6, "F": 7,
        }
        return valence_dict.get(atom.upper(), 0)

    @staticmethod
    def calculate_charge(filename):
        logging.info(f"Calculating charge for file: {filename}")
        pycmd.reinitialize()
        pycmd.load(filename)
        fragment_name = "all"
        fragment_atoms = pycmd.get_model(fragment_name)
        positive_charge = 0
        negative_atoms = []

        for atom in fragment_atoms.atom:
            valence_electrons = ChargeCalculator.get_valence_electrons(atom.name)
            pycmd.select("selected_atom", f"{fragment_name} and index {atom.index}")
            pycmd.select("selected_atom", "selected_atom xt. 1")
            bond_count = pycmd.count_atoms("selected_atom") - 1

            valence_total = valence_electrons + bond_count
            if atom.name == "H":
                valence_total -= 2
            else:
                valence_total -= 8

            if valence_total < 0:
                negative_atoms.append([atom, valence_total])
            else:
                positive_charge += valence_total

        return positive_charge + ChargeCalculator.expand_negative_charges(fragment_name, negative_atoms)

    @staticmethod
    def expand_negative_charges(fragment_name, negative_atoms):
        total_negative = 0
        while negative_atoms:
            atom, charge = negative_atoms.pop(0)
            pycmd.select("selected_atom", f"{fragment_name} and index {atom.index}")
            pycmd.select("selected_atom", "nbr. selected_atom")
            selected_atoms = pycmd.get_model("selected_atom")

            found = False
            for i in range(len(selected_atoms.atom)):
                if found:
                    break
                for j in range(len(negative_atoms)):
                    if selected_atoms.atom[i].index == negative_atoms[j][0].index:
                        atom2, charge2 = negative_atoms.pop(j)
                        charge2 += 1
                        if charge2 < 0:
                            negative_atoms.append([atom2, charge2])
                        found = True
                        break
            if not found:
                total_negative += charge
            else:
                charge += 1
                if charge < 0:
                    negative_atoms.append([atom, charge])
        return total_negative


class FragmentExtractor:
    @staticmethod
    def extract_fragments(filename):
        logging.info(f"Extracting fragments from file: {filename}")
        pycmd.reinitialize()
        pycmd.load(filename)
        all_atoms = pycmd.get_model("all")
        visited = set()
        indexeses = ""

        for atom in all_atoms.atom:
            if atom.index not in visited:
                temp_select_name = f"temp_fragment_{atom.index}"
                pycmd.select(temp_select_name, f"bm. (id {atom.index})")
                fragment_atoms = pycmd.get_model(temp_select_name)

                indexeses += ","
                for fragment_atom in fragment_atoms.atom:
                    if fragment_atom.index not in visited:
                        visited.add(fragment_atom.index)
                        indexeses += f"{fragment_atom.index} "
        indexeses = indexeses[1:]
        return indexeses
