import os
import subprocess

class RestartCalculation:
    def __init__(self, dir, name):
        self.dir = dir
        self.name = name
        self.keep = [".xyz", ".inp", ".sh"]

    def _read_files(self, extension):
        """Helper method to read files with a specific extension."""
        files_content = ""
        with open(os.path.join(self.dir, f"{self.name}{extension}"), 'r') as f:
            files_content = f.read()
        return files_content

    def get_error(self):
        """Read error files."""
        return self._read_files("_err.err")

    def get_output(self):
        """Read output files."""
        return self._read_files("_out.out")

    def get_input(self):
        """Read input files."""
        return self._read_files(".inp")

    def get_sh(self):
        """Read shell script files."""
        return self._read_files(".sh")

    def löschen(self):
        """Delete all files that do not match the allowed extensions."""
        for root, _, files in os.walk(self.dir):
            for file in files:
                if not file.endswith(tuple(self.keep)):
                    os.remove(os.path.join(root, file))

    def extrect_current(self):
        # extrahiert akuellen speicher und zeit verbrauch
        sh_content = self.get_sh().splitlines()
        self.memory = None
        self.time = None
        self.cpus = None
        self.scratch = None
        self.error_file = None
        self.output_file = None
        for line in sh_content:
            if "#SBATCH --mem=" in line:
                self.memory = int(line.split("=")[1].strip())
            if "#SBATCH --time=" in line:
                self.time = line.split("=")[1].strip()
            if "#SBATCH --cpus-per-task=" in line:
                self.cpus = int(line.split("=")[1].strip())
            if "#SBATCH --gres=scratch:" in line:
                self.scratch = int(line.split(":")[1].strip())
            if "#SBATCH --error=" in line:
                self.error_file = line.split("=")[1].strip().replace("_err.err", "")
            if "#SBATCH --output=" in line:
                self.output_file = line.split("=")[1].strip().replace("_out.out", "")
            
    def correct_error(self):
        # find error file
        if "JOB EFFICIENCY REPORT" in self.get_output():
            error_files = self.get_error().splitlines()
            for line in error_files:
                if "Out of memory" in line: # TODO: check if this is the right line
                    self.update_mem()
                if "Out of time" in line: # TODO: check if this is the right line
                    self.update_time()
            self.löschen()
            # subprocess.run(["sbatch", sh_path])
    
    def update_mem(self):
        """Update the memory attribute."""
        npros_options = [48, 24, 12, 8, 6, 4, 2, 1]
        mem_options = [180, 370, 750, 1500]
        diff = float('inf')
        ind = None
        for index, mem in enumerate(mem_options):
            if abs(self.memory - mem) < diff:
                diff = abs(self.memory - mem)
                ind = index
        if ind != len(mem_options) - 1:
            self.memory = mem_options[ind+1]
        else:
            diff = float('inf')
            ind = None
            for npros in npros_options:
                if abs(self.cpus - npros) < diff:
                    diff = abs(self.cpus - npros)
                    ind = npros
            if ind != npros_options[-1]:
                self.cpus = npros_options[ind+1]
                self.update_inp()
            else:
                print("No more options available for memory and cpus.")
    
    def update_time(self):
        # in sekunden umwandeln
        time = self.time.split(":")
        time = int(time[0]) * 3600 + int(time[1]) * 60 + int(time[2])
        # zeit * 2
        time = time * 2
        # in stunden umwandeln
        hours = time // 3600
        minutes = (time % 3600) // 60
        seconds = time % 60
        # in hh:mm:ss umwandeln
        self.time = f"{hours}:{minutes:02}:{seconds:02}"
    
    def update_inp(self):
        """Update the input file with new memory and cpus."""
        inp_file_full = self.get_input()
        inp_file = inp_file_full[0].splitlines()
        inp_file_content = ""
        for line in inp_file:
            if "nprocs" in line:
                inp_file_content += f"  nprocs {self.cpus}\n"
            elif "nprocs_group" not in line:
                inp_file_content += line + "\n"

if __name__ == "__main__":
    dir = "/path/to/directory"
    name = "calculation_name"
    restart_calc = RestartCalculation(dir, name)
    restart_calc.extrect_current()
    restart_calc.correct_error()