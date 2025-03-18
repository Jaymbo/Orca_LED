import sys
from pathlib import Path

BASE_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(BASE_PATH / 'scripts'))

from pipeline import ORCAInputFileCreator

file = BASE_PATH / "tests/BrBr_FCCH.xyz"
fragment_input = ""
header_input = """! DLPNO-CCSD(T1) aug-cc-pVTZ DEFGRID3 RIJCOSX def2/J aug-cc-pVQZ/C verytightscf TightPNO LED

%maxcore 160000
%scf maxiter 999 end
%basis auxj "autoaux" auxc "autoaux" autoauxlmax true end
%mdci TCutPNO 1e-6 printlevel 3 end"""

orca_input_creator = ORCAInputFileCreator(str(file), fragment_input, header_input)
orca_input_creator.create_inp_files()