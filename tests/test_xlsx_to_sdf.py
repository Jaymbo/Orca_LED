import sys
from pathlib import Path

BASE_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(BASE_PATH / 'scripts'))

from xlsx_to_sdf import SdfXyzMerger

# Beispielaufruf
merger = SdfXyzMerger("/lustre/work/ws/ws1/tu_zxofv28-my_workspace/tests/input.sdf", "/lustre/work/ws/ws1/tu_zxofv28-my_workspace/tests/input.xyz","/lustre/work/ws/ws1/tu_zxofv28-my_workspace/tests")
merger.run()