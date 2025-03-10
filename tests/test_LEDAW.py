import sys
sys.path.append('/lustre/work/ws/ws1/tu_zxofv28-my_workspace/LEDAW')
from LEDAW.ledaw_package.nbody_engine import engine_LED_N_body
name = "test8"
main_filenames = [f'/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}/{name}/{name}.out', 
                  f'/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}/subsys_001/subsys_001.out', 
                  f'/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}/subsys_002/subsys_002.out']
alternative_filenames = ['' for _ in range(len(main_filenames))]
method = "DLPNO-CCSD(T)"
conversion_factor = 627.5095
LEDAW_output_path = f'/lustre/work/ws/ws1/tu_zxofv28-my_workspace/{name}'
engine_LED_N_body(main_filenames, alternative_filenames, conversion_factor, method, LEDAW_output_path)