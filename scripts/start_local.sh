#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=187gb
#SBATCH --time=20:00:00

name=$2

workspace_directory=$1
orca=/opt/bwhpc/common/chem/orca/6.0.0_shared_openmpi-4.1.6_avx2/orca

echo $name
module load chem/orca/6.0.1
module load mpi/openmpi/4.1
module list


echo "ausführen"
$orca $workspace_directory/$name/$name.inp > $workspace_directory/$name/$name.out
