#!/bin/bash
#SBATCH --job-name=run
#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --account=e89-camc
#SBATCH --partition=standard
#SBATCH --qos=standard

# Load the VASP module
module load PrgEnv-cray
module load vasp/6
module load craype-network-ofi
module load cray-mpich

# Avoid any unintentional OpenMP threading by setting OMP_NUM_THREADS
export OMP_NUM_THREADS=1

# Launch the code.
srun --exclusive --distribution=block:block --hint=nomultithread vasp_std
