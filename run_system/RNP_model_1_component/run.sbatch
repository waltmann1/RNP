#!/bin/bash

#SBATCH --job-name=CAnoRNA22
#SBATCH --time=48:00:00
#SBATCH --partition=beagle3
#SBATCH --account=pi-gavoth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100
#SBATCH --exclude=beagle3-0001,beagle3-0029,beagle3-0008


##SBATCH --mail-type=FAIL,BEGIN,END
##SBATCH --mail-user=cwaltmann@uchicago.edu
#module load gsl
#conda activate openmscg
#export OPENMSCG_THREADS=4

module load openmpi-cuda-aware

mpirun -np 1 /home/cwaltmann/kokkos_noswitch/lammps-kk_hiv_cg/src/lmp_kokkos_cuda_mpi -in input.long1 -var SEED 12345 -k on g 1 -sf kk -pk kokkos newton off neigh half



