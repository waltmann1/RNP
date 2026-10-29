#!/bin/bash

#SBATCH --job-name=CAnoRNA22
#SBATCH --time=48:00:00
#SBATCH --partition=beagle3
#SBATCH --account=pi-gavoth
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:2
#SBATCH --exclude=beagle3-0001,beagle3-0029,beagle3-0008,beagle3-0017
#SBATCH --constraint=a100


##SBATCH --mail-type=FAIL,BEGIN,END
##SBATCH --mail-user=cwaltmann@uchicago.edu
#module load gsl
#conda activate openmscg
#export OPENMSCG_THREADS=4

module load openmpi-cuda-aware

export CUDA_VISIBLE_DEVICES=0,1


mpirun --cpu-set 0,1  --bind-to core -np 2 /beagle3/gavoth/cwaltmann/kokkos_double_gauss/lammps-kk_hiv_cg/src/lmp_kokkos_cuda_mpi -in input.long32 -var SEED 12345 -k on g 2 -sf kk -pk kokkos newton off neigh half

