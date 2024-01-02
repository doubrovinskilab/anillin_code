#!/bin/bash
#SBATCH -J GA_34
#SBATCH --partition=32GB
#SBATCH --output=log34.out
#SBATCH --error=log34.out
#SBATCH --nodes=1          # number of nodes requested by user
#SBATCH --ntasks=15        # number of total tasks
#SBATCH --time=7-23:00:00   # run time, format: D-H:M:S (max wallclock time)
#SBATCH -A Dobrovinski
#SBATCH --mail-user=mohamadibrahim.cheikh@utsouthwestern.edu
#SBATCH --mail-type=all

module purge
module add shared slurm openmpi/gcc/64/3.1.1 
export PATH=$PATH:/work/biophysics/s203893/petsc_2/valgrind-3.17.0/bin
export PATH=$PATH:/work/biophysics/s203893/petsc_2/petsc_install/bin

mpirun ./Petsc_IBM_Solver -in InputFiles/InputFile_34.txt -pc_type cholesky -use_stokes_pc -ksp_type minres -pc_factor_mat_solver_type mumps -recalculate_mat 50 -write_fluid_info 1 -write_solid_extra_info #-check_symmetry_A
