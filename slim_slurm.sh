#!/bin/bash
#SBATCH --job-name=sparse_moreplots
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --account=modelscape
#SBATCH --mem-per-cpu=200G

module load  arcc/1.0 slim/4.0.1 gcc/12.2.0 gsl/2.7.1 r

Rscript --vanilla /gscratch/emcfarl2/Slimulations/slim_working_script_5pops.R