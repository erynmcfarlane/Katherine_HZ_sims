#!/bin/bash
#SBATCH --job-name=sparse_moreplots
#SBATCH --nodes=1
#SBATCH --time=2-12:00:00
#SBATCH --account=evolgen
#SBATCH --mem-per-cpu=200G

module load arcc/1.0 gcc/14.2.0 gsl/2.7.1 r

for i in $(seq 1 100); do

 replicate_dir="/gscratch/erynmcfarlane/Slimulations/replicate_$i"
  mkdir -p $replicate_dir
  cd $replicate_dir
Rscript --vanilla /gscratch/erynmcfarlane/Slimulations/slim_working_script_5pops_50individuals.R > $replicate_dir/output.txt 2>&1

done