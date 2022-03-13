#!/bin/bash
#SBATCH --account GlobalDrivers
#SBATCH --job-name=GBM
##SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
##SBATCH --ntasks=


source ~/miniconda3/bin/activate R-env-4

echo Running script $i with starting point ${this_start[${i}]}
Rscript gbm.R $nmodels $this_start $var_type > log_"$SLURM_JOB_ID"_"$this_start"_"$nmodels".txt
