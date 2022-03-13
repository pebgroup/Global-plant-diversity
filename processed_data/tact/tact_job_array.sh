#!/bin/bash
# submit_array.sh

#SBATCH --account GlobalDrivers
#SBATCH --job-name=GD_TACT
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=145gb
#SBATCH --cpus-per-task 1
#SBATCH --time 90:00:00
#SBATCH --array=400-450


singularity exec tact.sif tact_build_taxonomic_tree goodsp.csv --output goodsp.taxonomy.tre &&

source ~/miniconda3/bin/activate tact

echo -e "\nrunning TACT\n"


singularity exec tact.sif  tact_add_taxa --backbone gbmb_matched_no_misplaced.tre --taxonomy goodsp.taxonomy.tre --output gbmb_matched_no_misplaced_$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID.tacted --verbose
