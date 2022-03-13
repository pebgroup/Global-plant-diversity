#!/bin/bash
#SBATCH --account GlobalDrivers
#SBATCH --job-name=var_select_master
##SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
##SBATCH --ntasks=



# model setup
nmodels=5 # number of models per script: 
starts=(1 6 11 16 21 26 31 36 41 46 51 56 61 66 71 76 81 86 91 96)
var_type='mdr'

export nmodels
export var_type


for ((i=0; i<=19; i++)) do
  this_start=${starts[${i}]}
  export this_start
  echo Submitting bash script with starting point $this_start
  sbatch var_selection_sequential.sh
done
