#!/bin/bash
#SBATCH --account GlobalDrivers
#SBATCH --job-name=mod_master
##SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=melanie.tietje@bio.au.dk
#SBATCH --partition normal
#SBATCH --mem-per-cpu=1gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
##SBATCH --ntasks=



# setup for 32 jobs, runs in <1h
nmodels=32704 # number of models per script: 31*32704 + 1*32705
starts=(1 32705 65409 98113 130817 163521 196225 228929 261633 294337 327041 359745 392449 425153 457857 490561 523265 555969 588673 621377 654081 686785 719489 752193 784897 817601 850305 883009 915713 948417 981121 1013825)

export nmodels


for ((i=0; i<=31; i++)) do
  this_start=${starts[${i}]}
  export this_start
  echo Submitting bash script with starting point $this_start
  sbatch mod_selection_sequential.sh
done
