#!/bin/bash
#SBATCH --ntasks=1
module purge
source /usr/local/module/spartan_old.sh
module load MATLAB/2016a
matlab -nodesktop -nodisplay -nosplash< polar-plot.m
