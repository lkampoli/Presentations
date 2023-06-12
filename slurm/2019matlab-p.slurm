#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=1:0:0

module purge
module load matlab/2020a

time matlab -nodesktop -nodisplay -nosplash < tictoc.m
time matlab -nodesktop -nodisplay -nosplash < tictoc-p.m

