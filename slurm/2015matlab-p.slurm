#!/bin/bash
#SBATCH -p physical
#SBATCH --ntasks=8
module load MATLAB 
time matlab -nodesktop -nodisplay -nosplash < tictoc.m
time matlab -nodesktop -nodisplay -nosplash < tictoc-p.m

