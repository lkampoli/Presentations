#!/bin/bash
# One task, one core, default partition, ten minutes walltime
module purge
module load matlab/2020a
matlab -nodesktop -nodisplay -nosplash< polar-plot.m
