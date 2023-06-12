#!/bin/bash
#SBATCH --nodes 1
#SBATCH --account hpcadmingpgpu
# Use a project ID that has gpgpu access.
#SBATCH --partition gpgpu
#SBATCH --gres=gpu:p100:4
#SBATCH --time 01:00:00
#SBATCH --cpus-per-task=24

module purge
module load fosscuda/2019b
module load tensorflow/2.1.0-python-3.7.4
module load Tensorflow/1.8.0-intel-2017.u

cd benchmarks/scripts/tf_cnn_benchmarks
python tf_cnn_benchmarks.py --num_gpus=4 --batch_size=64 --model=resnet50 --variable_update=parameter_server

