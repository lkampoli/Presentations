#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=3-18:00:00

##SBATCH --mail-type=begin        # send email when job begins
##SBATCH --mail-type=end          # send email when job ends
##SBATCH --mail-type=fail         # send email if job fails
##SBATCH --mail-user=shanrahan@student.unimelb.edu.au #Email Address

#SBATCH --account punim0394
#SBATCH --qos=gpgpumse
#SBATCH --partition extremecfd
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

module purge

#module load python/3.7.4
#module load scipy-bundle/2019.03
module load singularity/3.5.3

CONTAINER=$HOME/tf_21.12-py3.simg
# cd ./main/NS

# Launch multiple process python code
# echo "Searching for mentions"
# time srun -n 1 python3 PINNs_006.py >>PINN_log_006o.txt
# singularity run --nv /usr/local/singularity_images/Tensorflow/nvidia-tensorflow-21.12-tf1-py3.sif python3 PINNs_008_6_Dim.py>>PINN_log_008_1335x301_E800.txt
