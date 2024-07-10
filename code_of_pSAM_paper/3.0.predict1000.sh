#!/bin/sh
#SBATCH --job-name=pred1000
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./3.0.predict1000.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
~/miniconda3/envs/DLstudy25/bin/python 3.0.predict1000.py ${locs}
