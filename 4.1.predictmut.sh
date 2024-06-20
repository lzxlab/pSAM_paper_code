#!/bin/sh
#SBATCH --job-name=pred.mut
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./4.1.predictmut.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
files=$2
~/miniconda3/envs/DLstudy25/bin/python 4.1.predictmut.py ${locs} ${files}
