#!/bin/sh
#SBATCH --job-name=calcCPT
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./7.0.predict.allpros.log
source /home/yukai/.bashrc
conda activate DLstudy25
pro=$1
locs=$2
cutoff=$3
~/miniconda3/envs/DLstudy25/bin/python 7.0.predict.allpros.py ${pro} ${locs} ${cutoff}
