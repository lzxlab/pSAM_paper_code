#!/bin/sh
#SBATCH --job-name=Alldata
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./2.1.alldata.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
~/miniconda3/envs/DLstudy25/bin/python 2.1.getModel4alldata.py ../0.data/1.dataset.subloc.${locs}
