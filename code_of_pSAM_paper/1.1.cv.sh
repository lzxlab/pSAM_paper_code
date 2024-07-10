#!/bin/sh
#SBATCH --job-name=CVs
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./1.1.cv.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
nums=$2
~/miniconda3/envs/DLstudy25/bin/python 1.1.cv4eachdataset.py ../0.data/1.dataset.subloc.${locs}.txt.train.cv${nums} ../0.data/1.dataset.subloc.${locs}.txt.valid.cv${nums}
