#!/bin/sh
#SBATCH --job-name=hs.score
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./4.0.hs.score.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
~/miniconda3/envs/DLstudy25/bin/python 4.0.hs.scores.py ${locs}
