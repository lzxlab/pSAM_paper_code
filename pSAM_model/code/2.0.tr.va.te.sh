#!/bin/sh
#SBATCH --job-name=Tr.Va.Te
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./2.0.tr.va.te.log
source /home/yukai/.bashrc
conda activate DLstudy25
locs=$1
~/miniconda3/envs/DLstudy25/bin/python 2.0.getModel.tr.va.te.py ../0.data/1.dataset.subloc.${locs}.txt.train ../0.data/1.dataset.subloc.${locs}.txt.valid ../0.data/1.dataset.subloc.${locs}.txt.test
