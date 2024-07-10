#!/bin/sh
#SBATCH --job-name=subloc.bestmodel
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./1.0.gridsearch.log
source /home/yukai/.bashrc
conda activate DLstudy25
models=$1
~/miniconda3/envs/DLstudy25/bin/python 1.0.gridsearch.model.cnn_pars.py ../0.data/1.dataset.subloc.Nucleus.txt.train ../0.data/1.dataset.subloc.Nucleus.txt.valid ../0.data/1.dataset.subloc.Nucleus.txt.test $models 2
