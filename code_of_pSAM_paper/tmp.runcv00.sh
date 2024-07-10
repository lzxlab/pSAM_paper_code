#!/bin/sh
#SBATCH --job-name=subloc.cv
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -p gpuPartition
#SBATCH --gres=gpu:tesla:1
#SBATCH -o ./1.0.runcv00.log
source /home/yukai/.bashrc
conda activate DLstudy25
loc=$1

/home/rstudio_cu02/anaconda3/envs/R4.0/bin/Rscript ./0.10.train.valid.test.R $loc
~/miniconda3/envs/DLstudy25/bin/python ./1.1.gridsearch.model.cnn.py ../0.data/1.dataset.subloc.${loc}.txt.train ../0.data/1.dataset.subloc.${loc}.txt.valid ../0.data/1.dataset.subloc.${loc}.txt.test CNNBiLSTMAtten1 2


