#!python
#hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions
#Cytosol
import sys
nls2site = {}
with open("../7.eachpro/"+sys.argv[1]) as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpdata = data[1].split(';')
		nls2site[data[0]] = []
		for i in tmpdata:
			nls2site[data[0]].append(i.split('...')[0]+'-'+i.split('...')[1])

def str2region(str1, regions):
	#example, "fafsfasvfg", [2-3, 5-7]
	str1list = list(str1)
	movesites = []
	for i in regions:
		sites = i.split('-')
		sites = [int(j) for j in sites]
		movesites += range(max(sites[0]-8, 1), min(sites[1]+8, len(str1list)))
	tmp_list = []
	for i in range(len(str1list)):
		if i+1 in movesites:
			continue
		else:
			tmp_list.append(str1list[i])
	return("".join(tmp_list))

def randomtrunc(str1, length):
	minimal = 0
	maximun = len(str1) - length - 1
	import random
	start = random.randint(minimal, maximun)
	return(str1[:start]+str1[(start+length):])


out = open("../6.panStrunc/"+sys.argv[1]+".4predict.seq", "w")
with open("../4.panMut/hs.protein.score."+sys.argv[2]) as f:
	#out.write(f.readline().strip('\n')+'\tNLSregions\tTruncatedSequence\tRandomTruncated\n')
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		if len(data[2]) < 50:
			continue
		if data[0] in nls2site:
			seq = str2region(data[2], nls2site[data[0]])
			ranseq = randomtrunc(data[2], 20)
			out.write(line+'\t'+';'.join(nls2site[data[0]])+'\t'+seq+'\t'+ranseq+'\n')


out.close()


import tensorflow as tf
#tf.config.gpu.set_per_process_memory_growth(True)
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
	tf.config.experimental.set_memory_growth(gpu, True)

import sys
from functions import read_data, one_hot_label, roc_file
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2

####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

model1 = model_from_json(open('../2.TrainValidTest/model.'+sys.argv[2]+'.alldata.json').read())
model1.load_weights('../2.TrainValidTest/model'+sys.argv[2]+'.alldata.h5')

max_len = 2000
##human.dataset
all_data, all_label, all_line = [], [], []
with open("../6.panStrunc/"+sys.argv[1]+".4predict.seq") as f:
	tittle = f.readline().strip('\n')
	for line in f:
		data = line.strip('\n').split('\t')
		if len(data[5]) < max_len:
			all_line.append(line.strip('\n'))
			all_data.append(data[5])


train_data = all_data
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score = model1.predict(train_ohe, batch_size = 1)

all_data = []
with open("../6.panStrunc/"+sys.argv[1]+".4predict.seq") as f:
	tittle = f.readline().strip('\n')
	for line in f:
		data = line.strip('\n').split('\t')
		if len(data[6]) < max_len:
			all_data.append(data[6])


train_data = all_data
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score1 = model1.predict(train_ohe, batch_size = 1)


out = open("../6.panStrunc/"+sys.argv[1]+".4predict.results", "w")
out.write(tittle+'\tScore2\tScore3'+'\n')
for i in range(len(all_line)):
	out.write(all_line[i]+'\t'+str(score[i][1])+'\t'+str(score1[i][1])+'\n')


out.close()

