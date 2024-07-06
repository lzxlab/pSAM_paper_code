#!python
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
import re
####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

locs = "Nucleus"
subfiles = "mut.data"

model1 = model_from_json(open('../2.TrainValidTest/model.'+locs+'.alldata.json').read())
model1.load_weights('../2.TrainValidTest/model'+locs+'.alldata.h5')

uni2seq = {}
with open("../0.data/uniprot.subloc.wholedataset.redeal.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		uni2seq[data[0]] = data[2]



print("-------------------Runing "+ subfiles + "of " + locs + " -----------------------")
max_len = 2000
##human.dataset
all_line, all_wt, all_mut = [], [], []
with open("../10.NucleRes/"+subfiles) as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] not in uni2seq:
			continue
		if len(uni2seq[data[0]]) < max_len:
			ref = re.split('\d+',data[7])[0]
			alt = re.split('\d+',data[7])[1]
			pos = int(re.findall('\d+',data[7])[0])
			if pos < len(uni2seq[data[0]]):
				if uni2seq[data[0]][pos-1] == ref:
					tmpmut = list(uni2seq[data[0]])
					tmpmut[pos-1] = alt
					all_line.append(line.strip('\n'))
					all_wt.append(uni2seq[data[0]])
					all_mut.append("".join(tmpmut))


train_data = all_wt
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_wt = model1.predict(train_ohe, batch_size = 1)

train_data = all_mut
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe, batch_size = 1)


out = open("../10.NucleRes/collected.mut.maf.score."+locs, "w")
for i in range(len(all_line)):
	out.write(all_line[i]+'\t'+str(score_wt[i][1])+'\t'+str(score_mut[i][1])+'\n')


out.close()


print("-------------------Finished "+sys.argv[1]+" -----------------------")
