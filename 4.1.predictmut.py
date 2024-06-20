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

locs = sys.argv[1]
subfiles = sys.argv[2]

model1 = model_from_json(open('../2.TrainValidTest/model.'+locs+'.alldata.json').read())
model1.load_weights('../2.TrainValidTest/model'+locs+'.alldata.h5')

uni2seq, uni2score = {}, {}
with open("../4.panMut/hs.protein.score."+locs) as f:
	for line in f:
		data = line.strip('\n').split('\t')
		uni2seq[data[0]] = data[2]
		uni2score[data[0]] = data[3]



print("-------------------Runing "+ subfiles + "of " + locs + " -----------------------")
max_len = 2000
##human.dataset
score_train, all_line, all_mut = [], [], []
with open("../4.panMut/"+subfiles) as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if len(uni2seq[data[9]]) < max_len:
			ref = data[6][2]
			alt = re.split('\d+',data[6])[1]
			pos = int(re.findall('\d+',data[6])[0])
			if pos < len(uni2seq[data[9]]):
				if uni2seq[data[9]][pos-1] == ref:
					if data[-1] in uni2score:
						tmpmut = list(uni2seq[data[9]])
						tmpmut[pos-1] = alt
						all_line.append(line.strip('\n'))
						score_train.append(uni2score[data[9]])
						all_mut.append("".join(tmpmut))


train_data = all_mut
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe, batch_size = 1)


out = open("../4.panMut/merged.mut.maf.score."+locs, "a")
for i in range(len(all_line)):
	out.write(all_line[i]+'\t'+str(score_train[i])+'\t'+str(score_mut[i][1])+'\n')


out.close()


print("-------------------Finished "+sys.argv[1]+" -----------------------")
