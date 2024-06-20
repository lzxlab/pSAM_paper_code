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
max_len = 2000
model1 = model_from_json(open('../2.TrainValidTest/model.'+locs+'json').read())
model1.load_weights('../2.TrainValidTest/model.'+locs+'.h5')

seqs, lines = [],  []
with open("../0.data/1.dataset.subloc."+locs+".txt.compare") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		seqs.append(data[3])
		lines.append(line.strip('\n'))


train_data = seqs
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe)


out = open("../3.existingtools/sel.protein.score."+locs, "w")
for i in range(len(lines)):
	out.write(lines[i]+'\t'+str(score_mut[i][1])+'\n')


out.close()

