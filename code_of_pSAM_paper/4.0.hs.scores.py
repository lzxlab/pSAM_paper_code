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

model1 = model_from_json(open('../2.TrainValidTest/model.'+locs+'.alldata.json').read())
model1.load_weights('../2.TrainValidTest/model'+locs+'.alldata.h5')

uni2seq, uni2genes = {}, {}
with open("../0.data/uniprot.fasta") as f:
	for line in f:
		if re.search("^>", line):
			tmp = "haha"
			data = line.strip('\n').split('OS=')[1].split(" ")[0]
			if len(line.strip('\n').split('GN=')) < 2:
				continue
			if data == "Homo":
				tmp = "lala"
				unis = line.strip('\n').split('|')[1]
				uni2genes[unis] = line.strip('\n').split('GN=')[1].split(" ")[0]
				uni2seq[unis] = ""
		else:
			if tmp == "haha":
				continue
			uni2seq[unis] += line.strip('\n')

max_len = 2000
##human.dataset
all_line, all_mut = [], []
for line in uni2seq:
	if len(uni2seq[line]) < max_len:
		all_line.append(line+'\t'+uni2genes[line]+'\t'+uni2seq[line])
		all_mut.append(uni2seq[line])


train_data = all_mut
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe, batch_size = 1)


out = open("../4.panMut/hs.protein.score."+locs, "w")
for i in range(len(all_line)):
	out.write(all_line[i]+'\t'+str(score_mut[i][1])+'\n')


out.close()


print("-------------------Finished "+sys.argv[1]+" -----------------------")
