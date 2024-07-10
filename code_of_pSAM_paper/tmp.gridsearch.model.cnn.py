#!python
import tensorflow as tf
#tf.config.gpu.set_per_process_memory_growth(True)
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
	tf.config.experimental.set_memory_growth(gpu, True)

import sys
from functions import read_data, one_hot_label
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from keras.layers.core import Dense, Dropout, Activation
from keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2

max_len = 2000
file1 = sys.argv[1]
locs = file1.split("1.dataset.subloc.")[1].split(".txt.")[0]
train_data, train_label = read_data(file1)
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
file2 = sys.argv[2]
vali_data, vali_label = read_data(file2)
vali_data = one_hot_label(vali_data)
vali_data = sequence.pad_sequences(vali_data, maxlen=max_len, padding='post',value=0)
vali_data = np.array(vali_data)
file3 = sys.argv[3]
test_data, test_label = read_data(file3)
test_data = one_hot_label(test_data)
test_data = sequence.pad_sequences(test_data, maxlen=max_len, padding='post',value=0)
test_data = np.array(test_data)

# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
val_ohe = np_utils.to_categorical(vali_data, num_classes=22)
test_ohe = np_utils.to_categorical(test_data, num_classes=22)
print(train_ohe.shape, val_ohe.shape, test_ohe.shape)
#all_data = all_data.astype(np.int)
label_train = np.array(train_label, dtype = int)
label_vali = np.array(vali_label, dtype = int)
label_test = np.array(test_label, dtype = int)
y_train = np_utils.to_categorical(label_train, num_classes=2)
y_vali = np_utils.to_categorical(label_vali, num_classes=2)
y_test = np_utils.to_categorical(label_test, num_classes=2)

####################Model building####################
def residual_block1(data, filters, sizes):
	bn1 = BatchNormalization()(data)
	act1 = Activation('relu')(bn1)
	conv1 = Conv1D(filters, sizes[0], padding='valid', kernel_regularizer=l2(0.001))(act1)
	#bottleneck convolution
	bn2 = BatchNormalization()(conv1)
	act2 = Activation('relu')(bn2)
	conv2 = Conv1D(filters, sizes[1], padding='valid', kernel_regularizer=l2(0.001))(act2)
	#third convolution
	bn3 = BatchNormalization()(conv2)
	act3 = Activation('relu')(bn3)
	conv3 = Conv1D(filters, sizes[2], padding='valid', kernel_regularizer=l2(0.001))(act3)
	#skip connection
	return conv3

def residual_block2(data, filters, sizes):
	bn1 = BatchNormalization()(data)
	act1 = Activation('relu')(bn1)
	conv1 = Conv1D(filters, sizes[0], padding='valid', kernel_regularizer=l2(0.001))(act1)
	#bottleneck convolution
	bn2 = BatchNormalization()(conv1)
	act2 = Activation('relu')(bn2)
	conv2 = Conv1D(filters, sizes[1], padding='valid', kernel_regularizer=l2(0.001))(act2)
	#third convolution
	bn3 = BatchNormalization()(conv2)
	act3 = Activation('relu')(bn3)
	conv3 = Conv1D(filters, sizes[2], padding='valid', kernel_regularizer=l2(0.001))(act3)
	#forth convolution
	bn4 = BatchNormalization()(conv3)
	act4 = Activation('relu')(bn4)
	conv4 = Conv1D(filters, sizes[3], padding='valid', kernel_regularizer=l2(0.001))(act4)
	#skip connection
	return conv4

if sys.argv[4] == "CNN1":
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	densenodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in densenodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block1(x_input, co, (4,8,12))
					cnn2 = residual_block1(x_input, co, (12,4,8))
					cnn3 = residual_block1(x_input, co, (8,12,4))
					con123 = Concatenate()([cnn1, cnn2, cnn3])
					bn11 = BatchNormalization()(con123)
					act11 = Activation('relu')(bn11)
					maxpool11 = MaxPooling1D(1979)(act11)
					drop = Dropout(dr)(maxpool11)
					# softmax classifier
					fla = Flatten()(drop)
					dense = Dense(do)(fla)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "CNN2":
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	densenodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in densenodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block2(x_input, co, (1, 2, 4, 8))
					cnn2 = residual_block2(x_input, co, (2, 8, 1, 4))
					cnn3 = residual_block2(x_input, co, (4, 1, 8, 2))
					cnn4 = residual_block2(x_input, co, (8, 4, 2, 1))
					con1234 = Concatenate()([cnn1, cnn2, cnn3, cnn4])
					bn11 = BatchNormalization()(con1234)
					act11 = Activation('relu')(bn11)
					maxpool11 = MaxPooling1D(1989)(act11)
					drop = Dropout(dr)(maxpool11)
					# softmax classifier
					fla = Flatten()(drop)
					dense = Dense(do)(fla)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "BiLSTM":
	repeatime = int(sys.argv[5])
	embnodes = [64]
	drate = [0.5]
	lstmnodes = [128]
	for co in embnodes:
		for dr in drate:
			for do in lstmnodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,))
					emb = Embedding(22, co, input_length=max_len)(x_input)
					bi_rnn = Bidirectional(LSTM(do, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(emb)
					x = Dropout(dr)(bi_rnn)
					# softmax classifier
					x_output = Dense(2, activation='softmax')(x)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_data, y_train, batch_size=64, epochs=50, validation_data=(vali_data, y_vali),callbacks=cb)
					score = model1.predict(train_data)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(vali_data)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_data)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "CNNBiLSTM1":
	from tensorflow.keras.layers import Dense, Dropout, Activation
	from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
	from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
	from tensorflow.keras.models import Model
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	lstmnodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in lstmnodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block1(x_input, co, (4,8,12))
					cnn2 = residual_block1(x_input, co, (12,4,8))
					cnn3 = residual_block1(x_input, co, (8,12,4))
					con123 = Concatenate()([cnn1, cnn2, cnn3])
					bn11 = BatchNormalization()(con123)
					act11 = Activation('relu')(bn11)
					bi_rnn = Bidirectional(LSTM(do, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)
					drop = Dropout(dr)(bi_rnn)
					# softmax classifier
					dense = Dense(do)(drop)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "CNNBiLSTM2":
	from tensorflow.keras.layers import Dense, Dropout, Activation
	from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
	from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
	from tensorflow.keras.models import Model
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	lstmnodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in lstmnodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block2(x_input, co, (1, 2, 4, 8))
					cnn2 = residual_block2(x_input, co, (2, 8, 1, 4))
					cnn3 = residual_block2(x_input, co, (4, 1, 8, 2))
					cnn4 = residual_block2(x_input, co, (8, 4, 2, 1))
					con1234 = Concatenate()([cnn1, cnn2, cnn3, cnn4])
					bn11 = BatchNormalization()(con1234)
					act11 = Activation('relu')(bn11)
					bi_rnn = Bidirectional(LSTM(do, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)
					drop = Dropout(dr)(bi_rnn)
					# softmax classifier
					dense = Dense(do)(drop)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "CNNBiLSTMAtten1":
	from tensorflow.keras.layers import Dense, Dropout, Activation
	from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
	from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
	from tensorflow.keras.models import Model
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	lstmnodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in lstmnodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block1(x_input, co, (4,8,1))
					cnn2 = residual_block1(x_input, co, (1,4,8))
					cnn3 = residual_block1(x_input, co, (8,1,4))
					con123 = Concatenate()([cnn1, cnn2, cnn3])
					bn11 = BatchNormalization()(con123)
					act11 = Activation('relu')(bn11)
					bi_rnn = Bidirectional(LSTM(do, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)
					drop = Dropout(dr)(bi_rnn)
					atten = Attention(causal=True)([drop,drop])
					# softmax classifier
					dense = Dense(do)(atten)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
elif sys.argv[4] == "CNNBiLSTMAtten2":
	from tensorflow.keras.layers import Dense, Dropout, Activation
	from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
	from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
	from tensorflow.keras.models import Model
	repeatime = int(sys.argv[5])
	cnnnodes = [64]
	drate = [0.5]
	lstmnodes = [128]
	for co in cnnnodes:
		for dr in drate:
			for do in lstmnodes:
				for i in range(repeatime):
					out1 = open("../1.cvroc/"+locs+"loc"+"_"+sys.argv[4]+"_"+str(co)+"_"+str(dr)+"_"+str(do)+"_"+str(i)+".roc", "w")
					x_input = Input(shape=(2000,22))
					cnn1 = residual_block2(x_input, co, (1, 2, 4, 8))
					cnn2 = residual_block2(x_input, co, (2, 8, 1, 4))
					cnn3 = residual_block2(x_input, co, (4, 1, 8, 2))
					cnn4 = residual_block2(x_input, co, (8, 4, 2, 1))
					con1234 = Concatenate()([cnn1, cnn2, cnn3, cnn4])
					bn11 = BatchNormalization()(con1234)
					act11 = Activation('relu')(bn11)
					bi_rnn = Bidirectional(LSTM(do, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)
					drop = Dropout(dr)(bi_rnn)
					atten = Attention(causal=True)([drop,drop])
					# softmax classifier
					dense = Dense(do)(drop)
					x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)
					model1 = Model(inputs=x_input, outputs=x_output)
					model1.summary()
					model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
					cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
					history1 = model1.fit(train_ohe, y_train, batch_size=64, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)
					score = model1.predict(train_ohe)
					t_label = [int(i) for i in train_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc1 = roc_auc_score(t_label,t_score)
					score = model1.predict(val_ohe)
					t_label = [int(i) for i in vali_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc2 = roc_auc_score(t_label,t_score)
					score = model1.predict(test_ohe)
					t_label = [int(i) for i in test_label]
					t_score = [i[1] for i in score]
					print(roc_auc_score(t_label,t_score))
					roc3 = roc_auc_score(t_label,t_score)
					out1.write(locs+"loc"+"\t"+sys.argv[4]+"\t"+str(co)+"\t"+str(dr)+"\t"+str(do)+"\t"+str(i)+'\t'+str(roc1)+'\t'+str(roc2)+'\t'+str(roc3)+'\n')
					out1.close()
else:
	print("---------------------------------Wrong Input---------------------------")
