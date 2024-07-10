#!python
import tensorflow as tf
import sys
from functions import read_data, one_hot_label, roc_file
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from sklearn.metrics import roc_auc_score
from keras.regularizers import l1, l2
from shutil import copyfile
####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation, Attention
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

print(sys.argv[1]+'\n'+sys.argv[2]+'\n'+sys.argv[3])
pro = sys.argv[1]
locs = sys.argv[2]

model1 = model_from_json(open('../2.TrainValidTest/model.'+locs+'.alldata.json').read())
model1.load_weights('../2.TrainValidTest/model'+locs+'.alldata.h5')

max_len = 2000
cutoff = float(sys.argv[3])
##human.dataset
all_data, all_line = [], []
with open("../7.eachpro/inuloc.dataset.scores") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if len(data[4]) < max_len:
			if data[0] == pro:
				for i in range(len(data[4])):
					seque = data[2]
					all_line.append(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+str(i+1))
					all_data.append(data[4][:(i+1)])


train_data = all_data
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score = model1.predict(train_ohe)

out = open("../7.eachpro/"+locs+"/"+pro+".CPT.scores", "w")
for i in range(len(all_line)):
	out.write(all_line[i]+'\t'+str(score[i][1])+'\n')


out.close()
i
##############calculate delta and select high region##############
pro2pos2score = {}
pro2pos2trunc = {}
pro2pos2delta = {}
with open("../7.eachpro/"+locs+"/"+pro+".CPT.scores") as f:
	#f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] not in pro2pos2score:
			pro2pos2score[data[0]] = {}
		pro2pos2score[data[0]][int(data[4])] = float(data[5])


def calc_trunc(p2p2s, uniid, w=16):
	pos = p2p2s[uniid].keys()
	score = p2p2s[uniid].values()
	trunc, delta = {}, {}
	trunc[uniid] = {}
	delta[uniid] = {}
	for i in pos:
		lowlimit = int(max(1,i-(w/2)))
		higlimit = int(min(len(pos)+1,i+(w/2)+1))
		res = [p2p2s[uniid][j] for j in range(lowlimit, higlimit)]
		Strunc = sum(res)/(higlimit-lowlimit)
		trunc[uniid][i] = Strunc
	for i in pos:
		lowlimit = int(max(1,i-(w/2)))
		higlimit = int(min(len(pos)+1,i+(w/2)+1))
		delta[uniid][i] = trunc[uniid][higlimit-1] - trunc[uniid][lowlimit]
	return(trunc, delta)


for pro in pro2pos2score:
	tmp_trunc, tmp_delta = calc_trunc(pro2pos2score, pro, w=16)
	pro2pos2trunc = {**pro2pos2trunc, **tmp_trunc}
	pro2pos2delta = {**pro2pos2delta, **tmp_delta}


out = open("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores", "w")
with open("../7.eachpro/"+locs+"/"+pro+".CPT.scores") as f:
	#out.write(f.readline().strip('\n')+'\tStrunc\tDelta\n')
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		out.write(line+'\t'+str(pro2pos2trunc[data[0]][int(data[4])])+'\t'+str(pro2pos2delta[data[0]][int(data[4])])+'\n')

out.close()

copyfile("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores", "../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores.tmp")
###########write a circle to find all regions
allvalues, alllines = [], []
finalscore = 0
with open("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		allvalues.append(float(data[7]))
		alllines.append('\t'.join(data[:5]))
		finalscore = float(data[5])

tmp_max = 0
while sorted(allvalues, reverse = True)[0] > 0.1:
	pos = []
	for i in range(len(allvalues)):
		if allvalues[i] > 0.1:
			pos.append(i+tmp_max)
	
	if (len(seque) - (max(pos)+8) < 20):
		break
	tmp_max = max(pos)
	all_line, all_data = [], []
	for i in range(max(pos)+8, len(seque)):
		all_line.append(alllines[i])
		all_data.append(seque[max(pos)+8:i])
	
	train_data = all_data
	train_data = one_hot_label(train_data)
	train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
	train_data = np.array(train_data)
	train_ohe = np_utils.to_categorical(train_data, num_classes=22)
	score = model1.predict(train_ohe)
	finalscore = score[-1][1]
	out = open("../7.eachpro/"+locs+"/"+pro+".CPT.scores.tmp", "w")
	for i in range(len(all_line)):
		out.write(all_line[i]+'\t'+str(score[i][1])+'\n')
	
	out.close()
	
	pro2pos2score = {}
	pro2pos2trunc = {}
	pro2pos2delta = {}
	with open("../7.eachpro/"+locs+"/"+pro+".CPT.scores.tmp") as f:
		#f.readline()
		for line in f:
			data = line.strip('\n').split('\t')
			if data[0] not in pro2pos2score:
				pro2pos2score[data[0]] = {}
			pro2pos2score[data[0]][int(data[4])] = float(data[5])
	
	def calc_trunc_circ(p2p2s, uniid, w=16):
		pos = p2p2s[uniid].keys()
		score = p2p2s[uniid].values()
		trunc, delta = {}, {}
		trunc[uniid] = {}
		delta[uniid] = {}
		for i in pos:
			lowlimit = int(max(min(pro2pos2score[uniid].keys()),i-(w/2)))
			higlimit = int(min(max(pos)+1,i+(w/2)+1))
			res = [p2p2s[uniid][j] for j in range(lowlimit, higlimit)]
			Strunc = sum(res)/(higlimit-lowlimit)
			trunc[uniid][i] = Strunc
		for i in pos:
			lowlimit = int(max(min(pro2pos2score[uniid].keys()),i-(w/2)))
			higlimit = int(min(max(pos)+1,i+(w/2)+1))
			delta[uniid][i] = trunc[uniid][higlimit-1] - trunc[uniid][lowlimit]
		return(trunc, delta)
	
	
	for pro in pro2pos2score:
		tmp_trunc, tmp_delta = calc_trunc_circ(pro2pos2score, pro, w=16)
		pro2pos2trunc = {**pro2pos2trunc, **tmp_trunc}
		pro2pos2delta = {**pro2pos2delta, **tmp_delta}
	
	out = open("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores.regions", "w")
	allvalues = []
	with open("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores.tmp") as f:
		for line in f:
			line = line.strip('\n')
			data = line.split('\t')
			if data[0] in pro2pos2trunc:
				if int(data[4]) in pro2pos2trunc[data[0]]:
					data[-2] = str(pro2pos2trunc[data[0]][int(data[4])])
					data[-1] = str(pro2pos2delta[data[0]][int(data[4])])
					allvalues.append(pro2pos2delta[data[0]][int(data[4])])
					out.write('\t'.join(data)+'\n')
				else:
					out.write(line+'\n')
	out.close()
	copyfile("../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores.regions", "../7.eachpro/"+locs+"/"+pro+".CPT.delta.scores.tmp")











