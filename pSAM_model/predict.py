#!python
import tensorflow as tf
import sys
import os
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

##please modify this path before run
path_of_pSAM="/home/zhengyq/data/program/pSAM"

##load parameter and 
inputfasta=sys.argv[1]
taskDir = sys.argv[2]
cuthigh, cutmed, cutlow = 0.736568511, 0.475949377, 0.327701151
model1 = model_from_json(open(path_of_pSAM+'/model.alldata.json').read())
model1.load_weights(path_of_pSAM+'/model.alldata.h5')

def readFa(fa):
	with open(fa,'r') as FA:
		seqName,seq='',''
		while 1:
			line=FA.readline()
			line=line.strip('\n')
			if (line.startswith('>') or not line) and seqName:
				yield((seqName,seq))
			if line.startswith('>'):
				seqName = line[1:]
				seq=''
			else:
				seq+=line
			if not line:break

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


fastafile = readFa(inputfasta)
max_len = 2000
for line in fastafile:
	data = list(line)
	if len(data[1]) > max_len:
		out = open(taskDir+"/merged.data.txt", "w")
		out.write("ERROR: At least one of your input FASTA sequence contains moren than 2000 AAs!")
		out.close()
		sys.exit(0)

##human.dataset
#fastafile = readFa(taskDir+"/inputFile.fasta")
for line in fastafile:
	all_data, all_line = [], []
	data = list(line)
	tmpid = data[0].split(' ')[0]
	if len(data[1]) < max_len:
		for i in range(len(data[1])):
			seque = data[1]
			all_line.append(data[0].split(' ')[0]+'\t'+str(i+1))
			all_data.append(data[1][:(i+1)])
	train_data = all_data
	train_data = one_hot_label(train_data)
	train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
	train_data = np.array(train_data)
	# One hot encoding of sequences
	train_ohe = np_utils.to_categorical(train_data, num_classes=22)
	score = model1.predict(train_ohe)
	
	out = open( taskDir+"/"+tmpid+".res.CPT.scores", "w")
	for i in range(len(all_line)):
		out.write(all_line[i]+'\t'+str(score[i][1])+'\n')
	out.close()
	##############calculate delta and select high region##############
	pro2pos2score = {}
	pro2pos2trunc = {}
	pro2pos2delta = {}
	with open( taskDir+"/"+tmpid+".res.CPT.scores") as f:
		#f.readline()
		for line in f:
			data = line.strip('\n').split('\t')
			if data[0] not in pro2pos2score:
				pro2pos2score[data[0]] = {}
			pro2pos2score[data[0]][int(data[1])] = float(data[2])
	
	for pro in pro2pos2score:
		tmp_trunc, tmp_delta = calc_trunc(pro2pos2score, pro, w=16)
		pro2pos2trunc = {**pro2pos2trunc, **tmp_trunc}
		pro2pos2delta = {**pro2pos2delta, **tmp_delta}
	
	
	out = open( taskDir+"/"+tmpid+".res.CPT.delta.scores", "w")
	with open( taskDir+"/"+tmpid+".res.CPT.scores") as f:
		#out.write(f.readline().strip('\n')+'\tStrunc\tDelta\n')
		for line in f:
			line = line.strip('\n')
			data = line.split('\t')
			out.write(line+'\t'+str(pro2pos2trunc[data[0]][int(data[1])])+'\t'+str(pro2pos2delta[data[0]][int(data[1])])+'\n')
	
	out.close()
	
	copyfile( taskDir+"/"+tmpid+".res.CPT.delta.scores",  taskDir+"/"+tmpid+".res.CPT.delta.scores.tmp")
	###########write a circle to find all regions
	allvalues, alllines = [], []
	finalscore = 0
	with open( taskDir+"/"+tmpid+".res.CPT.delta.scores") as f:
		for line in f:
			data = line.strip('\n').split('\t')
			allvalues.append(float(data[4]))
			alllines.append('\t'.join(data[:2]))
			finalscore = float(data[2])
	
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
		out = open( taskDir+"/"+tmpid+".res.CPT.scores.tmp", "w")
		for i in range(len(all_line)):
			out.write(all_line[i]+'\t'+str(score[i][1])+'\n')
		
		out.close()
		
		pro2pos2score = {}
		pro2pos2trunc = {}
		pro2pos2delta = {}
		with open( taskDir+"/"+tmpid+".res.CPT.scores.tmp") as f:
			#f.readline()
			for line in f:
				data = line.strip('\n').split('\t')
				if data[0] not in pro2pos2score:
					pro2pos2score[data[0]] = {}
				pro2pos2score[data[0]][int(data[1])] = float(data[2])
		
		for pro in pro2pos2score:
			tmp_trunc, tmp_delta = calc_trunc_circ(pro2pos2score, pro, w=16)
			pro2pos2trunc = {**pro2pos2trunc, **tmp_trunc}
			pro2pos2delta = {**pro2pos2delta, **tmp_delta}
		
		out = open( taskDir+"/"+tmpid+".res.CPT.delta.scores.regions", "w")
		allvalues = []
		with open( taskDir+"/"+tmpid+".res.CPT.delta.scores.tmp") as f:
			for line in f:
				line = line.strip('\n')
				data = line.split('\t')
				if data[0] in pro2pos2trunc:
					if int(data[1]) in pro2pos2trunc[data[0]]:
						data[-2] = str(pro2pos2trunc[data[0]][int(data[1])])
						data[-1] = str(pro2pos2delta[data[0]][int(data[1])])
						allvalues.append(pro2pos2delta[data[0]][int(data[1])])
						out.write('\t'.join(data)+'\n')
					else:
						out.write(line+'\n')
		out.close()
		copyfile( taskDir+"/"+tmpid+".res.CPT.delta.scores.regions",  taskDir+"/"+tmpid+".res.CPT.delta.scores.tmp")


import re

def disorder(disorderfile):
	disorder_score = []
	with open(disorderfile) as f:
		for line in f:
			line = line.strip('\n')
			data = re.split(r' ',line)
			if re.search('^#',line):
				continue
			else:
				disorder_score.append(data[-1])
	#print(disorder_score)
	f.close()
	return disorder_score

def hydropathy_charge_polar(fastaFile):
	with open(fastaFile) as f:
		f.readline()
		AA_info = f.readline().strip('\n')
	if 'U' in AA_info :
		aa_info = AA_info.replace('U','A')
	elif 'X' in AA_info:
		aa_info = AA_info.replace('X','A')
	else:
		aa_info = AA_info
	
	hydropathy = []
	polar = []
	charge = []
	amino_hydropathy = {'A':1.8,'R':-4.5, 'N':-3.5, 'D':-3.5,'C':2.5,'Q':-3.5, 'E':-3.5, 'G':-0.4,'H':-3.2, 'I':4.5,'L':3.8,'K':-3.9, 'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
	polar_AA = ('G','F','N','Q','S','T','Y','H','K','R','D','E')
	positive_charge_AA = ('R','H','K')
	negative_charge_AA = ('D','E')

	for i in aa_info:
		if i in amino_hydropathy:
			hydropathy.append(amino_hydropathy[i])
		else:
			hydropathy.append(0)
	#hydropathy.append('hydropathy')
	for j in aa_info:
		if j in polar_AA:
			polar.append(1)
		else:
			polar.append(0)
	for i in aa_info:
		if i in positive_charge_AA:
			charge.append(1)
		elif i in negative_charge_AA:
			charge.append(-1)
		else:
			charge.append(0)

	return hydropathy,polar,charge

def e_b_rsa(eb_rsa_file):
	e_b = []
	rsa = []		##rsa: Relative Surface Accessibility
	ABC = []
	with open(eb_rsa_file) as f:
		for line in f:
			line = line.strip('\n')
			data = '\t'.join(line.split())
			data = data.split('\t')
			
			if re.search('^#',line):
				continue
			else:
				while True:
					if len(data) == 10 :
						break
					else:
						data.insert(0,'')
				#print(data)
				rsa_info = data[4]
				rsa.append(rsa_info)
				e_b_info = data [0]
				if e_b_info == 'E':								   
					e_b.append(1)
				else:
					e_b.append(0)
	#return e_b,rsa
				abc = [data[7],data[8],data[9]]
				if abc.index(max(abc)) == 0 :
					ABC.append('A')
				elif abc.index(max(abc)) == 1 :
					ABC.append('B')
				else:
					ABC.append('C')
	return e_b,rsa,ABC


fastafile = readFa(taskDir+"/inputFile.fasta")
for line in fastafile:
	data = list(line)
	tmpid = data[0].split(' ')[0]
	seque = data[1]
	with open(taskDir+"/"+tmpid+'.res.CPT.fasta','w',newline='\n') as f:
		f.write('>'+tmpid+'\n'+seque+'\n')
	infile = taskDir+"/"+tmpid+'.res.CPT.fasta'
	outdir = taskDir+"/"
	os.system("export IUPred_PATH='/var/www/html/software/iupred/'\nchmod 777 '"+infile+"'\n/var/www/html/software/iupred/iupred '"+infile+"' long > '"+outdir+tmpid+".res.CPT.iupred.txt'")
	os.system("/var/www/html/software/netsurfp-1.0/netsurfp -i '"+infile+"' -t FASTA  -a -k -o '"+outdir+tmpid+".res.CPT.netsurfp.txt'")
	k_disorder = disorder(outdir+tmpid+".res.CPT.iupred.txt")
	k_exposed_buried,k_surface_acc,k_second_str_ABC = e_b_rsa(outdir+tmpid+".res.CPT.netsurfp.txt")
	k_hydropathy,k_polar,k_charge = hydropathy_charge_polar(infile)
        
	out = open(outdir+tmpid+".res.CPT.otherPred.scores",'w',newline='\n')
	out.write('\t'.join([','.join(list(map(str,k_disorder))),','.join(list(map(str,k_exposed_buried))),','.join(list(map(str,k_surface_acc))),','.join(list(map(str,k_second_str_ABC))),','.join(list(map(str,k_hydropathy))),','.join(list(map(str,k_polar))),','.join(list(map(str,k_charge)))]))
	out.close()
	
	



def list2regions(allpos, allval):
	tmpvalue = [float(i) for i in allval]
	if max(tmpvalue) <= 0.1:
		tmp_res = []
	else:
		sigpos = []
		for i in range(len(allpos)):
			if float(allval[i]) > 0.1:
				sigpos.append(allpos[i])
		sigpos= [int(i) for i in sigpos]
		sigpos = sorted(list(set(sigpos)), reverse = False)
		mins, maxs = [], []
		for i in range(len(sigpos)):
			if i == 0:
				mins.append(sigpos[i])
			else:
				if sigpos[i] - sigpos[i-1] == 1:
					tmpmax = sigpos[i]
				else:
					maxs.append(sigpos[i-1])
					mins.append(sigpos[i])
		maxs.append(sigpos[-1])
		tmp_res = []
		for i in range(len(mins)):
			tmp_res.append(str(mins[i])+'...'+str(maxs[i]))
	return tmp_res


out = open(taskDir+"/merged.data.txt", "w")
fastafile = readFa(taskDir+"/inputFile.fasta")
for line in fastafile:
	data = list(line)
	if len(data[1]) < max_len:
		tmpid = data[0].split(' ')[0]
		seqs = data[1]
		poss, scoress = [], []
		with open( taskDir+"/"+tmpid+".res.CPT.delta.scores") as f:
			for line in f:
				data = line.strip('\n').split('\t')
				finalscore = data[2]
				poss.append(data[1])
				scoress.append(data[3])
		deltascores = []
		
		with open(taskDir+"/"+tmpid+".res.CPT.otherPred.scores") as f:
			otherPredScores = f.readline().strip('\n')
			
		if os.path.exists(taskDir+"/"+tmpid+".res.CPT.delta.scores.regions"):
			with open( taskDir+"/"+tmpid+".res.CPT.delta.scores.regions") as f:
				for line in f:
					data = line.strip('\n').split('\t')
					deltascores.append(data[4])
			regions = list2regions(poss, deltascores)
			if float(finalscore) > cuthigh:
				potent = "High"
			elif float(finalscore) > cutmed:
				potent = "Medium"
			elif float(finalscore) > cutlow:
				potent = "Low"
			else:
				potent = "Nonucleus"
			out.write(tmpid+'\t'+seqs+'\t'+finalscore+'\t'+";".join(regions)+'\t'+','.join(scoress)+'\t'+','.join(deltascores)+'\t'+otherPredScores+'\t'+potent+'\n')
		else:
			with open( taskDir+"/"+tmpid+".res.CPT.delta.scores") as f:
				for line in f:
					data = line.strip('\n').split('\t')
					deltascores.append(data[4])
			regions = list2regions(poss, deltascores)
			if float(finalscore) > cuthigh:
				potent = "High"
			elif float(finalscore) > cutmed:
				potent = "Medium"
			elif float(finalscore) > cutlow:
				potent = "Low"
			else:
				potent = "Nonucleus"
			out.write(tmpid+'\t'+seqs+'\t'+finalscore+'\t'+";".join(regions)+'\t'+','.join(scoress)+'\t'+','.join(deltascores)+'\t'+otherPredScores+'\t'+potent+'\n')

out.close()


def removetmpfiles():
	filePath = taskDir+"/"
	name = os.listdir(filePath)
	for i in name:
		tmppath = filePath+'{}'
		path = tmppath.format(i)
		print(path)
		if 'res.CPT' in i:
			os.remove(path)

removetmpfiles()

import time 

with open(taskDir+"/runInfo.txt","a") as f:
	f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

with open(taskDir+"/runInfo.txt") as f:
	inputEmail = f.readline().strip('\n')
if inputEmail != '':
	os.system('echo "Dear user,\n\nThanks for your interest in iNuLoC!\n\nYour task has compeleted, and you can access the result with the following link(s):\n\nhttp://inuloc.omicsbio.info/action.php?taskid='+taskDir+'\n\nYours sincerely,\nZexian Liu Ph.D. Associate Professor,\nSun Yat-sen University Cancer Center\n\nBuilding 2#20F, 651 Dongfeng East Road,\n\nGuangzhou 510060, P. R. China\nTel/Fax: +86-20-87342025\nPersonal website: http://lzx.cool" | mail -s "pNuLoC prediction result" '+inputEmail)