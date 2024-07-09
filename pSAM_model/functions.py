#!some functions to deal with the origin data
#deal the protein2label file
import numpy as np

def read_data(filename):
	all_data, all_label = [], []
	with open(filename) as f:
		f.readline()
		for line in f:
			line = line.strip('\n')
			data = line.split('\t')
			all_data.append(data[2])
			all_label.append(data[3])
	return all_data, all_label

def one_hot_label(aa_line):
	cor = {}
	AA_list_k=["G","A","L","I","V","P","F","M","W","S","Q","T","C","N","Y","D","E","K","R","H","-"]
	m = 1
	for i in AA_list_k:
		cor[i] = m
		m += 1
	one_hot_data_list = []
	for pep in aa_line:
		one_hot_list = []
		for aa in pep:
			if aa in AA_list_k:
				one_hot_list.append(cor[aa])
			else:
				one_hot_list.append(cor['-'])
		one_hot_data_list.append(one_hot_list)
	#one_hot_data_list = np.array(one_hot_data_list)
	return one_hot_data_list




def roc_file(score, label, num = 1):
	output = open('../3.model/roc_data.'+str(num)+'.txt',"w")
	output.write('score\tlabel\tsp\tsn\tpr\trc\tac\tmcc\n')
	data, pre = [], []
	for i in score:
		data.append(float(i))
	for i in label:
		pre.append(int(i))
	start = 0
	for i in data:
		tp, tn, fp, fn = 0, 0, 0, 0
		for j in range(len(data)):
			if data[j] >= i:
				if pre[j] == 1.:
					tp += 1
				else:
					fp += 1
			else:
				if pre[j] == 1.:
					fn += 1
				else:
					tn += 1
		sp = tn/(tn + fp)
		sn = tp/(tp + fn)
		if (tp + fp) == 0:
			pr = 1
		else:
			pr = tp/(tp + fp)
		rc = sn
		ac = (tp+tn)/(tp+tn+fp+fn)
		low = (tp+fn)*(tp+fp)*(tn+fp)*(tn+fn)
		if low == 0:
			mcc = 1
		else:
			mcc = (tp*tn-fp*fn)/(low**0.5)
		output.write(str(i)+'\t'+str(pre[start])+'\t'+str(sp)+'\t'+str(sn)+'\t'+str(pr)+'\t'+str(rc)+'\t'+str(ac)+'\t'+str(mcc)+'\n')
		start += 1
	output.close()
