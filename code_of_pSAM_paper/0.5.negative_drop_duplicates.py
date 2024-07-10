#!python

import os
import re

filepath = '../0.data'
filelist = os.listdir(filepath)
duplicates_files = []
negative_files = []

for file in filelist:
	filename = file.split('_')
	if len(filename) > 1:
		if filename[1] == 'duplicates':
			duplicates_files.append(file)
		if filename[1] == 'negative':
			negative_files.append(file)

for dup in duplicates_files:
	local1 = dup.split('_')[0]
	local1 = local1[6:]
	id_dup = []
	with open(os.path.join(filepath, dup)) as f1:
		for line1 in f1:
			line1 = line1.strip('\n')
			data1 = line1.split('\t')
			id1 = data1[0]
			id_dup.append(id1)
		id_dup = list(set(id_dup))
	for neg in negative_files:
		local2 = neg.split('_')[0]
		local2 = local2[8:]
		if local1 == local2:
			with open(os.path.join(filepath, neg)) as f2:
				for line2 in f2:
					line2 = line2.strip('\n')
					data2 = line2.split('\t')
					id2 = data2[0]
					org = data2[1]
					loc = data2[2]
					seq = data2[3]
					if id2 not in id_dup:
						out = open("../0.data/"+"0.5"+str(local1)+"_negative_drop_predict.txt", "a")
						out.write(str(id2)+'\t'+str(org)+'\t'+str(loc)+'\t'+str(seq)+'\n')
						out.close()
