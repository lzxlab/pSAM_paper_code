#!python

import os
import re

filepath = '../0.data'
filelist = os.listdir(filepath)
path_files = []

regex1 = re.compile('Centroso', re.IGNORECASE)
regex2 = re.compile('Cytoske', re.IGNORECASE)
regex3 = re.compile('Cytoso', re.IGNORECASE)
regex4 = re.compile('Endoso', re.IGNORECASE)
regex5 = re.compile('Endoplasmic', re.IGNORECASE)
regex6 = re.compile('Golgi', re.IGNORECASE)
regex7 = re.compile('Lysoso', re.IGNORECASE)
regex8 = re.compile('Membr', re.IGNORECASE)
regex9 = re.compile('Mitocho', re.IGNORECASE)
regex10 = re.compile('Nucle', re.IGNORECASE)
regex11 = re.compile('Peroxi', re.IGNORECASE)
regex12 = re.compile('Secre', re.IGNORECASE)

for file in filelist:
	filename = file.split('.')
	if filename[-1] == 'tsv':
		path_files.append(file)
		
for tsv in path_files:
	organism = tsv.split('_')[0]
	with open(os.path.join(filepath, tsv)) as f:
		for line in f:
			line = line.strip('\n')
			data = line.split('\t')
			gene = data[1]
			loc = data[3]
			score = data[4]
			if float(score) > 4.0:
				if re.match(regex1, loc) != None:
					out1 = open("../0.data/0.3.1Centrosome_predict_data.txt", "a")
					out1.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex2, loc) != None:
					out2 = open("../0.data/0.3.2Cytoskeleton_predict_data.txt", "a")
					out2.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex3, loc) != None:
					out3 = open("../0.data/0.3.3Cytosol_predict_data.txt", "a")
					out3.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex4, loc) != None:
					out4 = open("../0.data/0.3.4Endosome_predict_data.txt", "a")
					out4.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex5, loc) != None:
					out5 = open("../0.data/0.3.5ER_predict_data.txt", "a")
					out5.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex6, loc) != None:
					out6 = open("../0.data/0.3.6Golgi_predict_data.txt", "a")
					out6.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex7, loc) != None:
					out7 = open("../0.data/0.3.7Lysosome_predict_data.txt", "a")
					out7.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex8, loc) != None:
					out8 = open("../0.data/0.3.8Membrane_predict_data.txt", "a")
					out8.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex9, loc) != None:
					out9 = open("../0.data/0.3.9Mitochondrion_predict_data.txt", "a")
					out9.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex10, loc) != None:
					out10 = open("../0.data/0.3.10Nucleus_predict_data.txt", "a")
					out10.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex11, loc) !=None:
					out11 = open("../0.data/0.3.11Peroxisome_predict_data.txt", "a")
					out11.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')
				if re.match(regex12, loc) != None:
					out12 = open("../0.data/0.3.12Secreted_predict_data.txt", "a")
					out12.write(str(gene)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(score)+'\n')

out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
out6.close()
out7.close()
out8.close()
out9.close()
out10.close()
out11.close()
out12.close()
