#!python

import re

filename1 = "../0.data/uniprot-proteome.fasta"
filename2 = "../0.data/uniprot.subloc.wholedataset.redeal.txt"

regex = re.compile('OS=.+OX')
seq = ''
file1 = {}
model_organism = ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Drosophila melanogaster", "Saccharomyces cerevisiae", "Arabidopsis thaliana", "Caenorhabditis elegans"]

with open(filename1) as f1:
	for line1 in f1:
		if line1.startswith('>'):
			name1 = line1.split("|")[1]
			info = line1.split("|")[2]
			org = regex.findall(info)
			for i in org:
				organism = i[3:-3].split()
				organism = organism[0]+" "+organism[1]
			file1[name1] = organism

output = open("../0.data/0.1uniprot.subloc.model_organism.txt", "w")
with open(filename2) as f2:
	for line2 in f2:
		line2 = line2.strip('\n')
		data = line2.split('\t')
		name2 = data[0]
		loc = data[1]
		seq = data[2]
		for name1, organism in file1.items():
			if name1 == name2:
				#if organism in model_organism:
				output.write(str(name2)+'\t'+str(organism)+'\t'+str(loc)+'\t'+str(seq)+'\n')
output.close

