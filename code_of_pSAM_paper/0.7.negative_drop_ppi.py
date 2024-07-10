import re
import os

filepath = '../0.data'
fasta = "../0.data/uniprot-proteome_(taxonomy_Eukaryota+[2759]_)-filtered-reviewed--.fasta"
filelist = os.listdir(filepath)
negative_files = []
ppi_files = []

for file in filelist:
	filename = file.split('_')
	if len(filename) > 1:
		if filename[1] == 'negative' and filename[2] == 'drop':
			negative_files.append(file)
		if filename[1] == 'positive' and filename[2] == 'ppi.txt':
			ppi_files.append(file)

print(negative_files)
print(ppi_files)

regex1 = re.compile('OS=.+OX')
regex2 = re.compile('GN=.+PE')

for ppi in ppi_files:
	local = ppi.split('_')[0]
	local = local[5:]
	dict = {}
	with open(os.path.join(filepath, ppi)) as f1:
		for line1 in f1:
			line1 = line1.strip('\n')
			data1 = line1.split('\t')
			string_id = data1[0]
			gene1 = data1[1]
			org1 = data1[2]
			dict[gene1] = org1
	pro_id = []
	with open(fasta) as f2:
		for line2 in f2:
			if line2.startswith('>'):
				id2 = line2.split("|")[1]
				info = line2.split("|")[2]
				org2 = regex1.findall(info)
				gene2 = regex2.findall(info)
				for x in org2:
					org2 = x[3:-3]
				for y in gene2:
					gene2 = y[3:-3]
				for gene1, org1 in dict.items():
					if gene1==gene2 and org1==org2:
						pro_id.append(id2)
	ppi_id = list(set(pro_id))
	for neg in negative_files:
		local2 = neg.split('_')[0]
		local2 = local2[5:]
		if local == local2:
			with open(os.path.join(filepath, neg)) as f3:
				for line3 in f3:
					line3 = line3.strip('\n')
					data3 = line3.split('\t')
					id3 = data3[0]
					org3 = data3[1]
					loc = data3[2]
					seq = data3[3]
					if id3 not in ppi_id:
						out = open("../0.data/"+"0.7"+str(local2)+"_negative_drop_predict_ppi.txt", "a")
						out.write(str(id3)+'\t'+str(org3)+'\t'+str(loc)+'\t'+str(seq)+'\n')
						out.close()
