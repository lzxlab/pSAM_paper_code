import os
import re

file1 = "../0.data/protein.physical.links.v11.5.txt"
file2 = "../0.data/protein.info.v11.5.txt"
file3 = "../0.data/species.v11.5.txt"

protein1 = []
protein2 = []
gene_name = {}
organism = {}

with open(file1) as f1:
	first_line1 = f1.readline()
	for line1 in f1:
		line1 = line1.strip('\n')
		data1 = line1.split()
		pro1 = data1[0]
		pro2 = data1[1]
		p = pro1
		num = p.split('.')[0]
		if num in ["9606", "10090", "10116", "7227", "4932", "3702", "6239"]:
			protein1.append(pro1)
			protein2.append(pro2)
print(protein1)
print(protein2)

with open(file2) as f2:
	first_line2 = f2.readline()
	for line2 in f2:
		line2 = line2.strip('\n')
		data2 = line2.split('\t')
		pro = data2[0]
		gene = data2[1]
		p = pro
		num = p.split('.')[0]
		if num in ["9606", "10090", "10116", "7227", "4932", "3702", "6239"]:
			gene_name[gene] = pro
print(gene_name)

with open(file3) as f3:
	first_line3 = f3.readline()
	for line3 in f3:
		line3 = line3.strip('\n')
		data3 = line3.split('\t')
		num = data3[0]
		org = data3[3]
		if num in ["9606", "10090", "10116", "7227", "4932", "3702", "6239"]:
			organism[num] = org
print(organism)

filepath = "../0.data"
filelist = os.listdir(filepath)
fasta = "../0.data/uniprot-proteome_(taxonomy_Eukaryota+[2759]_)-filtered-reviewed--.fasta"
positive_files = []
regex = re.compile('GN=.+PE')

for file in filelist:
	filename = file.split('_')
	if len(filename) > 1:
		if filename[1] == 'postive':
			positive_files.append(file)
print(positive_files)

for pos in positive_files:
	local = pos.split('_')[0]
	local = local[7:]
	with open(os.path.join(filepath, pos)) as f4:
		for line4 in f4:
			line4 = line4.strip('\n')
			data4 = line4.split('\t')
			pro_id1 = data4[0]
			org = data4[1]
			loc = data4[2]
			seq = data4[3]
			with open(fasta) as f5:
				for line5 in f5:
					if line5.startswith('>'):
						pro_id2 = line5.split('|')[1]
						info = line5.split('|')[2]
						gene = regex.findall(info)
						for x in gene:
							gene = x[3:-3]
						if pro_id1 == pro_id2:
							if gene in gene_name.keys():
								p1 = gene_name[gene]
								print(p1)
								p2 = p1
								org_num = p2.split('.')[0]
								org2 = organism[org_num]
								print(org2)
								if org == org2:
									if p1 in protein1:
										index = protein1.index(p1)
										ppi = protein2[index]
										for x, y in gene_name.items():
											if y == ppi:
												g = x
										out = open("../0.data/"+"0.6"+str(local)+"_positive_ppi.txt", "a")
										out.write(str(ppi)+'\t'+str(g)+'\t'+str(org)+'\n')
										out.close()
									if p1 in protein2:
										index = protein2.index(p1)
										ppi = protein1[index]
										for x, y in gene_name.items():
											if y == ppi:
												g = x
										out = open("../0.data/"+"0.6"+str(local)+"_positive_ppi.txt", "a")
										out.write(str(ppi)+'\t'+str(g)+'\t'+str(org)+'\n')
										out.close()

