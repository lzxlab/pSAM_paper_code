import re

file1 = "../0.data/protein.physical.links.v11.5.txt"
protein1 = []
protein2 = []
scores = []
with open(file1) as f1:
	for line1 in f1:
		line1 = line1.strip('\n')
		data1 = line1.split()
		pro1 = data1[0]
		pro2 = data1[1]
		score = data1[2]
		num = pro1.split('.')[0]
		if num in ["9606", "10090", "10116", "7227", "4932", "3702", "6239"]:
			protein1.append(pro1)
			protein2.append(pro2)
			scores.append(score)
print(len(protein1))
print(len(protein2))
print(len(scores))

file2 = file2 = "../0.data/protein.info.v11.5.txt"
gene_pro = {}
with open(file2) as f2:
	first_line2 = f2.readline()
	for line2 in f2:
		line2 = line2.strip('\n')
		data2 = line2.split('\t')
		pro = data2[0]
		gene = data2[1]
		num = pro.split('.')[0]
		if num in ["9606", "10090", "10116", "7227", "4932", "3702", "6239"]:
			gene_pro[pro] = gene
print(gene_pro)

file3 = "../0.data/species.v11.5.txt"
organism = {}
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

fasta = "../0.data/uniprot-proteome_(taxonomy_Eukaryota+[2759]_)-filtered-reviewed--.fasta"
pro_id = {}
with open(fasta) as f4:
	for line4 in f4:
		if line4.startswith('>'):
			print(line4)
			id_f = line4.split("|")[1]
			gene_f = line4.split("GN=")
			if len(gene_f) > 1:
				gene_f = gene_f[1].split(" ")[0]
				pro_id[gene_f] = id_f
print(pro_id)

for p1, p2 in zip(protein1, protein2):
	org1 = p1.split('.')[0]
	p1_org = organism[org1]
	p1_gene =gene_pro[p1]
	org2 = p2.split('.')[0]
	p2_org = organism[org2]
	p2_gene = gene_pro[p2]
	index = protein1.index(p1)
	score = scores[index]
	if p1_gene in pro_id.keys() and p2_gene in pro_id.keys():
		p1_id = pro_id[p1_gene]
		p2_id = pro_id[p2_gene]
		out = open('../0.data/0.8.ppi.txt', 'a')
		out.write(str(p1_org)+'\t'+str(p1)+'\t'+str(p1_gene)+'\t'+str(p1_id)+'\t'+str(p2_org)+'\t'+str(p2)+'\t'+str(p2_gene)+'\t'+str(p2_id)+'\t'+str(score)+'\n')
out.close()
