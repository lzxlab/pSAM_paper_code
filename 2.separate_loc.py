#!python

filename = "../0.data/0.1uniprot.subloc.model_organism.txt"
loc = []
with open(filename) as f:
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		loc_info = data[2].split(',')
		for i in loc_info:
			loc.append(str(i))
loc = list(set(loc))

for i in loc:
	with open(filename) as f:
		for line in f:
			line = line.strip('\n')
			data = line.split('\t')
			protein = data[0]
			organism = data[1]
			localization = data[2]
			local = localization.split(',')
			seq = data[3]
			if i in local:
				output1 = open("../0.data/"+"0.2"+str(i)+"_"+"positive_data"+".txt", "a")
				output1.write(str(protein)+'\t'+str(organism)+'\t'+str(localization)+'\t'+str(seq)+'\n')
			else:
				output2 = open("../0.data/"+"0.2"+str(i)+"_"+"negative_data"+".txt", "a")
				output2.write(str(protein)+'\t'+str(organism)+'\t'+str(localization)+'\t'+str(seq)+'\n')
output1.close()
output2.close()
