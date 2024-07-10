#!python
#sys.argv[1] = hs.cpt.delta.method1.sigdelta.sites
import random
import sys
uni2region = {}

with open("../7.eachpro/"+sys.argv[1]) as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		tmpuni = data[0]
		if tmpuni in uni2region:
			uni2region[tmpuni].append(int(data[4]))
		else:
			uni2region[tmpuni] = [int(data[4])]

out = open("../7.eachpro/"+sys.argv[1]+".regions", "w")
for uni in uni2region:
	mins, maxs = [], []
	for i in range(len(sorted(uni2region[uni]))):
		if i == 0:
			mins.append(sorted(uni2region[uni])[i])
		else:
			if sorted(uni2region[uni])[i] - sorted(uni2region[uni])[i-1] == 1:
				tmpmax = sorted(uni2region[uni])[i]
			else:
				maxs.append(sorted(uni2region[uni])[i-1])
				mins.append(sorted(uni2region[uni])[i])
	maxs.append(sorted(uni2region[uni])[-1])
	tmp_res = []
	for i in range(len(mins)):
		tmp_res.append(str(mins[i])+'...'+str(maxs[i]))
	out.write(uni+'\t'+';'.join(tmp_res)+'\n')

out.close()















