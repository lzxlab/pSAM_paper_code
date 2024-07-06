#!python
import sys
AA_list_k=["G","A","L","I","V","P","F","M","W","S","Q","T","C","N","Y","D","E","K","R","H","-"]

regions = sys.argv[1]

uni2seq = {}
with open("../4.panMut/hs.protein.score.Nucleus") as f:
    for line in f:
        data = line.strip('\n').split('\t')
        uni2seq[data[0]] = data[2]

out = open(regions+".fasta", "w")
with open(regions) as f:
    for line in f:
        data = line.strip('\n').split('\t')
        tmpdata = data[1].split(';')
        for i in tmpdata:
            seq = uni2seq[data[0]]
            start = int(i.split('...')[0]) - 1
            ends = int(i.split('...')[1])
            out.write(">"+data[0]+"_"+i+'\n')
            tmpseq = list(seq[start:ends])
            resseq = []
            for i in tmpseq:
                if i in AA_list_k:
                    resseq.append(i)
                else:
                    resseq.append("A")
            out.write("".join(resseq)+'\n')

out.close()










