import os
import re

filepath = '../0.data'
fasta = "../0.data/uniprot-proteome_(taxonomy_Eukaryota+[2759]_)-filtered-reviewed--.fasta"
filelist = os.listdir(filepath)
predict_files = []
negative_files = []
organisms = {"Homo sapiens": "human", "Mus musculus": "mouse", "Rattus norvegicus": "rat", "Drosophila melanogaster": "fly",
             "Saccharomyces cerevisiae": "yeast", "Arabidopsis thaliana": "arabidopsis", "Caenorhabditis elegans": "worm"}

regex1 = re.compile('OS=.+OX')
regex2 = re.compile('GN=.+PE')

for file in filelist:
    filename = file.split('_')
    if len(filename) > 1:
        if filename[1] == 'predict':
            predict_files.append(file)
        if filename[1] == 'negative':
            negative_files.append(file)
print(negative_files)
print(predict_files)

for neg in negative_files:
    local1 = neg.split('_')[0]
    local1 = local1[8:]
    with open(os.path.join(filepath, neg)) as f1:
        for line1 in f1:
            line1 = line1.strip('\n')
            data1 = line1.split('\t')
            id1 = data1[0]
            org1 = data1[1]
            loc1 = data1[2]
            seq1 = data1[3]
            with open(fasta) as f2:
                for line2 in f2:
                    if line2.startswith('>'):
                        id2 = line2.split("|")[1]
                        info = line2.split("|")[2]
                        org = regex1.findall(info)
                        gene1 = regex2.findall(info)
                        for x in org:
                            org2 = x[3:-3]
                        for y in gene1:
                            gene1 = y[3:-3]
                        if id1 == id2:
                            print("111")
                            for pre in predict_files:
                                local2 = pre.split('_')[0]
                                local2 = local2[6:]
                                if local1 == local2:
                                    with open(os.path.join(filepath, pre)) as f3:
                                        for line3 in f3:
                                            line3 = line3.strip('\n')
                                            data3 = line3.split('\t')
                                            gene2 = data3[0]
                                            org3 = data3[1]
                                            loc2 = data3[2]
                                            score = data3[3]
                                            if organisms[org2] == org3:
                                                if gene1 == gene2:
                                                    print("11111")
                                                    out = open("../0.data/"+"0.4"+str(local1)+'_'+"duplicates_data.txt", "a")
                                                    out.write(str(id1)+'\t'+str(gene2)+'\t'+str(org1)+'\t'+str(loc2)+'\t'+str(seq1)+'\t'+str(score)+'\n')
                                                    out.close()


