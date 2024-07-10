#!python

with open("../tmp/seq.txt") as f:
    for line in f:
        data = line.strip('\n').split('\t')
        seqid = data[0]
        seqpep = data[1]

allclass = ["Origin"]
aachange = ["No"]
allseqs = [seqpep]

for i in range(1, 28):
    allclass.append("Single-mut")
    aachange.append(seqpep[i]+str(i)+"-"+"A")
    tmpseq = list(seqpep)
    tmpseq[i] = "A"
    allseqs.append("".join(tmpseq))

for i in range(1, 27):
    for j in range(i+1, 28):
        allclass.append("Pair-mut")
        aachange.append(seqpep[i]+str(i)+"-"+"A; "+seqpep[j]+str(j)+"-"+"A")
        tmpseq = list(seqpep)
        tmpseq[i] = "A"
        tmpseq[j] = "A"
        allseqs.append("".join(tmpseq))

tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["H", "R", "K"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("Positive-mut")
aachange.append("HRK-A")
allseqs.append("".join(newseq))

tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["D", "E"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("Negative-mut")
aachange.append("DE-A")
allseqs.append("".join(newseq))


tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["F", "A", "L", "M", "I", "W", "P", "V"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("Hydrophobic-mut")
aachange.append("FALMIWPV-A")
allseqs.append("".join(newseq))


tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i not in ["F", "A", "L", "M", "I", "W", "P", "V"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("Hydrophilic-mut")
aachange.append("notFALMIWPV-A")
allseqs.append("".join(newseq))


tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["S", "T", "Y"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("STY-mut")
aachange.append("STY-A")
allseqs.append("".join(newseq))

out = open("../tmp/seq.mut.txt", "w")
for i in range(len(allclass)):
    out.write(allclass[i]+'\t'+aachange[i]+'\t'+allseqs[i]+'\n')

out.close()

alllines = []
allseqs = []
with open("../tmp/seq.mut.txt") as f:
    for line in f:
        line = line.strip('\n')
        data = line.split('\t')
        alllines.append(line)
        allseqs.append(data[2])



#!python
import tensorflow as tf
#tf.config.gpu.set_per_process_memory_growth(True)
import sys
from functions import read_data, one_hot_label, roc_file
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2
import re
####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

model1 = model_from_json(open('../2.TrainValidTest/model.Mitochondrion.alldata.json').read())
model1.load_weights('../2.TrainValidTest/modelMitochondrion.alldata.h5')

max_len = 2000

train_data = allseqs
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe, batch_size = 1)


out = open("../tmp/seq.mut.score.txt", "w")
for i in range(len(alllines)):
	out.write(alllines[i]+'\t'+str(score_mut[i][1])+'\n')


out.close()



##################
#!python

with open("../tmp/seq.txt") as f:
    for line in f:
        data = line.strip('\n').split('\t')
        seqid = data[0]
        seqpep = data[1]

allclass = ["Origin"]
aachange = ["No"]
allseqs = [seqpep]

for i in range(1, 28):
    if seqpep[i] in ["H", "R", "K"]:
        allclass.append("Single-mut")
        aachange.append(seqpep[i]+str(i)+"-"+"A")
        tmpseq = list(seqpep)
        tmpseq[i] = "A"
        allseqs.append("".join(tmpseq))

for i in range(1, 28):
    if seqpep[i] in ["H", "R", "K"]:
        allclass.append("Single-mut")
        aachange.append(seqpep[i]+str(i)+"-"+"E")
        tmpseq = list(seqpep)
        tmpseq[i] = "E"
        allseqs.append("".join(tmpseq))

for i in range(1, 27):
    for j in range(i+1, 28):
        if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"]:
            allclass.append("Pair-mut")
            aachange.append(seqpep[i]+str(i)+"-"+"A; "+seqpep[j]+str(j)+"-"+"A")
            tmpseq = list(seqpep)
            tmpseq[i] = "A"
            tmpseq[j] = "A"
            allseqs.append("".join(tmpseq))

for i in range(1, 27):
    for j in range(i+1, 28):
        if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"]:
            allclass.append("Pair-mut")
            aachange.append(seqpep[i]+str(i)+"-"+"E; "+seqpep[j]+str(j)+"-"+"E")
            tmpseq = list(seqpep)
            tmpseq[i] = "E"
            tmpseq[j] = "E"
            allseqs.append("".join(tmpseq))

for i in range(1, 26):
    for j in range(i+1, 27):
        for k in range(j+1, 28):
            if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"] and seqpep[k] in ["H", "R", "K"]:
                allclass.append("Triple-mut")
                aachange.append(seqpep[i]+str(i)+"-"+"A; "+seqpep[j]+str(j)+"-"+"A; "+seqpep[k]+str(k)+"-"+"A")
                tmpseq = list(seqpep)
                tmpseq[i] = "A"
                tmpseq[j] = "A"
                tmpseq[k] = "A"
                allseqs.append("".join(tmpseq))

for i in range(1, 26):
    for j in range(i+1, 27):
        for k in range(j+1, 28):
            if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"] and seqpep[k] in ["H", "R", "K"]:
                allclass.append("Triple-mut")
                aachange.append(seqpep[i]+str(i)+"-"+"E; "+seqpep[j]+str(j)+"-"+"E; "+seqpep[k]+str(k)+"-"+"E")
                tmpseq = list(seqpep)
                tmpseq[i] = "E"
                tmpseq[j] = "E"
                tmpseq[k] = "E"
                allseqs.append("".join(tmpseq))

for i in range(1, 25):
    for j in range(i+1, 26):
        for k in range(j+1, 27):
            for h in range(k+1, 28):
                if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"] and seqpep[k] in ["H", "R", "K"] and seqpep[h] in ["H", "R", "K"]:
                    allclass.append("Four-mut")
                    aachange.append(seqpep[i]+str(i)+"-"+"A; "+seqpep[j]+str(j)+"-"+"A; "+seqpep[k]+str(k)+"-"+"A; "+seqpep[h]+str(h)+"-"+"A")
                    tmpseq = list(seqpep)
                    tmpseq[i] = "A"
                    tmpseq[j] = "A"
                    tmpseq[k] = "A"
                    tmpseq[h] = "A"
                    allseqs.append("".join(tmpseq))

for i in range(1, 25):
    for j in range(i+1, 26):
        for k in range(j+1, 27):
            for h in range(k+1, 28):
                if seqpep[i] in ["H", "R", "K"] and seqpep[j] in ["H", "R", "K"] and seqpep[k] in ["H", "R", "K"] and seqpep[h] in ["H", "R", "K"]:
                    allclass.append("Four-mut")
                    aachange.append(seqpep[i]+str(i)+"-"+"E; "+seqpep[j]+str(j)+"-"+"E; "+seqpep[k]+str(k)+"-"+"E; "+seqpep[h]+str(h)+"-"+"E")
                    tmpseq = list(seqpep)
                    tmpseq[i] = "E"
                    tmpseq[j] = "E"
                    tmpseq[k] = "E"
                    tmpseq[h] = "E"
                    allseqs.append("".join(tmpseq))


tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["H", "R", "K"]:
        newseq.append("A")
    else:
        newseq.append(i)

allclass.append("Positive-mut")
aachange.append("HRK-A")
allseqs.append("".join(newseq))

tmpseq = list(seqpep)
newseq = []
for i in tmpseq:
    if i in ["H", "R", "K"]:
        newseq.append("E")
    else:
        newseq.append(i)

allclass.append("Positive-mut")
aachange.append("HRK-E")
allseqs.append("".join(newseq))


out = open("../tmp/seq.mut.HRK.txt", "w")
for i in range(len(allclass)):
    out.write(allclass[i]+'\t'+aachange[i]+'\t'+allseqs[i]+'\n')

out.close()

alllines = []
allseqs = []
with open("../tmp/seq.mut.HRK.txt") as f:
    for line in f:
        line = line.strip('\n')
        data = line.split('\t')
        alllines.append(line)
        allseqs.append(data[2])



#!python
import tensorflow as tf
#tf.config.gpu.set_per_process_memory_growth(True)
import sys
from functions import read_data, one_hot_label, roc_file
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2
import re
####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

model1 = model_from_json(open('../2.TrainValidTest/model.Mitochondrion.alldata.json').read())
model1.load_weights('../2.TrainValidTest/modelMitochondrion.alldata.h5')

max_len = 2000

train_data = allseqs
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
score_mut = model1.predict(train_ohe, batch_size = 1)


out = open("../tmp/seq.mut.HRK.score.txt", "w")
for i in range(len(alllines)):
	out.write(alllines[i]+'\t'+str(score_mut[i][1])+'\n')


out.close()







