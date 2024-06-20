#!python

import os

filepath = "../0.data"
filelist = os.listdir(filepath)
positive_files = []
negative_files = []
negative_predict = []
negative_predict_ppi = []

for file in filelist:
	if file.startswith('0'):
		filename = file.split('_')
		if len(filename) > 1:
			if len(filename)==3 and  filename[1]=="postive" and filename[2]=="data.txt":
				positive_files.append(file)
			if len(filename)==3 and filename[1]=="negative" and filename[2]=="data.txt":
				negative_files.append(file)
			if len(filename)==4 and filename[2]=="drop" and filename[3]=="predict.txt":
				negative_predict.append(file)
			if len(filename)==5 and filename[3]=="predict" and filename[4]=="ppi.txt":
				negative_predict_ppi.append(file) 

positive_files.sort()
negative_files.sort()
negative_predict.sort()
negative_predict_ppi.sort()
print(positive_files)
print(negative_files)
print(negative_predict)
print(negative_predict_ppi)

out = open("statistics.txt", "a")

for file1 in positive_files:
	loc1 = file1.split('_')[0]
	loc1 = loc1[8:]
	print(loc1)
	with open(os.path.join(filepath, file1)) as f1:
		lines1 = f1.readlines()
		count1 = len(lines1)
		for file2 in negative_files:
			loc2 = file2.split('_')[0]
			loc2 = loc2[8:]
			if loc1==loc2:
				print(loc2)
				with open(os.path.join(filepath, file2)) as f2:
					lines2 = f2.readlines()
					count2 = len(lines2)
		for file3 in negative_predict:
			loc3 = file3.split('_')[0]
			loc3 = loc3[6:]
			if loc1==loc3:
				print(loc3)
				with open(os.path.join(filepath, file3)) as f3:
					lines3 = f3.readlines()
					count3 = len(lines3)
		for file4 in negative_predict_ppi:
			loc4 = file4.split('_')[0]
			loc4 = loc4[6:]
			if loc1==loc4:
				print(loc4)
				with open(os.path.join(filepath, file4)) as f4:
					lines4 = f4.readlines()
					count4 = len(lines4)
		if loc1 == 'Nucleus':
			out.write('Nucleus'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Nucleus'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Nucleus'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Nucleus'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'Peroxisome':
			out.write('Peroxisome'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Peroxisome'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Peroxisome'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Peroxisome'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'Secreted':
			out.write('Secreted'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Secreted'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Secreted'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Secreted'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'entrosome':
			out.write('Centrosome'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Centrosome'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Centrosome'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Centrosome'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		#if loc1 == 'ytoskeleton':
		#	out.write('Cytoskeleton'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Cytoskeleton'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Cytoskeleton'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Cytoskeleton'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'ytosolskel':
			out.write('Cytosol'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Cytosol'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Cytosol'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Cytosol'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'ndosome':
			out.write('Endosome'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Endosome'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Endosome'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Endosome'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'R':
			out.write('ER'+'\t'+'positive'+'\t'+str(count1)+'\n'+'ER'+'\t'+'negative'+'\t'+str(count2)+'\n'+'ER'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'ER'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'olgi':
			out.write('Golgi'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Golgi'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Golgi'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Golgi'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'ysosome':
			out.write('Lysosome'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Lysosome'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Lysosome'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Lysosome'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'embrane':
			out.write('Membrane'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Membrane'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Membrane'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Membrane'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
		if loc1 == 'itochondrion':
			out.write('Mitochondrion'+'\t'+'positive'+'\t'+str(count1)+'\n'+'Mitochondrion'+'\t'+'negative'+'\t'+str(count2)+'\n'+'Mitochondrion'+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+'Mitochondrion'+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
'''
		out.write(str(loc1)+'\t'+'positive'+'\t'+str(count1)'\n'+str(loc2)+'\t'+'negative'+'\t'+str(count2)+'\n'+str(loc3)+'\t'+'negative_predict'+'\t'+str(count3)+'\n'+str(loc4)+'\t'+'negative_predict_ppi'+'\t'+str(count4)+'\n')
'''
out.close()
