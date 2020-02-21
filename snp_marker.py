#!/usr/bin/env python3
#Zimpel et al., 2020. "Global distribution and evolution of Mycobacterium bovis lineages".
#Author: Naila C. Soler Camargo

#The snp_marker.py algorithm takes as input two files:
#i) A file in a txt format containing the name of the strain and the associated cluster separated by a tab.
#ii) A multifasta file containing the snps of each strain.
#Expected output: a table containing the cluster, snp and position of the snp in the fasta file.

import glob
from Bio import SeqIO

#output
output = open('snp_markers.xlsx', 'w')

#initializing lineage database

lb1 = []
lb2 = []
lb3 = []
lb4 = []
lb5 = []
lb6 = []
lb7 = []

a, b, c, d, e, f, g = 0, 0, 0, 0, 0, 0, 0

for file in glob.glob('*.txt'):

	reading = open(file, 'r')
	reading = reading.read()
	reading = reading.split('\n')
	

for file in glob.glob('*.fasta'):

	for seq_record in SeqIO.parse(file, "fasta"):

		for line in reading:

			if seq_record.id in line:

				column = line.split('\t')
				
				if column[1] == '1':

					lb1.append('')
					lb1[a] = seq_record.id
					a+=1
					lb1.append('')
					lb1[a] = str(seq_record.seq)
					a+=1

				if column[1] == '2':

					lb2.append('')
					lb2[b] = seq_record.id
					b+=1
					lb2.append('')
					lb2[b] = str(seq_record.seq)
					b+=1

				if column[1] == '3':

					lb3.append('')
					lb3[c] = seq_record.id
					c+=1
					lb3.append('')
					lb3[c] = str(seq_record.seq)
					c+=1

				if column[1] == '4':

					lb4.append('')
					lb4[d] = seq_record.id
					d+=1
					lb4.append('')
					lb4[d] = str(seq_record.seq)
					d+=1


				if column [1] == '5':

					lb5.append('')
					lb5[e] = seq_record.id
					e+=1
					lb5.append('')
					lb5[e] = str(seq_record.seq)
					e+=1

				if column[1] == '6':
					
					lb6.append('')
					lb6[f] = seq_record.id
					f+=1
					lb6.append('')
					lb6[f] = str(seq_record.seq)
					f+=1

				if column[1] == '7':

					lb7.append('')
					lb7[g] = seq_record.id
					g+=1
					lb7.append('')
					lb7[g] = str(seq_record.seq)
					g+=1


#searching for snp markers by position

num_snps = len(lb1[1])
p = 0

while p < num_snps:
	
	#LB1

	l = 1
	current = []
	a = 0
	while l < len(lb1):

		seq = lb1[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb1)/2): 			

		flag = 0
		for file in glob.glob('*.fasta'):
			for seq_record in SeqIO.parse(file, "fasta"):

				if seq_record.id not in lb1:
				
					sequence = str(seq_record.seq)

					if sequence[p] == current[0]:

						flag = 1

			if flag == 0:

				line = 'cluster1' + '\t' + current[0] + '\t' + str(p+1) + '\n'
				output.writelines(line)

	
	#LB2
	
	l = 1
	current = []
	a = 0
	while l < len(lb2):

		seq = lb2[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb2)/2): 			

		flag = 0
		for file in glob.glob('*.fasta'):

			for seq_record in SeqIO.parse(file, "fasta"):

				if seq_record.id not in lb2:
				
					sequence = str(seq_record.seq)

					if sequence[p] == current[0]:

						flag = 1

			if flag == 0:

				line = 'cluster2' + '\t' + current[0] + '\t' + str(p+1) + '\n'
				output.writelines(line) 

	
	#LB3
	
	l = 1
	current = []
	a = 0
	while l < len(lb3):

		seq = lb3[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb3)/2):

		flag = 0
		for file in glob.glob('*.fasta'):

			for seq_record in SeqIO.parse(file, "fasta"):

				if seq_record.id not in lb3:
				
					sequence = str(seq_record.seq)

					if sequence[p] == current[0]:

						flag = 1

			if flag == 0:

				line = 'cluster3' + '\t' + current[0] + '\t' + str(p+1) + '\n'
				output.writelines(line) 

	
	#LB4
	
	l = 1
	current = []
	a = 0
	while l < len(lb4):

		seq = lb4[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb4)/2): 

		flag = 0
		for file in glob.glob('*.fasta'):

			for seq_record in SeqIO.parse(file, "fasta"):

				if seq_record.id not in lb4:
				
					sequence = str(seq_record.seq)

					if sequence[p] == current[0]:

						flag = 1

			if flag == 0:

				line = 'cluster4' + '\t' + current[0] + '\t' + str(p+1) + '\n' 
				output.writelines(line)

	
	#LB5
	
	l = 1
	current = []
	a = 0
	while l < len(lb5):

		seq = lb5[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb5)/2):

		flag = 0
		for file in glob.glob('*.fasta'):

			for seq_record in SeqIO.parse(file, "fasta"):

				if seq_record.id not in lb5:
				
					sequence = str(seq_record.seq)

					if sequence[p] == current[0]:

						flag = 1

			if flag == 0:

				line = 'cluster5' + '\t' + current[0] + '\t' + str(p+1) + '\n' 
				output.writelines(line)

	#LB6

	l = 1
	current = []
	a = 0
	while l < len(lb6):
        
		seq = lb6[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb6)/2):
        
		flag = 0
		for file in glob.glob('*.fasta'):
            
			for seq_record in SeqIO.parse(file, "fasta"):
                
				if seq_record.id not in lb6:
                    
					sequence = str(seq_record.seq)
                    
					if sequence[p] == current[0]:
                        
						flag = 1
        
			if flag == 0:
                
				line = 'cluster6' + '\t' + current[0] + '\t' + str(p+1) + '\n'
				output.writelines(line)

	#LB7         

	l = 1
	current = []
	a = 0
	while l < len(lb7):
        
		seq = lb7[l]
		current.append('')
		current[a] = seq[p]
		a+=1
		l+=2

	if current.count(current[0]) == (len(lb7)/2):
    
		flag = 0
		for file in glob.glob('*.fasta'):
            
			for seq_record in SeqIO.parse(file, "fasta"):
                
				if seq_record.id not in lb7:
                    
					sequence = str(seq_record.seq)
                    
					if sequence[p] == current[0]:
                        
						flag = 1
        
			if flag == 0:
                
				line = 'cluster7' + '\t' + current[0] + '\t' + str(p+1) + '\n'
				output.writelines(line)

            

	print('Checking SNPs at position ', p+1)
	p+=1


output.close()

