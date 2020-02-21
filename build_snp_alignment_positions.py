#!/usr/bin/env python3

#Zimpel et al., 2020. "Global distribution and evolution of Mycobacterium bovis lineages".
#Author: Aureliano Guedes

#instructions to generate the snps files are described in the "Snps_pipeline" file.
#input: vcf files
#reference: fasta format
#outputs: snp matrix and tsv file containing SNPs positions

import pandas as pd
from Bio import SeqIO
import sys
import os
import argparse
from argparse import RawTextHelpFormatter

def parser_cli():
    parser = argparse.ArgumentParser(description= 
                                     '''Build trimmed alignment based in reference genome in fasta format and list of SNPs tables.
                                     Also, generates a table mapping the nt for all SNPs positions for each genome.
                                     USAGE: 
                                    1- Using directory filled of table:
                                        $ build_snp_alignment -f reference_genome.fasta -o trimed_aln.fa -d PATH_DIR_TO_TABLES/ > SNP_Matrix.tsv
                                    2- Passing list of table manually:
                                        $ build_snp_alignment -f reference_genome.fasta -o trimed_aln.fa -s table1.tsv table2.tsv table3.tsv > SNP_Matrix.tsv
                                    3- Passing files in a directory as list of table:
                                        $ build_snp_alignment -f reference_genome.fasta -o trimed_aln.fa -s $( \ls -1rtd PATH_DIR_TO_TABLES/* | xargs ) > SNP_Matrix.tsv
                                    4- Passing directory and list:
                                        $ build_snp_alignment -f reference_genome.fasta -o trimed_aln.fa -d PATH_DIR_TO_TABLES/ -s table1.tsv table2.tsv > SNP_Matrix.tsv
                                    5- Passing files in a directory as list of table and add other files:
                                        $ build_snp_alignment -f reference_genome.fasta -o trimed_aln.fa -s $( \ls -1rtd PATH_DIR_TO_TABLES/* | xargs ) table1.tsv table2.tsv > SNP_Matrix.tsv
                                        Author: Aureliano Guedes''', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-s", "--SNP",
                        nargs='+',
                        help= '''Insert a list of tables''')
    parser.add_argument("-o", "--output",
                        required=True,
                        help= '''Insert output name.''')
    parser.add_argument("-d", "--directory", help="Insert directory path contain only SNPs table.")
    parser.add_argument("-f", "--fastafile",
                        required=True,
                        help= '''Insert reference genome fasta file.''')
    args = parser.parse_args()
    return args 

args = parser_cli()
genomeref_file = args.fastafile
output = args.output
table_list = args.SNP
table_list = list() if table_list is None else table_list
directory = args.directory


if (directory and not directory.endswith('/')):
    directory = directory + '/'

if (not directory and not table_list):
    print ("You should give me a list of file or directory contain thoses files")
    sys.exit(1)
else:
    if (directory):
        if ( os.listdir(directory) ):
            fromdir = [x for x in os.listdir(directory) if not x.startswith('.')]
        else:
            fromdir = list()
            directory = './'
    else:
        fromdir = list()
        directory = './'
    lista = fromdir + table_list
    if (not lista):
        print ("No file found.")


def build_aln(alphabet, df, name):
    '''This function changes the original nucleotide by the SNP.'''
    for i in df.index:
        alphabet[df.iloc[i, 0] - 1 ] = df.iloc[i, 1]
    return alphabet


#OLD
#######################################
#'''Get all tables with SNPs'''
# import list of tables
#lista = [x for x in os.listdir("SNP_matrix/") if not x.startswith('.')]
# set directory of of lists (to path)
#directory = './SNP_matrix/'
#######################################

#######################################
'''Load all tables'''
#Build pandas object with first table
bigt = pd.read_csv(directory+lista[0], sep='\t', header=None)
#concatenate the others tables
for i in lista[1:]:
    tmp = pd.read_csv(directory+i, sep='\t', header=None)
    bigt = pd.concat([bigt, tmp])
#Get non-redundant list of SNPs positions
keep = [x -1 for x in sorted(list(set(bigt[2].tolist())))]
#######################################

#######################################
'''Load reference genome'''
#Load reference genome fasta
seq = SeqIO.to_dict(SeqIO.parse(genomeref_file, "fasta"))
#Extract sequence
genome_name = list(seq.keys())[0]
genome_alphabet = list(seq[genome_name].seq)
#Check length of genome
genome_size = len(genome_alphabet)
#######################################

#######################################
'''Process all data and build trimmed alignment
     keep olny positions which was map some mutation 
     in any genome.'''
#Create filehandle to output
#ftout=open()
sys.stdout.write("genome\t"+"\t".join([str(i) for i in keep])+"\n")
fout=open(output, "w")
# get each genome per iteration 
for i in lista:
    # open table and use in function which returns full genome with SNPs
    path = directory + i
    df = pd.read_csv(path, sep='\t', header=None, index_col=False)
    tmp = genome_alphabet[:]
    g = build_aln(tmp, df[[2,4]], i)
    #Trim the genome keep only positions where have SNP at last one genome of our list
    aln = list()
    for x in keep:
        aln.append(g[x])
    #Save result in fasta format
    sys.stdout.write(i+"\t"+"\t".join(aln)+"\n")
    fout.write(">"+i + "\n" + "".join(aln)+"\n")
#    fout.write(build_aln(tmp, df[[2,4]], i))
#ADD reference genome
aln = list()
for i in keep:
    aln.append(genome_alphabet[i])
fout.write(">"+genome_name + "\n" + "".join(aln)+"\n")
fout.close()
#######################################

