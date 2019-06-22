#recieve a text file with a list of genome files in fasta aa sequences 
#select a random genome from the list
#calculate total number of aa sequences in genome grep -c "^>"
#specify average BGC length in gene number, std 2 genes
#select a random start positon from the total of aa sequences in the genome
#decide a random length of genes normal distribution centeres in average BGC lenghth
#return a genome fragment file sequences of contiguos genes in fasta aa format
#repeat process for specified N number of genomes, 2 files per genomes
#return N*2 files with faa sequences and concatenate in single file

#with the resulting file do pfam histogram analisis
#do pfam histogram analisis with total pfam in 2000 random genome

import argparse
import subprocess
from Bio import SeqIO
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_genomes", dest="gnm_list", required=True,
            help="text file with list of full path of genome files in fasta format with aa sequences",
            metavar="FILE")
parser.add_argument("-n", "--number", dest="n_fragments", type=int, required=False, default=None,
            help="number of random fragments to extract", metavar="INT")
parser.add_argument("-o", "--output", dest="output", required=False, default="output",
                      help="Output fasta file name", metavar="FILE")
parser.add_argument("-s","--size", dest="size", type=int, default=25,
            help="expected fragment size in number of fasta sequences (default 25)", metavar="INT")
parser.add_argument("--std", default=3, 
            help="std deviation average numger of genes (default 3)", metavar="INT")
options = parser.parse_args()

genomes=[]
with open(options.gnm_list,"r") as glist:
    for line in glist:
        line=line.rstrip()
        genomes.append(line)
if not options.n_fragments:
    n_fragments=len(genomes)
else:
    n_fragments=options.n_fragments

fragments=[]

for i in range(n_fragments):
    genome=np.random.choice(genomes)
    cmd=["grep","-c",'^>', genome]
    proc=subprocess.Popen(cmd, stdout=subprocess.PIPE)
    genomesize,_=proc.communicate()
    genomesize=int(genomesize)
    np_len=np.random.normal(options.size, options.std, 1)
    start=np.random.randint(0, genomesize-int(np_len[0]))
    length=int(np_len[0])
    with open(genome, "r") as faa:
        index=-1
        for record in SeqIO.parse(faa, "fasta"):
            index+=1
            if index < start:
                continue
            elif index >= start+length:
                break
            else: #return sequences
                fragments.append(record)
outname=options.output+".fas"
SeqIO.write(fragments, outname, "fasta")

