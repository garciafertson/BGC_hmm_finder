#set log ratio analysis
#this needs to sample genome fragments from Genome files,
#alternatively take a pre processed set, prepare preprocessed set variable samples
#select a random genome pfam.csv, set random start uniform from sequence
#select random size normal mean and sd per sample
#combine into one csv
#log ratio analysis, ratio of obsbgc/background against rand/background

import pickle
import numpy as np
from collections import defaultdict
from collections import Counter
import csv
import math

import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

def random_csv(genomes_csv,nfrags,npfams,dirdata):
    genomes=[]
    with open(genomes_csv,"r") as glist:
        next(glist,None)
        genomes = [line.rstrip() for line in glist]
    fragments=[]
    for i in range(nfrags):
        genome=np.random.choice(genomes)
        name=dirdata+"pfam_csv_genomes/"+genome
        genomesize = len(open(name).readlines())
        np_len=np.random.normal(npfams, 3, 1)
        start=np.random.randint(0, genomesize-int(np_len[0]))
        length=int(np_len[0])
        name=dirdata+"pfam_csv_genomes/"+genome
        with open(name, "r") as csv:
            index=-1
            for record in csv:
                index+=1
                if index < start:
                    continue
                elif index >= start+length:
                    break
                else: #return sequences
                    fragments.append(record)
    with open("tmp.csv", "w") as ran_csv:
        ran_csv.write("tmp pfam sequence fragments in csv format\n")
        for record in fragments:
            ran_csv.write("%s" %record)


def histogram(infile):
    with open(infile) as domtbl:
        reader=csv.reader(domtbl, delimiter=',')
        next(reader, None)
        counter=Counter(line[6] for line in reader if line[6])
    return(counter)


def logratio(infile,output,dirdata):
    #needs csv file
    #load all genomes Pfam Background frequency
    #####Create .csv random sample observation, (aproximately) same size as csv
    size = len(open(infile).readlines())
    with open(infile) as domtbl:
        reader=csv.reader(domtbl, delimiter=",")
        next(reader, None)
        bgc_obs=Counter(line[0] for line in reader if line[0])
    nfragmets=len(bgc_obs)
    npfams=size/len(bgc_obs)
    genomes_csv=dirdata+"names_pfam_csv_genomes.txt"
    nreps=1
    with open(name,'rb') as back: background=pickle.load(back)
    ##**implement repetitions on random observation, currently only 1 
    ##repeat for i in range(100) and  list values_ran_logratio.extend(tmplist)
    random_csv(genomes_csv, nfragments, npfams, dirdata)
    RANobs=histogram("tmp.csv") #counter
    total=float(sum(RANobs.values()))
    random={}
    for pfam in RANobs:
        random[pfam]=RANobs[pfam]/total
    name=dirdata+"back_freq_dict2000genomes.pkl" #dict backfreq
    ran_logratio={}
    for pfam in random:
        if pfam in background:
            ratio=random[pfam]/background[pfam]
            ran_logratio[pfam]=(math.log(ratio, 2))
    values_ran_logratio=list(ran_logratio.values())
    
    ##logratio in obseved pfam bgctype 
    BGCobs=histogram(infile) #counter
    total=float(sum(BGCobs.values()))
    observations={}        
    for pfam in BGCobs:
        observations[pfam]=BGCobs[pfam]/total
    obs_logratio={}
    for pfam in observations:
        if pfam in background:
            ratio=observations[pfam]/background[pfam]
            obs_logratio[pfam]=(math.log(r/b, 2))
    values_obs_logratio=list(obs_logratio.values())
    v=np.array(values_obs)
    per=np.percentile(v,80) #80th pecentile: representative PFAMs for BGCtype
    rep_PFAM= [pf for pf in obs_logratio.keys() if obs_logratio[pf] > per ]
    name=dirdata+output+"PF_log.pkl"
    
    with open(name,"wb") as pflog: pickle.dump(rep_PFAM, pflog) #pkl with repPFAMs
    #plot_hist(values_obs_logratio, values_ran_logratio,output)
    




