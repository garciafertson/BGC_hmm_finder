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

def plot_logratio(values, values_bg, outname):
    kwargs=dict(density=True, stacked=True, alpha=0.5)
    #plot histogram with bgc type observed values
    plt.hist(values, bins=50, **kwargs, color='g')
    #plot density distribution with random logratio
    #plt.hist(values_bg,bins=50, **kwargs) #bins=max(values)
    density = stats.gaussian_kde(values_bg)
    xmin,xmax=plt.xlim()
    n, x, _ = plt.hist(values_bg, bins=np.linspace(xmin, xmax, 50),
                   histtype=u'step', density=True)
    plt.plot(x, density(x))
    #set cutoff top 90 percentile of random logratio
    v=np.array(values_bg)
    per=np.percentile(v,90)
    plt.vlines(per,ymin=0,ymax=1,color='r',
            lw=0.3, label="90th percentile")
    #fit normal density function to random logratio
    mu,sd= stats.norm.fit(values_bg)
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu, sd)
    plt.plot(x, p, 'k', linewidth=1.5)
    #plot axis and titles
    plt.axis([min(values)-1, max(values)+1, 0, 0.5])
    plt.title("BGC type: "+ outname )
    plt.xlabel("log-ratio observed/background")
    plt.ylabel("density")
    figure_name=outname+"logratio.png"
    plt.savefig(figure_name,dpi=300, figzise=(6.42,4))
    plt.close()

def random_csv(nfrags,npfams,dirdata):
    #the function create a tmp.csv file with a sequence of pfam
    #from random genomes with similar characteristics as 
    #the set of observed bgc: number of observation "nfrags" and
    # number of pfams per observation "npfams")
    #in folder "data/pfam_csv_genome/" exist the csv genome files
    genomes=[]
    genomes_csv=dirdata+"names_pfam_csv_genomes.txt"
    with open(genomes_csv,"r") as glist:
        genomes = [line.rstrip() for line in glist]
    fragments=[]
    for i in range(nfrags):
        genome=np.random.choice(genomes) #choice 1 random genome from list
        name=dirdata+"pfam_csv_genomes/"+genome #open file 
        genomesize = len(open(name).readlines()) #lines minus header
        np_len=np.random.normal(npfams, 3, 1) #(mean, sd, ntrials)
        start=np.random.randint(0, genomesize-int(np_len[0]))
        length=int(np_len[0]) #pfam length
        with open(name, "r") as csv:
            next(csv, None) #skip header
            index=-1
            for record in csv:
                index+=1
                if index < start:
                    continue
                elif index >= start+length:
                    break
                else: #return sequences
                    fragments.append(record)
    with open("tmp.csv", "w") as rand_csv:
        rand_csv.write("tmp pfam sequence fragments in csv format\n")
        for record in fragments:
            rand_csv.write("%s" %record)

def histogram(infile):
    with open(infile) as domtbl:
        reader=csv.reader(domtbl, delimiter=',')
        next(reader, None)
        counter=Counter(line[6] for line in reader if line[6])
    return(counter)

def logratio(infile,output,dirdata,perc):
    #needs csv file
    #load all genomes Pfam Background relative frequency
    name=dirdata+"back_freq_dict2000genomes.pkl" #dict back relative freq of pfams
    with open(name,'rb') as back: background=pickle.load(back)
    #####Create .csv random sample observation, (aproximately) same size as csv
    size = len(open(infile).readlines())
    with open(infile) as domtbl:
        reader=csv.reader(domtbl, delimiter=",")
        next(reader, None)
        bgc_obs=Counter(line[0] for line in reader if line[0])
    nfragments=len(bgc_obs)
    npfams=size/len(bgc_obs)
    ##**implement repetitions on random observation, currently only 1 
    ##repeat for i in range(100) and  list values_ran_logratio.extend(tmplist)
    nreps=10
    values_random_logratio=[]
    for i in range(nreps): #sample nreps random sets
        random_csv(nfragments, npfams, dirdata) #sample random set
        RANobs=histogram("tmp.csv")             #counter
        total=float(sum(RANobs.values()))       #total number of pfams
        random={}
        for pfam in RANobs:
            random[pfam]=RANobs[pfam]/total #relative freq of pfams in random set
        ran_logratio={}
        for pfam in random: #logratio random observes agains background expected
            if pfam in background:
                ratio=random[pfam]/background[pfam]
                ran_logratio[pfam]=(math.log(ratio, 2))
        values_random_logratio.extend(list(ran_logratio.values()))
        print(i)
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
            obs_logratio[pfam]=(math.log(ratio, 2))
    #set log ratio 80% quantile cutoff
    values_obs_logratio=list(obs_logratio.values())
    v=np.array(values_random_logratio)
    per=np.percentile(v,perc) #80th pecentile: representative PFAMs for BGCtype
    rep_PFAM= [pf for pf in obs_logratio.keys() if obs_logratio[pf] > per ]
    name=dirdata+"pfam_csv_bytype/"+output+"PF_log.pkl"
    #save list with BGCtype representative PFams
    with open(name,"wb") as pflog: pickle.dump(rep_PFAM, pflog) #pkl with repPFAMs 
    #plot values for logartio analysis
    plot_logratio(values_obs_logratio, values_random_logratio, output)

