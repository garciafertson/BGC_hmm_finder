'''the script takes in 1 the count list of pfams from 
background 2000 complete genomes selected at random
2 the count list obs pfams from random fragments from genomes
3 the count list obs pfams from bgc fragments

calculate relative frequency of pfams from all lists
calculate log of ratio between rel-freqs obs/bg only list of observed

estimate distribution of random log(obs/bg), observed and adjust normal,
calculate cutoff for pvalue <0.00001 1e-5

Retrieve list of pfams from BGC with log(obs/bg) values above cutoff 
from random log(obs/bg)'''

import argparse
import numpy as np
from collections import defaultdict
import csv
import math
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bglist", dest="background", required=True,
                      help="pfam count list from background genomes", metavar="FILE")
parser.add_argument("-i", "--bgc_obs", dest="bgc_obs", required=True,
                      help="pfam count list from bgc type", metavar="FILE")
parser.add_argument("-r", "--rand_obs", dest="rand_obs", required=True,
        help="pfam count list from random genome fragments", metavar="FILE")
parser.add_argument("-o", "--output", dest="output" , required=False,
        default="logratio_output", help="output root name ", metavar="STRING")
options = parser.parse_args()

background={}
bgc_obs={}
rand_obs={}

total_background=0
total_bgc_obs=0
total_rand_obs=0

with open(options.background, "r") as f:
    next(f)
    for line in f:
        line=line.rstrip()
        pfam,freq=line.split("\t")
        freq=int(freq)
        background[pfam]=freq
        total_background+=freq
    
with open(options.bgc_obs, "r") as f:
    next(f)
    for line in f:
        line=line.strip()
        pfam,freq=line.split("\t")
        freq=int(freq)
        bgc_obs[pfam]=freq
        total_bgc_obs+=freq
with open(options.rand_obs,"r") as f:
    next(f)
    for line in f:
        line=line.strip()
        pfam,freq=line.split("\t")
        freq=int(freq)
        rand_obs[pfam]=freq
        total_rand_obs+=freq

rand_logratio={}
for pfam in rand_obs.keys():
    if pfam in background.keys():
        b=background[pfam]/total_background
        r=rand_obs[pfam]/total_rand_obs
        rand_logratio[pfam]=(math.log(r/b, 2))    
bgc_logratio={}
for pfam in bgc_obs.keys():
    if pfam in background.keys():
        b=background[pfam]/total_background
        r=bgc_obs[pfam]/total_bgc_obs
        bgc_logratio[pfam]=(math.log(r/b, 2))


#pois=np.random.poisson(lam=1, size=5000)
#plot distriburion of rand_logratio
    #plt.hist(pois,bins=max(pois), **kwargs, color='g' ) 

values_bg=list(rand_logratio.values())
values=list(bgc_logratio.values())
kwargs=dict(density=True, stacked=True, alpha=0.5)

plt.hist(values_bg,bins=50, **kwargs) #bins=max(values)
plt.hist(values, bins=50, **kwargs, color='g')

v=np.array(values_bg)
per=np.percentile(v,99)
print(per)
mu,sd= stats.norm.fit(values_bg)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.norm.pdf(x, mu, sd)
plt.plot(x, p, 'k', linewidth=2)

plt.axis([min(values)-1, max(values)+1, 0, 0.5])
plt.title("Pfam histogram")
plt.xlabel("log-ratio obs/bckgnd ")
plt.ylabel("density")
plt.vlines(per,ymin=0,ymax=1,color='r',
        lw=0.3, label="95th percentile")
#plt.show()
figure_name=options.output+"logr_hist.png"
plt.savefig(figure_name)
plt.close()

#save tsv pfam with log ratio above cutoff
labels_perc= [label for label in bgc_logratio.keys() if bgc_logratio[label]>per ]
tablename=options.output+".logr_p01.csv"
with open(tablename, 'w') as out:
        for label in labels_perc:
            out.write("%s\t%s\t%s\t%s\n" %(label, bgc_logratio[label],bgc_obs[label],background[label] ))






