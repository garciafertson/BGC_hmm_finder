"""steps followed to build emission probability matrix from PFam.csv sequences
   recieve a folder with csv files containing pfam sequences
   estimate  emission matrix
   saves variable 
"""

from collections import defaultdict
from collections import Counter
import pickle # save emission matrix in python varable format
#import defaultdict # fill in pfam to local_pfam_id_to emmissionmatrix
import numpy as np
import pickle
import csv
import math

def encode_emission(pfamcounter, backdict, repbgc_list):
    #the size of matrix same as len(pfamdict)
    #create emission matrix and pfam index dict
    i=0
    index_dict={}
    bgc_freq=[]
    back_freq=[]
    total=float(sum(pfamcounter.values()))
    pfamdict={}
    pfamset=set(pfamcounter.keys())
    for pfam in pfamset:
        pfamdict[pfam]=pfamcounter[pfam]/total
    #**construct emission matrix around logratio overrepresented PFAMs (BGC-characteristic)
    #**if not in overrepresented set then in "other"(non BGC-characteristic) category
    for pfam in repbgc_list:
        if pfam in pfamcounter and pfam in backdict:
            index_dict[pfam]=i
            bgc_freq.append(pfamdict[pfam])
            back_freq.append(backdict[pfam])
            i+=1
    index_dict["other"]=i #len(back_freq)-1
    bgc_freq.append(1-sum(bgc_freq))
    back_freq.append(1-sum(back_freq))
    emission_matrix=np.array([back_freq,bgc_freq])
    return(index_dict, emission_matrix)

def hmm_multinomial_train(infile,output,dirdata):
    name=dirdata+"back_freq_dict2000genomes.pkl"#**probably need back counter 2000 genomes
    print(name)
    with open(name,"rb") as pfam_bg_freq: bg_freqdict=pickle.load(pfam_bg_freq)
    name=dirdata+"pfam_csv_bytype/"+output+"PF_log.pkl"
    with open(name,"rb") as rep_bgc: repbgc_list=pickle.load(rep_bgc)
    #calculate pfam count for input file
    with open(infile,"r") as domtbl:
        reader=csv.reader(domtbl, delimiter=',')
        next(reader, None) #skip header
        pfamcounter=Counter(line[6] for line in reader if line[6])
    #HMM model 2 states 0->nonBGC 1->BGC
    #emission matrix from frequencies of Pfam in BGC predicted by Antismash
    #transition matrix and start probability similar to ClusterFinder
    #observations len(index_dict) pfams in analysis
    index_dict,emission_matrix=encode_emission(pfamcounter,bg_freqdict, repbgc_list)
    transition_matrix= np.array([[0.9999, 0.00001],
                                 [0.05  , 0.95   ]])
    start_prob=np.array([0.7,0.3])
    #save in output folder files with .pkl extension
    #1 EM Emission numpy matrix #2 TM Transition numpy matrix
    #3 SP Start probability numpy vector #4 IndexDict index to Pfam dictionary
    EMname=dirdata+output+"EM.pkl"
    with open(EMname,"wb") as EM: pickle.dump(emission_matrix, EM)
    TMname=dirdata+output+"TM.pkl"
    with open(TMname,"wb") as TM: pickle.dump(transition_matrix,TM)
    SPname=dirdata+output+"SP.pkl"
    with open(SPname,"wb") as SP: pickle.dump(start_prob,SP)
    IDictname=dirdata+output+"IndexDict.pkl"
    with open(IDictname, "wb") as IDict: pickle.dump(index_dict, IDict)

    #perform logratio analysis if specified
#    if logratio:
#        overrepr_pfam=log_ratio()
        #sample random genome fragments, hmmscan, hist and freq of pfams
        #plot log ratio of Pfam distribution in BGC type against random sample
        #return PFAM list of overrepresented PFAMS
        #compare against detection rules of BGCtype and filter

        ##for detection grep list-of-overrepr-pfam in binned-contigs.fasta and return fasta id
        ##overlap in prediction between HMM and logratio


