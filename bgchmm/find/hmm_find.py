#This modules runs predict hidden state 0->nonBGC and 1->BGC
from hmmlearn import hmm
import numpy as np
import sys
import pickle

def hmm_multinomial_predict(infile, output, bgctype, cutoff,dirdata):
    #predict hidden state  per pfam-contig sequnce
    EMname=dirdata+bgctype+"EM.pkl"
    with open(EMname,"rb") as EM: emission_matrix=pickle.load(EM)
    TMname=bgctype+"TM.pkl"
    with open(TMname,"rb") as TM: transition_matrix=pickle.load(TM)
    SPname=output+"SP.pkl"
    with open(SPname, "rb") as SP: start_prob=pickle.load(SP)
    IDictname=output+"IndexDict.pkl"
    with open(IDictname, "rb") as IDict: index_dict=pickle.load(IDict)
    model=hmm.MultinomialHMM(n_components=2)
    np.random.seed(42)
    model = hmm.MultinomialHMM(n_components=3)
    model.startprob_ = start_prob
    model.transmat_ = transition_matrix
    model.emissionprob_= emission_matrix
    #for each pfam-seq-contig predict BGC state and return contig seq id
 
