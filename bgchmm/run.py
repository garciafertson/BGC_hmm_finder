#!/usr/bin/env python
import subprocess
import os, sys
from bgchmm.train.hmm_train import hmm_multinomial_train
from bgchmm.train.log_ratio import logratio 
from bgchmm.find.hmm_find import hmm_multinomial_predict

class Run:
    def __init__(self, args):
        self.args=args
        self.dirdata= os.path.join(os.path.dirname(os.path.realpath(__file__)))+"/../data/"
    def find(self): 
        infile=self.args.input
        bgctype=self.args.type
        cutoff=self.args.cutoff
        output=sel.args.output
        dirdata=self.dirdata
        hmm_multinomial_predict(pfam_seqs, output, bgctype, cutoff, dirdata)

    def train(self):
        infile=self.args.infile
        output=self.args.output
        perc=self.args.perc_cutoff
        dirdata=self.dirdata
        logratio(infile,output,dirdata,perc)
        hmm_multinomial_train(infile,output,dirdata)

    def main(self):
        if self.args.subparser_name=="find":
            self.find()
            print("predicting %s contigs in %s\n" %(self.args.type, self.args.input))
        elif self.args.subparser_name=="train":
            self.train()
            print("calculating %s emmision matrix\n" %self.args.output)



