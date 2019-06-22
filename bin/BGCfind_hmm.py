#!/usr/bin/env python
################################################
#BGCfind_HMM
#subparser 1 TRAIN train a BGC type
#subparser 2 FIND find a BGC type
################################################
__author__= "Fernando Garcia Guevara"
__credits__= "Fernando Garcia Guevara"
__email__= "garciafertson near gmail.com"
__status__= "Development"

import argparse
import sys
import os
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
#print (sys.path)

#import module for running program
import bgchmm
from bgchmm.run import Run


def print_header():
    bgchmm.version()
def phelp():
    print("""
                    BGCfind
    This program uses a HMM model to predict a particular BGC type 
    in contigs from a binned genome.

    for more information type:

    BGCfind_hmm.py find -h

    or

    BGCfind_hmm.py train -h

    """)

##the next section add parsers and subparser to the program

if __name__ == '__main__':
    parser= argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='bgchmm %s' % bgchmm.__version__)
    subparser= parser.add_subparsers(help="sub-comand help:", dest='subparser_name')
    
    ##### PARSER 1 "find" subparser######
    find_parser = subparser.add_parser( 'find',
        description='Predict presence of specified BGC type in multifasta file, return headers',
        help='Predict presence of specifief BGC type in contigs from bin',
        epilog='''
    #####################################################################\n
    for running call:\n
    BGCfind_hmm.py find --input binned_contigs.fasta --type BGC_type \n
    #####################################################################\n
     ''')
    input_options = find_parser.add_argument_group('input options')
    input_options.add_argument('--input',"-i",
            dest='input',
            metavar='FILE.FASTA',
            help="",
            required=True)
    input_options.add_argument('--type','-t',
            dest='type',
            metavar='BGC_TYPE',
            help="Specify BGC type, from the following list:\
            PUFA acyl_amino_acids amglyccycl aminocoumarin arylpolyene\
            bacteriocin blactam bottromycin butyrolactone cyanobactin\
            ectoine furan fused glycocin head_to_tail hserlactone \
            indole ladderane lantipeptide lassopeptide linaridin melanin\
            microcin microviridin nrps nucleoside oligosaccharide other\
            otherks pbde phenazine phosphoglycolipid phosphonate ppysks\
            proteusin resorcinol sactipeptide siderophore t1pks t2pks\
            t3pks terpene thiopeptide transatpks",
            required=True)

    output_options = find_parser.add_argument_group('output options')
    output_options.add_argument('--cutoff','-c',
            dest='cutoff',
            metavar='FLOAT',
            help="specify cutoff value for HMM model",
            default=0.05,
            type=float)
    output_options.add_argument('--output', '-o',
            dest='output',
            metavar='PREFIX',
            help='specify the prefix for the name of output files')


    ######### PARSER 2 construct hmm emission matrices
    train_parser= subparser.add_parser('train',
            description="recieve a folder with csv files of Pfam sequences \
            from a BGC type and outputs a numpy emission matrix for HMM model",
            help="specify a folder containign csv with Bfam sequences of a BGC type",
            epilog='''
    #################################################################\n
    for running call:\n
    BGCfind_hmm.py train --folder path/to/bgc_type/files/pfam.csv
    ''')
    input_train = train_parser.add_argument_group('input train')
    input_train.add_argument('--infile', '-i',
            dest="infile",
            metavar='FOLDER',
            help="Path to folder containig the PFAM sequences from BGC type \
            in csv format",
            type=str,
            required=True)
    logratio_train=train_parser.add_argument_group('logratio')
    logratio_train.add_argument('--logratio', '-L',
            dest='logratio',
            help="run log ratio analysis of Pfam",
            action="store_true")
    output_train=train_parser.add_argument_group('output')
    output_train.add_argument('--output', '-o',
            dest="output",
            metavar="PREFIX",
            help='specify the prefix for output file names',
            required=True)

    #check whether --help is needed
    if (len(sys.argv)==1 or sys.argv[1]== '-h' or sys.argv[1]== '--help'):
        phelp()
    #call Run module passing the arguments in here
    else:
        args=parser.parse_args()
        Run(args).main()

        
