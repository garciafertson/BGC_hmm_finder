#!/bin/bash
#$ -N hmmscan_prot2fa
#$ -cwd
#$ -l h_vmem=8G
#$ -j Y
#$ -pe mpi 12

. ~/.profile

	/usr/bin/hmmscan --domtblout $1.hmmscan \
		--noali \
		--cut_tc \
		-o $1.hmmscan.log \
		--cpu 12 \
		$PFAM_A $1

