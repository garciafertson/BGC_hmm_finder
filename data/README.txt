The folder "pfam_csv_bytype" contain the pfam sequences in csv format, from Antismash DB 2 
with the following format:
contig_id,locus_tag,protein_id,gene_start,gene_end,gene_strand,pfam_id,domain_start,domain_end,evalue,bitscore

The remaining folders with bgc-type-names contain:
1Emission probability matrix for that type
2Transition probability matrix(fixed for all types)
and overestimated frequency of Pfam BGC positive
3start probability vector
4Python dictionary for PFAM emmission matrix, (allows transformation
from PFAM bgc list into numpy vector of observations)

