
#!/bin/bash

### Name of the job
#$ -N gene_ssa_sr2_cl1
### Declare job is non-rerunable
#$ -r n
### Export all environment variables to batch job
#$ -V
#$ -o /space1/home/tjkwon/ConservedRegulatoryProgram/results/gene_ssa_sr2_cl1.out
#$ -e /space1/home/tjkwon/ConservedRegulatoryProgram/results/gene_ssa_sr2_cl1.err
### E-mail notification on job abort
#$ -m a
#$ -M tjkwon@cmmt.ubc.ca

echo $HOSTNAME

/space/devel/oPOSSUM3/scripts/standalone/opossum_gene_ssa.pl -s human -g /space1/home/tjkwon/ConservedRegulatoryProgram/muscle_ref_validated_combined.ensembl -bnr 2000 -d /space1/home/tjkwon/ConservedRegulatoryProgram/results -o gene_ssa_sr2_cl1.results.txt -h gene_ssa_sr2_cl1.hits.txt -l gene_ssa_sr2_cl1.log

exit 0