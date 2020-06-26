##############################################
# Compute neoantigen quality
#
# Directory structure:
#
# data:
#   neoantigen-data file and iebd-epitope file.
#
# alignments:
#    precomputed blastp alignments for all neoantigens, split into files for each sample.
#    blastp -query <alignments/neoantigens_samplename.fasta> -db data/iedb.fasta -outfmt 5 -evalue 100000000  -gapopen 11 -gapextend 1 > <alignments/neoantigens_samplename.xml>
#
# src: 
#   source code folder
#
# output:
#   source code output folder
##############################################

# fitness model paramaters
a=26.
k=1.

python src/main.py data/SupplementaryTable1.txt alignments $a $k > output/neontigenQuality.txt


