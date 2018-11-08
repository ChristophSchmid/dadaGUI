#$ -S /bin/sh
#$-cwd 
#$-e log
#$-o log
# name for job
#$ -N phylotree.R
#
# pe request for slots, means processors, with
# minimum to maximum.
# e.g. 2    means two processors
#      2-   means 2 upto number of processor restricted by pe
#      2-8  means 2 upto 8 processors
#
# The number of slots / processors restriction can be retrieved by
# qmon under Parallel Environment Configuration choosing scampi and
# there by the number of Slots.
#
# ATTENTION: If the number of slots ist exceeded by single or
# maximum number, the job does not run.
#
#$ -pe hmp 20
#
# V1.00 written by Christoph Schmid
# ---------------------------

#check input arguments by amount (min. 3)
[ $# -lt 2 ] && { printf "qsub_phylotree_DADA2.sh V1.00\nUsage: qsub -cwd qsub_phylotree_DADA2.sh -i <input file> -o <output path>\n"; exit 1; }

#parse arguments using getopts
while getopts ':i:o:d:btp' OPTION ; do
  case "${OPTION}" in
    i) INPUT=$OPTARG;;
    o) OUTPUT=$OPTARG;;
		*) echo "Unknown command $OPTION\nUsage: qsub -cwd qsub_phylotree_DADA2.sh -i <input file> -o <output path>\n"
	esac
done

#execute taxonomy.R R-Script with corresponding command line arguments
Rscript /project/genomics/Christoph/DADA2/phylotree.R -i $INPUT -o $OUTPUT

#run raxml-ng software to compute phylogenetic tree from starting NJ tree and sequence alignment
#(LD_PRELOAD necessary to use locally installed libstdc++.so.6 version!)
LD_PRELOAD='/project/genomics/Christoph/DADA2/raxml_ng/libstdc++.so.6.0.25' /project/genomics/Christoph/DADA2/raxml_ng/raxml-ng-mpi --msa $OUTPUT'/sequenceAlignment.fasta' --tree $OUTPUT'/startingTreeNJ.tre' --model GTR+G+I