#$ -S /bin/sh
#$-cwd 
#$-e log
#$-o log
# name for job
#$ -N inference.R
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
# V1.1 version control added, m option removed
# ---------------------------

#check input arguments by amount (min. 3, max 6)
[ $# -lt 3 ] && { printf "qsub_inference_DADA2.sh V1.00\nUsage: qsub -cwd qsub_inference_DADA2.sh -i <input file> -o <output path> -p <number of plots> <v / s / c options>\n"; exit 1; }

#parse arguments using getopts
VERSION=1.6.0
while getopts 'i:p:o:v:sc' OPTION ; do
  case $OPTION in
    i) INPUT=$OPTARG;;
    p) PLOTS=$OPTARG;;
    o) OUTPUT=$OPTARG;;
    v) VERSION=$OPTARG;;
    #m) MERGEOFF="-m";;
		s) SEQOFF="-s";;
		c) CHIMOFF="-c";;
		*) echo "Unknown command $OPTION"
	esac
done

#execute taxonomy.R R-Script with corresponding commands
Rscript /project/genomics/Christoph/DADA2/inference.R -i $INPUT -p $PLOTS -o $OUTPUT -V $VERSION  ${SEQOFF+"$SEQOFF"} ${CHIMOFF+"$CHIMOFF"}