#$ -S /bin/sh
#$-cwd 
#$-e log
#$-o log
# name for job
#$ -N taxonomy.R
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
#$ -pe hmp 12
#
# V1.00 written by Christoph Schmid
# V1.1 version control added
# ---------------------------

#check input arguments by amount (min. 3)
[ $# -lt 3 ] && { printf "qsub_taxonomy_DADA2.sh V1.00\nUsage: qsub -cwd qsub_taxonomy_DADA2.sh -i <input file> -o <output path> -d <databasename>\n"; exit 1; }

#parse arguments using getopts
VERSION=1.6.0
while getopts ':i:o:d:V:btp' OPTION ; do
  case "${OPTION}" in
    i) INPUT=$OPTARG;;
    d) DATABASE=$OPTARG;;
    o) OUTPUT=$OPTARG;;
    V) VERSION=$OPTARG;;
    b) BIOM="-b";;
    t) TREE="--tree";;
    p) PHYLOSEQ="--noPS";;
		*) echo "Unknown command $OPTION\nUsage: qsub -cwd qsub_taxonomy_DADA2.sh -i <input file> -o <output path> -d <databasename> -V <dada2 version> -b -t -p\n"
	esac
done

#execute taxonomy.R R-Script with corresponding commands
Rscript /project/genomics/Christoph/DADA2/taxonomy.R -i $INPUT -o $OUTPUT -d $DATABASE -V $VERSION ${BIOM+"$BIOM"} ${TREE+"$TREE"} ${PHYLOSEQ+"$PHYLOSEQ"}

#create biom file
if [ "$BIOM" ]
then
  export PATH=/opt/qiime/1.9.1/bin/:$PATH
  cd "$OUTPUT"
  biom convert -i seqTabClean_forBiom.txt -o seqTabClean.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
  
fi
