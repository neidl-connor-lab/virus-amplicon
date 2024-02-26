#!/bin/bash -l

# qsub options
#$ -l h_rt=6:00:00
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup -----------------------------------------------------------------------
# functions
mesg () { echo -e "[MSG] $@"; }
err () { echo -e "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}

# pre-set variables
LOFREQ="pipeline/lofreq"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index path and prefix
  -h show this message and exit"

# parsing arguments
while getopts ":hf:b:" opt
do
  case ${opt} in
    f ) FASTA="${OPTARG}"
      ;;
    b ) BOWTIE="${OPTARG}"
      ;;
    h ) echo "${HELP}" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## print job info for output log -----------------------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
# load bowtie2
module load bowtie2/2.4.2
checkcmd "Loading bowtie2" 

# load samtools module 
module load htslib/1.18
module load samtools/1.18
checkcmd "Loading samtools"

# check reference FASTA
if [ -z "$FASTA" ]
then
  err "No reference FASTA provided"
elif [ -f "$FASTA" ]
then
  mesg "Valid reference FASTA file: $FASTA"
else
  err "Invalid reference FASTA file: $FASTA"
fi

# check that bowtie prefix was provided
if [ -z "$BOWTIE" ]
then
  err "No bowtie2 index name provided"
else
  mesg "Bowtie2 index prefix: $BOWTIE"
fi

# make bowtie index output directory if necessary
if [ ! -d "$(dirname $BOWTIE)" ]
then
  mesg "Creating bowtie2 index directory: $(dirname $BOWTIE)"
  mkdir -p "$(dirname $BOWTIE)"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## unpack LoFreq if necessary --------------------------------------------------
# check for lofreq tarball
if [ -f "$LOFREQ.tar.bz2" ]
then 
  mesg "LoFreq tarball detected. Expanding."
  tar -xf "$LOFREQ.tar.bz2" -C "$(dirname $LOFREQ)"
  checkcmd "LoFreq decompression"
  # remove tarball
  rm "$LOFREQ.tar.bz2"
fi

# check that the lofreq executable works
if [ ! -z "$($LOFREQ/lofreq version 2> /dev/null)" ]
then
  mesg "LoFreq is ready to go!"
else
  err "Problem with LoFreq: $LOFREQ"
fi
echo ""

## make index ------------------------------------------------------------------
mesg "Creating Bowtie2 index: $BOWTIE"

# build command and execute
CMD="bowtie2-build --quiet '$FASTA' '$BOWTIE'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Bowtie2 index"
echo ""

## all done --------------------------------------------------------------------
mesg "Setup complete!"
module list
