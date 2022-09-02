#!/bin/bash -l

# qsub options
#$ -l h_rt=24:00:00
#$ -l mem_total=252G
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup --------------------------------------------------------
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
BDIR="pipeline/bowtie"
SDIR="pipeline/samtools"
LDIR="pipeline/lofreq"

# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index name
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

## print job info for output log --------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs -------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

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

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## make index ---------------------------------------------------
mesg "STEP 1: CREATE INDEX"
mesg "Creating Bowtie2 index: $BDIR/$BOWTIE"

# load bowtie2
module load bowtie2

# create directory for bowtie2 if it doesn't already exist
if [ ! -d "$BDIR" ]
then
  mkdir -p "$BDIR"
fi

# build command and execute
CMD="bowtie2-build --quiet '$FASTA' '$BDIR/$BOWTIE'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Bowtie2 index"
echo ""

## install samtools ---------------------------------------------
# create directory for samtools if it doesn't already exist
if [ ! -d "$SDIR" ]
then
  mkdir -p "$SDIR"
fi

# quick check for samtools; if installed skip this step
if [ ! -f "$SDIR/bin/samtools" ]
then
  # only print message once we know we need to install samtools from source
  mesg "STEP 2: INSTALL SAMTOOLS"
  mesg "Installing SAMtools v1.15.1 from source"

  # get the absolute path for installation
  PFIX="$(pwd)/$SDIR"

  # pull tarball, expand, and enter directory
  wget --quiet "https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2"
  tar -xf "samtools-1.15.1.tar.bz2"
  rm "samtools-1.15.1.tar.bz2"
  cd "samtools-1.15.1"

  # configure and install
  ./configure --quiet --prefix "$PFIX"
  make --quiet
  make install --quiet
  checkcmd "SAMtools installation"

  # clean up by removing samtools directory
  cd ..
  rm -r "samtools-1.15.1"
  cd "$CDIR"
fi

## all done -----------------------------------------------------
mesg "Setup complete!"
module list
