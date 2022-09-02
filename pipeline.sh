#!/bin/bash -l

# qsub options
#$ -l h_rt=24:00:00
#$ -pe omp 4
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
SAMTOOLS="pipeline/samtools/bin/samtools"
LOFREQ="pipeline/lofreq/lofreq"

# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -p PRIMERS -i INDEX -f FASTA -o ODIR -s SAMPLE -x R1 [-y R2]
Please submit the job from the pipeline directory!

arguments:
  -p primer BED file
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -h print this message and exit
"

# parsing arguments
while getopts ":hp:i:f:o:s:x:y:" opt 
do 
  case ${opt} in 
    p ) BED="${OPTARG}"
      ;;
    i ) IDX="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    s ) ID="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
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

# double-check that samtools exists
if [ ! -f "$SAMTOOLS" ]
then
  err "SAMtools not installed. Please re-run setup.sh"
fi

# double-check that lofreq exists
if [ ! -f "$LOFREQ" ]
then
  err "LoFreq missing from pipeline directory."
fi

# BED file
if [ -z "$BED" ]
then
  err "No primer BED file provided"
elif [ -f "$BED" ]
then
  mesg "Primer BED file: $BED"
else
  err "Invalid primer BED file: $BED"
fi

# bowtie2 index should have 6 files
if [ -z "$IDX" ]
then
  err "No bowtie2 index provided"
elif [ "$(ls -1 ${IDX}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Bowtie2 index: $IDX"
else
  err "Invalid bowtie2 index; run setup.sh first: $IDX"
fi

# reference sequence
if [ -z "$REFSEQ" ]
then
  err "No reference FASTA provided"
elif [ -f "$REFSEQ" ]
then
  mesg "Reference FASTA: $REFSEQ"
else
  err "Invalid reference FASTA: $REFSEQ"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided"
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# sample ID
if [ -z "$ID" ]
then
  err "No sample ID provided"
else
  mesg "Using sample ID: $ID"
fi

# R1/R0 FASTQ file
if [ -z "$R1" ]
then
  err "No FASTQ file provided (-x)"
elif [ -f "$R1" ]
then
  mesg "First FASTQ file: $R1"
else
  err "Invalid first FASTQ file: $R1"
fi

# check for R2 FASTQ file
if [ -z "$R2" ]
then 
  mesg "Only 1 FASTQ file detected. Running as unpaired."
elif [ -f "$R2" ]
then
  mesg "Second FASTQ file: $R2"
else
  err "Invalid second FASTQ file: $R2"
fi

# done checking inputs!
mesg "Done checking inputs!"
VAR="$ODIR/$ID"
echo ""

## alignment ----------------------------------------------------
mesg "STEP 1: ALIGN TO GENOME"

# load bowtie2 
module load bowtie2

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 4 -x '$IDX' -U '$R1' > '$VAR.sam'"
else # paired
  CMD="bowtie2 --threads 4 -x '$IDX' -1 '$R1' -2 '$R2' > '$VAR.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Alignment"

# compress SAM to BAM
CMD="$SAMTOOLS view --threads 4 -b -h '$VAR.sam' > '$VAR-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$VAR.sam"
echo ""

## soft-clip primers --------------------------------------------
mesg "STEP 2: CLIP PRIMERS"

# soft-clipping primers with samtools ampliconclip
CMD="$SAMTOOLS ampliconclip --threads 4 -b '$BED' '$VAR-raw.bam' -o '$VAR-clipped.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Primer clipping"
rm "$VAR-raw.bam"
echo ""

## process BAM --------------------------------------------------
mesg "STEP 3: PROCESS BAM"

# sort BAM
CMD="$SAMTOOLS sort --threads 4 '$VAR-clipped.bam' > '$VAR-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$VAR-clipped.bam"

# score indels to get final BAM
CMD="$LOFREQ indelqual --dindel --ref '$REFSEQ' '$VAR-sorted.bam' > '$VAR.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$VAR-sorted.bam"

# index final BAM
CMD="$SAMTOOLS index '$VAR.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage -------------------------------------------
mesg "STEP 4: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="$SAMTOOLS depth --threads 4 -a -H '$VAR.bam' > '$VAR.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus -------------------------------------------
mesg "STEP 5: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="$SAMTOOLS consensus --use-qual --min-depth 10 --call-fract 0.5 --output '$VAR-tmp.fa' '$VAR.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"

# update consensus header
mesg "Updating consensus header"
echo ">$ID" > "$VAR.fa"
cat "$VAR-tmp.fa" | grep "^[^>]" >> "$VAR.fa"
rm "$VAR-tmp.fa"
echo ""

## quantify SNVs ------------------------------------------------
mesg "STEP 6: QUANTIFY SNVs"

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="$LOFREQ call-parallel --pp-threads 4 --call-indels --min-cov 10 --ref '$REFSEQ' '$VAR.bam' > '$VAR.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"
echo ""

## package version ----------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
echo "LoFreq"
$LOFREQ version
echo ""
$SAMTOOLS --version
echo "" 
