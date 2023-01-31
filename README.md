# virus-amplicon
v1.0 by Jackie T.

## Requirements

Computing cluster with a [module](https://www.bu.edu/tech/support/research/software-and-programming/software-and-applications/modules/) for [Bowtie2](https://doi.org/10.1038%2Fnmeth.1923) and [SAMtools](http://www.htslib.org/download/).

| Package   | Version |
| :-------- | :------ |
| bowtie2   | 2.4.2   |
| samtools  | 1.15.1  |

## Quick start

### 1. Download repo 

Cloning this repo will create a directory holding the pipeline script and assorted files. Move to the directory where you'd like to save the pipeline; here, we're using `workflows/`.

```
cd workflows/
git clone https://github.com/neidl-connor-lab/virus-amplicon.git
cd virus-amplicon
```

### 2. Run `setup.sh` to set up LoFreq and make an index

Create an alignment index from a reference sequence FASTA file, and write it to `pipeline/indices/`. Most viruses have an official NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/). This script will also set up LoFreq in your `pipeline/` directory. 

Run `setup.sh` with the `-h` flag to view the full list of options.

```
./setup.sh -h
```


| Flag | Argument                                  |
| :--- | :---------------------------------------- |
| `-P` | SCC project                               |
| `-N` | job name                                  |
| `-f` | genome reference FASTA file               |
| `-b` | bowtie2 index ID                          |

```
usage: qsub -P PROJECT -N JOBNAME setup.sh -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index name
  -h show this message and exit
```

Here is an example where we download the [SARS-CoV-2 RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) and then run `setup.sh`. If this is your first time running `setup.sh` in the directory, it will install SAMtools in the directory as well!

```
# download and decompress the reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz
gunzip GCA_009858895.3_ASM985889v3_genomic.fna.gz

# move the reference into your pipeline directory for safekeeping
mv GCA_009858895.3_ASM985889v3_genomic.fna pipeline/

# submit the indexing job
qsub -P test-project \
     -N test-index \
     setup.sh \
     -f pipeline/GCA_009858895.3_ASM985889v3_genomic.fna \
     -b sarscov2

# wait until the job is done
qstat -u $(whoami)

# check out the pipeline directory!
ls pipeline/*
```

If you would like to make another index, just run `setup.sh` again!

### 3. Run pipeline

> You _must_ run this script from the directory where the `pipeline/` directory was created in step 2.

View pipeline options and required arguments by running `pipeline.sh` with the `-h` flag.

```
./pipeline.sh -h
```

The help message indicates the required arguments and how to pass them:

| Flag | Argument                                  |
| :--- | :---------------------------------------- |
| `-P` | SCC project                               |
| `-N` | job name                                  |
| `-b` | path to primer file (`pipeline/bedfiles`) |
| `-i` | path to index created in step 2           |
| `-f` | reference FASTA file                      |
| `-o` | output directory                          |
| `-s` | sample ID to use as output file prefix    |
| `-x` | path to R1 or unpaired FASTQ file         |
| `-y` | path to R2 FASTQ file (paired-read only)  |

```
usage: qsub -P PROJECT -N JOBNAME ./pipeline.sh -p PRIMERS -i INDEX -f FASTA -o ODIR -s SAMPLE -x R1 [-y R2]
Please submit the job from the pipeline directory!

arguments:
  -b primer BED file
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -h print this message and exit
```

Here is an example where we're using the index we created in step 2. Paired FASTQ input files are in `input-files/`, and output files should go to `output-files/`. We're sThe example job is named `test-job`, and the project allocation used is `test-project`. The job output will be written to a file named `log-test-job.qlog`.

```
qsub -P test-project \
     -N test-job \
     pipeline.sh \ 
     -b primers.bed \
     -i pipeline/indices/sarscov2 \
     -f pipeline/GCA_009858895.3_ASM985889v3_genomic.fna \
     -o output-files/ \
     -s test-sample \
     -x input-files/r1.fq.gz \
     -y input-files/r2.fq.gz
```

## Pipeline steps

### 1. Align to reference

Raw FASTQ files are aligned to the previously-constructed reference using [Bowtie2](https://doi.org/10.1038%2Fnmeth.1923). All output files have the sample ID as a prefix. 

| Flag        | Meaning                |
| :---------- | :--------------------- |
| `--threads` | parallelize this job   |
| `-x`        | index path and prefix  |
| `-1`        | (paired) R1 FASTQ file |
| `-2`        | (paired) R2 FASTQ file |
| `-U`        | (unpaired) FASTQ file  |
| `*.sam`     | uncompressed alignment |

```
# paired
bowtie2 --threads 4 \
        -x pipeline/bowtie/index \
        -1 paired-r1.fq.gz \
        -2 paired-r2.fq.gz > \
        odir/id.sam

# unpaired
bowtie2 --threads 4 \
        -x pipeline/bowtie/index \
        -U unpaired.fq.gz > \
        odir/id.sam
```

The uncompressed SAM output is then compressed to BAM format.

```
# compress
samtools view --threads 4 \
              -b \
              -h \
              odir/id.sam > \
              odir/id-raw.bam
```

### 2. Soft-clip primers

Aligned reads in the `*-raw.bam` file are clipped using `samtools ampliconclip` to remove primers found in the BED file provided with the `-b` option.

| Flag        | Meaning              |
| :---------- | :------------------- |
| `--threads` | parallelize this job |
| `-b`        | primer BED file      |
| `-o`        | output BAM file      |
| `*-raw.bam` | input alignment      |

```
samtools ampliconclip --threads 4 \
                      -b primers.bed \
                      -o odir/id-clipped.bam \
                      odir/id-raw.bam
```

### 3. Process alignment

These steps are get the alignment file ready for coverage, SNV, and consensus calling. We use SAMtools to sort the BAM, LoFreq to score insertions and deletions, and then SAMtools again to index the alignment.

| Flag            | Meaning                      |
| :-------------- | :--------------------------- |
| `--threads`     | parallelize this job         |
| `*-clipped.bam` | alignment from step 2        |
| `*-sorted.bam`  | sorted alignment             |
| `--dindel`      | algorithm for scoring indels |
| `--ref`         | reference FASTA file         |
| `id.bam`        | final alignment file         |

```
samtools sort --threads 4 \
              odir/id-clipped.bam > \
              odir/id-sorted.bam

lofreq indelqual --dindel \
                 --ref reference.fa \
                 odir/id-sorted.bam > \
                 odir/id.bam

samtools index odir/id.bam
```

### 4. Calculate coverage

The `samtools depth` command calculates the aligned read depth for each nucleotide of the genome used for alignment. The output is an easy-to-analyze TSV table.

| Flag        | Meaning                 |
| :---------- | :---------------------- |
| `--threads` | parallelize this job    |
| `-a`        | include all nucleotides |
| `-H`        | include a file header   |
| `id.bam`    | final aligment file     |
| `id.tsv`    | coverage table          |

```
samtools depth --threads 4 \
               -a \
               -H \
               odir/id.bam > \
               odir/id.tsv
```

### 5. Assemble consensus

The `samtools consensus` command assembles a consensus by examining the reads aligned to each nucleotide in the reference sequence and calling the most frequent allele.

| Flag           | Meaning                      |
| :------------- | :--------------------------- |
| `--threads`    | parallelize this job         |
| `--use-qual`   | use quality scores           |
| `--min-depth`  | minimum aligned read depth   |
| `--call-fract` | minimum nucleotide frequency |
| `--output`     | output FASTA file            |
| `id.bam`       | final aligment file          |

```
samtools consensus --threads 4 \
                   --use-qual \
                   --min-depth 10 \
                   --call-fract 0.5 \
                   --output odir/id.fa\
                   odir/id.bam
```

### 6. Quantify SNVs

Use LoFreq to make a detailed table of consensus and sub-consensus SNVs.

| Flag            | Meaning                    |
| :-------------- | :------------------------- |
| `--pp-threads`  | parallelize this job       |
| `--call-indels` | include indels in output   |
| `--min-cov`     | minimum aligned read depth |
| `--ref`         | reference FASTA file       |
| `id.bam`        | final alignment file       |
| `id.vcf`        | output SNV table           |

```
lofreq call-parallel --pp-threads 4 \
                     --call-indels \
                     --min-cov 10 \
                     --ref reference.fa \
                     odir/id.bam > \
                     odir/id.vcf
```
