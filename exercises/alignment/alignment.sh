# create a project directory, fill in below
export PROJ_DIR=/insert/here
mkdir -p $PROJ_DIR

# I Chose Thermus thermophilus, because it's cool and nice sized genome

## Genome reference
# https://www.ncbi.nlm.nih.gov/genome/genomes/461

# chose one, 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/604/845/GCA_900604845.1_TTHNAR1

# Get genome files
mkdir $PROJ_DIR/genome && cd $PROJ_DIR/genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/604/845/GCA_900604845.1_TTHNAR1/*

# unpack genome DNA sequence
gunzip ./GCA_900604845.1_TTHNAR1_genomic.fna.gz 
ln -s GCA_900604845.1_TTHNAR1_genomic.fna Thermus_thermophilus_TTHNAR1.fa

# build genome index
bowtie2-build Thermus_thermophilus_TTHNAR1.fa Thermus_thermophilus_TTHNAR1


## Sample Reads
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5324768

# from "Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life."
# Nat Microbiol. 2017 Nov;2(11):1533-1542. doi: 10.1038/s41564-017-0012-7. Epub 2017 Sep 11.

# get reads from SRA
cd $PROJ_DIR
export SRR=SRR5324768
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ${SRR}

## Alignment
cd $PROJ_DIR
mkdir alignment

# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# http://www.htslib.org/doc/samtools.html
bowtie2 --threads $(expr $SLURM_CPUS_ON_NODE - 1) -x genome/Thermus_thermophilus_TTHNAR1 \
                                                  -1 fastq/${SRR}_pass_1.fastq.gz \
                                                  -2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
                                                  --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
    | samtools view -hb - | samtools sort -l 5 -o alignment/${SRR}.bam
samtools index alignment/${SRR}.bam

# text view alignment with 
# samtools tview alignment/SRR5324768.bam genome/Thermus_thermophilus_TTHNAR1.fa
# type ? for help, q to quit

## Variant Calls with GATK