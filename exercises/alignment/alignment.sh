# create a project directory, fill in below
export PROJ_DIR=/project/eeb723-seqaln

# I chose Thermus thermophilus, because it's cool and has a nice sized genome

## Genome reference
# https://www.ncbi.nlm.nih.gov/genome/genomes/461

# chose one, 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/604/845/GCA_900604845.1_TTHNAR1


###############
## Genome setup
###############
# Get genome files
mkdir -p $PROJ_DIR/genome && cd $PROJ_DIR/genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/604/845/GCA_900604845.1_TTHNAR1/*

# unpack genome DNA sequence & index
gunzip ./GCA_900604845.1_TTHNAR1_genomic.fna.gz 
ln -s GCA_900604845.1_TTHNAR1_genomic.fna Thermus_thermophilus_TTHNAR1.fa

# samtools index genome
samtools faidx Thermus_thermophilus_TTHNAR1.fa

# picard dictionary for genome
picard CreateSequenceDictionary \
    REFERENCE=Thermus_thermophilus_TTHNAR1.fa \
    OUTPUT=Thermus_thermophilus_TTHNAR1.dict

# build bowtie2 genome index
bowtie2-build Thermus_thermophilus_TTHNAR1.fa Thermus_thermophilus_TTHNAR1


###############
## Get sequence reads
###############
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5324768

# from "Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life."
# Nat Microbiol. 2017 Nov;2(11):1533-1542. doi: 10.1038/s41564-017-0012-7. Epub 2017 Sep 11.

# get reads from SRA
cd $PROJ_DIR
export SRR=SRR5324768
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ${SRR}


############
## Alignment
############
cd $PROJ_DIR
mkdir -p alignment

# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# http://www.htslib.org/doc/samtools.html
bowtie2 -x genome/Thermus_thermophilus_TTHNAR1 \
        -1 fastq/${SRR}_pass_1.fastq.gz \
        -2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
        --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
    | samtools view -hb - | samtools sort -l 5 -o alignment/${SRR}.bam
samtools index alignment/${SRR}.bam

# text view alignment with 
# samtools tview alignment/SRR5324768.bam genome/Thermus_thermophilus_TTHNAR1.fa
# type ? for help, q to quit


##########################
## Variant Calls with GATK
##########################
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
# Warning! This step will take seveeral minutes (5-20) depending on your hardware

cd $PROJ_DIR
mkdir -p variants 

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   --reference genome/Thermus_thermophilus_TTHNAR1.fa \
   --sample-ploidy 1 \
   --input alignment/${SRR}.bam \
   --output variants/${SRR}.vcf
