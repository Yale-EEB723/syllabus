## *Syllabus*

# Comparative Genomics

	Yale University, Spring 2019
	E&EB 723

	Time: Thursday 1-3
	Location: ESC 100

	Instructor: Casey Dunn
	Prepend the subject line of all course related emails with "genomics: "
	Office hours: Thursdays 10AM-11:30AM


The field of evolutionary biology is increasingly drawing on genomic data. The field of genomic biology is becoming more evolutionary as genomes are sequenced for a broader diversity of organisms. This course focuses on the evolution of genome sequence and function at macroevolutionary timescales, with an emphasis on building practical computational skills for genomic and phylogenetic comparative analyses. There will be more focus on using phylogenies to understand genome evolution than on using genomes to build phylogenies.

## Technical details

### Git

This course is organized with [github education tools](https://education.github.com/guide). The course will make heavy use of [git](https://git-scm.com/) for sharing, communication, and collaboration. All students need to have a github account, preferably one that is registered to their Yale email address so that they get the [full academic features](https://education.github.com/discount).

### HPC

The class will be given student cluster accounts on one of the Yale High Performance Compute (HPC) clusters. Access to both compute and storage resources will last the course of the semester, after which data need to be copied elsewhere if they are important enough to be saved.

## Course format

Classes will consist of lectures, student led discussions, and computational labs.

## Course site

All materials for the course, including the syllabus, are available at the [course site](https://github.com/Yale-EEB723). The syllabus will be updated as the course progresses, please check it weekly. Please submit suggestions and corrections for the class via the [issue tracker](https://github.com/Yale-EEB723/syllabus/issues).

## Assignments

All assignments will be distributed and submitted via GitHub. Here are the
basic steps:

1. Go to the class repository for the assignment, for example
https://github.com/Yale-EEB723/finalproject . Click the "Fork" button
to create a fork. This creates your own copy of the repository with its own url,
for example https://github.com/YourGitHubUsername/finalproject .

2. Clone the repository to your laptop. On
your repository website, click on the green "Clone or download" button and copy
the link there. It will be something like `https://github.com/YourGitHubUsername/finalproject.git`

3. Open a Terminal window on your laptop. `cd` to a directory where you keep your
git repositories. Then clone your form of the repository (substitute the link
copied above): `git clone https://github.com/YourGitHubUsername/finalproject.git`

4. Edit the files and sage changes. Use `git add` to add any new changes.

5. Commit the changes with `git commit -am "my message"` where "my message" describes
the changes you made. commit often as you work.

6. Push the changes back to GitHub with `git push`.

7. Once the assignment is complete, return to the repo page at GitHub, for example
https://github.com/YourGitHubUsername/finalproject. Make sure all your changes are
reflected there. Save, commit, and push again if not.

8. Click the "New pull request" button. Then click the "Create pull request" button
and submit the pull request. This will notify the instructor that you have
submitted the assignment. Some assignments will be submitted multiple
times as sequential tasks are completed, which will require multiple pull
requests.


### Final Projects

Each student will work on a project, either in collaboration or individually. The project must relate to one or more themes covered in the course of the class. Final project plans will be presented in week 3 of the course. After the team and topic are set, fork the repository at https://github.com/Yale-EEB723/finalproject to create a repository for your project. Submit the project as a pull request. Your forked repository can be private if it includes unpublished original data, but if private all course members should be granted access so they can view it and provide feedback.

The final project can consider a research project already in progress (eg something that is part of thesis research), analysis of publicly available data, analysis of simulated data, development and testing of statistical methods or software, etc. Ideally each project will advance the existing research goals of each student, or advance an interesting topic identified in the course.

Here are some suggested final project ideas:

- A deep dive on a specific technical challenge of de novo genome sequencing and assembly,
eg repeats or heterozygosity

- Assembly and annotation of an original or publicly available de novo genome

- Examine the evolution of genome structure (eg synteny, size, intron distribution, etc)
with phylogenetic methods

- Explore the fit of models of evolution to genomic or functional genomic data

- Test phylogenetic hypotheses with genomic data

- Analyze one or more categories of functional genomic data in a phylogenetic context
to test hypotheses about the evolution of genome function

- Use comparative functional genomic and/or genomic data to identify genes that
may relate to specific phenotypes

- Compare within population genome variation to variation at broader phylogenetic scales

### Exercises

Several exercises will be assigned. They will usually be started in class, and then
due by the following week.

## Reading

Reading includes manuscripts, book chapters, online resources, and videos to be watched ahead of class. The readings will be posted by the Monday before each class. In most weeks, the 15-20 minute discussion of the reading will be led by a group of students. All students will get a chance to participate in these groups. A bibliogrpay at the end of this document includes a variety of references that readings can be selected from.


In addition to reading assigned for each class, the following will be used as references throughout the course:

- Haddock, SHD and CW Dunn (2011). Practical Computing for Biologists. [amazon](http://www.amazon.com/Practical-Computing-Biologists-Steven-Haddock/dp/0878933913/ref=sr_1_1) *I wrote this book with my colleague Steve Haddock as an introduction to general computing skills for biologists. If you are not already comfortable at the command line then you should get this book as a reference.*

- Whickham, H (2017). R for Data Science. http://r4ds.had.co.nz *This book is free online at the provided link. It is an excellent introduction to data analysis in R, and more broadly how to think about data structure and analysis. It presents a coherent introduction to the Tidyverse, a set of R packages for general data manipulation and analysis. Our R coding will follow conventions in this book.*


## Setting up your computer

In this course, you will perform some exercises and analyses on your own laptop
in class, and some on the cluster. Below are instructions on how to set up
your laptop.

Setup an account at [GitHub](https://github.com) using your educational email address.

Install [git](https://git-scm.com/downloads).

Install the [Atom](https://atom.io) text editor.

Install [Docker](https://www.docker.com/get-started).



## Schedule

### Week 1, January 17 - Intro to comparative genomics. Questions, methods, history

#### Reading

Discussion leader: Casey Dunn
- Felsenstein, J. 1985. Phylogenies and the Comparative Method. American Naturalist, 125:1–15. https://www.jstor.org/stable/2461605
- Dunn CW, Zapata F, Munro C, Siebert S, Hejnol A. 2018 Pairwise comparisons across species are problematic when analyzing functional genomic data. Proc. Natl. Acad. Sci., 115:E409–E417. https://www.doi.org/10.1073/pnas.1707515115


#### Agenda and notes

- Introductions
- Discussion of course goals and structure
  - Readings
  - Projects
  - Class formats
- Course logistics
  - Bring laptop to each class
  - github account required
  - YCRC account setup description, needed prior to week 3
- Review readings
- Overview of computational framework and tools

#### Exercises

First, confirm that docker is working by running a container:

    docker run -it rocker/rstudio /bin/bash

Next, we will walk through regular expressions in the exercises at
https://github.com/Yale-EEB723/syllabus/blob/master/regular_expressions.txt .


### Week 2, January 24 - Sequencing technology and applications

#### Reading

Discussion leader: Ian Gilman
- Goodwin et al. 2016. Coming of age: ten years of next-generation sequencing
technologies. Nature Reviews Genetics.  https://doi.org/10.1038/nrg.2016.49
*This review covers a lot of ground. Focus on the bits about Illumina, PacBio,
and Oxford Nanopore Technologies (ONT).*

- Practical Computing for Biologists Chapters 2-3. *This optional reading
provides background for, and builds on, the regular expressions exercises.*

#### Agenda and notes

- Sequencing technology and instruments
  - Conceptual overview
    - Single molecule vs. populations of molecules
    - Multiplexing
    - Sequencing overview
      - Sample preparation
      - Data acquisition
      - Data preprocessing
      - Base calling
      - Read processing (trimming, binning, etc) and export
      - Downstream analysis (application specific)
    - Tradeoffs
      - Cost (initial and realtime)
      - Read length
      - Error rate and error profile (base miscalls, phasing noise, homopolymer length, etc)
      - Throughput
      - Hands-on limitations (sample prep cost, instrument portability, ease of use, run time, etc)
  - Current sequencing technologies
    - Illumina
      - https://www.youtube.com/watch?v=fCd6B5HRaZ8
      - The recent shift to [reduced  colors](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology/2-channel-sbs.html)
    - PacBio
      - https://www.youtube.com/watch?v=NHCJ8PtYCFc
      - Very long molecules can be sequenced or the same molecule can be
      sequenced repeatedly with  [Circular Consensus Sequencing](http://files.pacb.com/software/smrtanalysis/2.2.0/doc/smrtportal/help/!SSL!/Webhelp/Portal_PacBio_Glossary.htm)
    - Oxford Nanopore
      - https://www.youtube.com/watch?v=GUb1TZvMWsw
      - https://www.youtube.com/watch?v=hs0FdiTHMbc

- Genome sequencing
  - Challenging factors
    - Large size
    - Repeats
    - Heterozygosity
    - Tissue limitation

- Take homes
  - Focus on inputs and outputs, not intermediates. For example, assembly quality
  is usually much more important than read quality.
  - Take a wholistic perspective on costs, including your time. Saving a bit of
  money on sequencing can sometimes incur large data analysis costs, for example.
  - Focus your time and resources on what differentiates your project from others.
  - Always, always be thinking about the end goal and evaluate intermediate
  decisions in terms of these final objectives.

#### Exercises

We will walk through regular expressions in the exercises at
https://github.com/Yale-EEB723/syllabus/blob/master/regular_expressions.txt .

To get the files for the exercises, make a local clone of the syllabus repository:

    git clone https://github.com/Yale-EEB723/syllabus.git


#### Other

- Hand out git chapter
- Ask students to prepare to present preliminary ideas on final projects next
week


### Week 3, January 31 - Practical computing skills

#### Reading

- Jain et al. Nanopore sequencing and assembly of a human genome with ultra-long
reads. Nature Biotech. https://doi.org/10.1038/nbt.4060. (Discussion leaders:
Edgar Benavides and Vincent Dimassa)

-  Practical Computing for Biologists Chapters git chapter draft
(to be provided as hard copy in previous week)

-  Practical Computing for Biologists Chapters chapters 4-6, 20  *This optional
reading provides background on working in bash and remote access to computers.*


#### Agenda and notes

- Quick (<1 minute) description of final project plan

- Working on your laptop
- Working on the cluster via your account
- Getting started with git

[Ben's cluster presentation](https://docs.google.com/presentation/d/1uA7LKhTbDO8YHD4u3b0DtOQRG90dpOlsauij9Bo_jvs/edit?usp=sharing)

### Week 4 - Genome assembly

#### Reading

- Sedlazeck et al. 2018. Piercing the dark matter: bioinformatics of long-range sequencing and mapping. Nature Reviews Genetics. https://doi.org/10.1038/s41576-018-0003-4 (Discussion leader: Ava Ghezlayagh )


#### Agenda

The agenda:
- Review paper
- Discuss final projects
- Walk through alignment exercise
- Discuss items below about assembly

#### Notes

We often want  to know the full sequence of a genome, but data are fragmented and redundant because:
- DNA isolation leads to fragmentation by mechanical and chemical processes
- Sample preparation leads to fragmentation by mechanical and chemical processes, often deliberately to adjust the length of input molecules
- Sequencing often doesn’t span the full length of a molecule, due to technical limitations, damaged template, damaged sequencing
- The same regions of the genome are sequenced multiple times because:
  - It is far easier to randomly sequence regions of the genome than to systematically tile sequencing effort across the genomes. To sequence with enough depth to ensure coverage everywhere, some places will have quite deep coverage
  - Sequencing is error prone. Sequencing each spot multiple times enables error detection and correction.

In the broad sense,
- Sequencing takes large biological molecules in and generates character strings in computer memory that are redundant overlapping estimates (reads) of subsequence of of the original molecules. Underlying this is a generative model, ie an idea of how molecular structure impacts observed read sequences.
- Assembly is in some respects the reverse process - it starts with reads and generates an estimate of the sequences of the
input biological molecules. It runs the generative model in reverse, and generates descriptions of molecules in computer memory
rather than


Common tasks:
- de novo assembly: going from reads (and sometimes additional new structural
  data) to a genome assembly without reference to an existing assembly.
- mapping: Tiling reads onto an existing reference genome sequence. Used to assess
how the reads cover the reference sequence (more on this later when we discuss functional genomics),
or to identify how the genome from which the reads are derived differs from the reference genome.
- reference based assembly: assembling a genome by reference to an existing genome, usually by
mapping new reads to an existing reference genome sequence. Does not require as much sequencing
depth as de novo assembly, also much easier computationally.


The challenge:
- Find similar sequences
- Categorize differences between similar sequences according to whether they are
  - Sample prep errors
  - Sequencing errors
  - Different regions of the genome (eg paralogy, repeats)
  - Different alleles
  - Mixtures (eg somatic variations)
- Estimate the sequence of the original molecules


The general assembly process (see Figure 1 from https://doi.org/10.1038/s41576-018-0003-4)
- Identify overlaps between reads
- Construct a string graph
  - Nodes are unambiguous sequences
  - Edges are possible connections between those sequences
  - Each path through the string graph is a possible assembly
- Contig construction
  - modify and traverse the string graph to derive contig sequences
    - Errors create bubbles with low coverage that can be popped
    - Some nodes cannot be combined with neighboring nodes because there isn't
    enough information to know which alternative path to take, and those nodes
    are emitted as contigs
    - Adjacent nodes (ie nodes connected directly by edges) can often be combined
  - Ends of contigs usually due to ambiguity of some sort, there are multiple
  paths and the assembler doesn't know which to take so it chops them all off
    - Sources of ambiguity include error, repeats not spanned by reads, or heterozygosity
    - As read quality and length improve, extent of contigs determined in larger part by heterozygosity
      - Makes heterozygous-aware assembly even more important
- Scaffolding
  - Physical ordering of contigs, sometimes introducing gaps
  - Usually based on additional structural information
    - HiC
    - Optical mapping

Error correcting can occur at multiple steps
- Improvements to base calling
- Use short high quality reads to correct lower quality long reads
  - Risks introducing errors, eg by mistaking one instance of a repeat for another similar instance of the same repeat
- Use short high quality reads to correct contigs derived from lower quality long reads
  - Same risks as above
- Use lower quality long reads to correct each other (requires greater depth)
  - In the last year this is where things have started to head

Characterizing an assembly
- contiguity (eg N50)
- completeness (eg BUSCO)
- correctness (eg base level, structural, phasing)

Genome assembly challenging factors
  - Large size
  - Repeats
  - Heterozygosity
  - There are tools for assessing all of these before attempting a full assembly, eg https://github.com/schatzlab/genomescope


Phasing
  - Collapse haplotypes into a single consensus. Can introduce many errors and fragment the assembly
  - Assemble into regions that are collapsed and unzipped
    - Can arbitrarily resolve into primary assembly and alternate haplotigs, or pseudophased diploid genome
    - Phase into two haploid genomes


Core algorithmic concepts
- Similarity and extension
  - Identify similarity be identification of similar seed sequences in different reads
    - Expensive because each read needs to be compared to every other read
    - Exact matches are really fast, but often need to allow for variation due to errors
    - A few methods
      - Canu uses MHAP, a kmer method
      - Falcon etc use DALIGNER, dynamic programming and kmers
  - Extension searches for regions of similarity beyond the seed.
    - It essentially sees if if is possible to zip the reads together starting at the seed
    - Extension is generally not as expensive is initial identification of similarity, because extension is
- k-mers
  - Short sequences of length k (often 15-70 nucleotides)
  - Very cheap to work with
    - Defined memory footprint
    - If short relative to frequency of errors, can focus on exact matches
    - Easy to code
  - Hash tables are sorted lists of k-mer sequences, often with a count of how many times the sequence exists
  - de Bruijn assembly
  - identification of similarity seeds


### Week 5 - Genome annotations

#### Reading

Yandell and Ence. 2012. A beginner's guide to eukaryotic genome annotation.
Nature Review Genetics. https://doi.org/10.1038/nrg3174 (Discussion leader: Diego and Jasmine)

#### Agenda

The agenda:
- Discuss assembly notes from last class
- Go over today's paper
- Identify two papers for next week
- Walk through forking, cloning, and editing final project git repository. Make
  preliminary edits and push them.
    - Add a few sentences for Goal and Data sections. Make a list of tools you
		  plan to use in the Methods section. Send a pull request when you have added these stubs.
- Explain how to finish alignment exercise from last class as an assignment
- Discuss genome annotation notes

#### Notes

A genome assembly is just a big fasta file. Much more interesting with annotations.

Annotation goals include:
- Identifying repeats
  - Biologically interesting
  - Technically important for understanding genome assembly and so that
  repeats can be masked for some downstream analyses that they can negatively
  impact
- Identifying protein coding genes
  - Build inventory of protein coding genes
  - Identify introns, exons, promoters
  - Predict mRNA structure
- Identifying other regions,
  - Noncoding RNA sequences
  - Promoters
  - etc...
- Understanding genome structure
  - Centromeres, telomeres, how scaffolds map to chromosomes

The set of annotations can be encoded in standard formats and loaded into graphical browsers
or interrogated computationally (as when comparing genome features across species)

Annotations can be made
- With *ab initio* methods based on the understanding of what these structures look like or
  particular properties they have
- Based on new evidence, like mRNA-seq
- Based on comparison to a reference of similar sequences
  - Looking for known repetitive DNA elements
  - Blasting known protein coding genes

Assessing annotations
  - Sensitivity = (True positives) / ( True positives + False negatives )
	  - BUSCO is a common tool for this
	- Specificity = (True positives) / ( True positives + False positives )
	- Accuracy = (Sensitivity + Specificity) / 2. If accuracy is high you know you have good sensitivity and specificity, but if it is low could be problems with sensitivity or specificity or both.


Repeats
- The extent and composition of repeats varies widely across species
- Many repeats are well-known conserved sequences
  - LINEs Long Interspersed Nuclear Elements. Retrotransposons about 7kb in length. They are
    transcribed and encode a reverse transcriptase that facilitates integration at new sites.
  - SINEs Short Interspersed Nuclear Elements. 100-700 bp retrotransposons. They do not have their
    on reverse transcriptase, so they are dependent on those of other elements
- Repeat annotation is complicated by the fact that different instances of the same repeat
  are not identical, and may be quite different
- Tools for identifying and masking repeats
  - RepeatMasker

Protein coding genes
- Multiple nested annotation steps
  - Gene finding
  - What is entire region needed for gene to function, including upstream upstream regions,
    transcribed region, and downstream regions
  - What region is transcribed
    - Identify promoter where transcription is initiated
      - For eukaryotic protein coding genes, the binding site of RNA Polymerase II
    - Identify the location where transcription is terminated and mRNA is polyadenylated
  - Introns and exons
  - Mature mRNA sequence prediction
  - Protein sequence prediction
- Commonly used tools
  - AUGUSTUS http://bioinf.uni-greifswald.de/augustus/ . Developed for human genome.
	- MAKER http://www.yandell-lab.org/software/maker.html . Annotation workflow that integrates multiple tools.
	- Funanotate https://funannotate.readthedocs.io/en/latest/ . Initially made for fungal genomes.



### Week 6 - Genome annotations

#### Reading

Discussion leaders: Elyse Parker and Dan MacGuigan

Holt and Yandell. 2011. MAKER2: an annotation pipeline and genome-database
management tool for second-generation genome projects. BMC Bioinformatics.
https://doi.org/10.1186/1471-2105-12-491

Kim et al. 2018. The genome of common long-arm octopus Octopus minor.
GigaScience. https://doi.org/10.1093/gigascience/giy119

MarkDown guide - https://guides.github.com/features/mastering-markdown/

#### Agenda

#### Notes

Success of *ab initio* methods
- Tests of *ab initio* methods work well in species with well annotated genomes, because the gene training datasets are so good. It is easy to find
genes when you know what they look like.
- *Ab initio* methods work less well alone in poorly studied species without good existing gene models to server as training datasets. Additional evidence, like RNAseq, greatly improves outcome in these projects.

Annotation Edit Distance (AED)
- A measure the distance between intron and exon coordinates of two annotations
for the same gene, where 0 is identical and 1 is completely different
- This distance indicate incongruence between annotation methods and is
interpreted as uncertainty about an annotation
- Genes with higher AED tend to
	- Have fewer PFAM domains
	- Change more in subsequent annotations that include additional information
	- Have less evidence of orthologs in closely related species


MAKER is now part of the GMOD family of genomic tools - http://gmod.org/wiki/Main_Page

Ian notes that `conda install -c bioconda maker` works well for installation

A detailed guide on how to use MAKER - https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
- Identify and mask repeats
- Build gene models based on RNAseq data and homology to genes of related species
- Train *ab initio* programs based on these evidence-based models
- Rerun MAKER with *ab initio* prediction
- Rerun MAKER again to refine everything


Mapping
- Mapping is a general task that comes up in many types of genomic analyses.
- It is a highly asymmetric sequence comparison, usually between reads and
reference sequences

- Burrows-Wheeler transform
  - Transforms the sequence to a sorted string that is easy to compress
  - Can work with compressed strings, which is more computationally efficient
  - Has the very special property that the original string can be recovered from
  the sorted form


### Week 7 - Functional genomics

Discussion leader: Andrew Verdegaal

RNAseq
- Randomly sequence RNA from a cell or collection of cells
- In nearly all cases mRNA is reverse transcribed into complimentary DNA (cDNA)
and that is analyzed with with short read DNA sequencing (eg Illumina)
- The frequency of reads for each gene are used as a proxy for abundance of transcripts of each gene
(often referred to as expression)
- If a suitable high quality reference is available, it can also be used as a proxy for abundance
of splice variants.
- Replaced Expressed Tags (ESTs), shotgun sequencing of cDNA with Sanger sequencing
- Full length sequencing of cDNA (eg PacBio isoseq, Oxford Nanopore cDNA sequencing) provides
a better understanding of slice variants than short reads, but lower throughput provides less
statistical power.
- What it does and doesn't do
	- It does not measure absolute transcript count
	- It does not measure the relative expression of different genes
	- It does come close to measuring the differential expression of the same gene
	in different tissues.
- Analysis steps
	- Map to reference, which can be gene models from genome or transcriptome assemblies.
	Many short read aligners can do this, eg Bowtie
	- Process mapping data to derive the counts of reads for transcripts and genes.

	ATACseq
	- https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142727.mb2129s109


### Week 8 - Phylogenetic methods

### Week 9 - Phylogenetic comparative methods

### Week 10 - Genome evolution

### Week 11 - Comparative functional genomics

### Week 12 - Project

### Week 13 - Project


## Bibliography

You can suggest references to add to this list via a pull request or the [issue tracker](https://github.com/Yale-EEB723/syllabus/issues). The intent of this bibliography is to serve as a resource for class participants in their own work and as a list of potential readings for class.

### Genome sequencing

Goodwin et al. 2016. Coming of age: ten years of next-generation sequencing
technologies. Nature Reviews Genetics.  https://doi.org/10.1038/nrg.2016.49

Heather et al. 2015. The sequence of sequencers: The history of sequencing DNA.
Molecular Cell. https://doi.org/10.1016/j.molcel.2015.05.004

Reuter et al. 2015. High-Throughput Sequencing Technologies. Molecular Cell.
https://doi.org/10.1016/j.molcel.2015.05.004

Shendure et al. 2017. DNA sequencing at 40: past, present and future. Nature.
https://doi.org/10.1038/nature24286

Vurture et al. 2017. GenomeScope: fast reference-free genome profiling from
short reads. Bioinformatics. https://doi.org/10.1093/bioinformatics/btx153

### Genome assembly

Alkan et al. 2011. Genome structural variation discovery and genotyping. Nature Reviews Genetics. https://www.nature.com/articles/nrg2958

Bradnam et al. 2013. Assemblathon 2: evaluating de novo methods of genome
assembly in three vertebrate species. GigaScience.
https://doi.org/10.1186/2047-217X-2-10

Gurevich et al. 2013. QUAST: quality assessment tool for genome assemblies
Bioinformatics. https://doi.org/10.1093/bioinformatics/btt086

Koren et al. 2018 De novo assembly of haplotype-resolved genomes with trio binning.
Nature Biotechnology. https://doi.org/10.1038/nbt.4277 *Good overview of phasing*

Paajanen et al. 2019. A critical comparison of technologies for a plant genome
sequencing project. GigaScience. https://doi.org/10.1093/gigascience/giy163

Schatz et al. 2010. Assembly of large genomes using second-generation sequencing.
Genome Research. https://doi.org/10.1101/gr.101360.109

Sedlazeck et al. 2018. Piercing the dark matter: bioinformatics of long-range
sequencing and mapping. Nature Reviews Genetics. https://doi.org/10.1038/s41576-018-0003-4

Sohn and Nam 2016. The present and future of de novo whole-genome assembly.
Briefings in Bioinformatics. https://doi.org/10.1093/bib/bbw096

### Example genome projects

Edgar et al. 2018. Single-molecule sequencing and optical mapping yields an
improved genome of woodland strawberry (Fragaria vesca) with chromosome-scale
contiguity. Gigascience. https://doi.org/10.1093/gigascience/gix124

Jain et al. Nanopore sequencing and assembly of a human genome with ultra-long
reads. Nature Biotech. https://doi.org/10.1038/nbt.4060 . (For an interesting
set of followup analyses see https://genomeinformatics.github.io/na12878update/ )

Kim et al. 2018. The genome of common long-arm octopus Octopus minor.
GigaScience. https://doi.org/10.1093/gigascience/giy119

Mohr et al. 2017. Improved de novo Genome Assembly: Linked-Read Sequencing
Combined with Optical Mapping Produce a High Quality Mammalian Genome at
Relatively Low Cost. BioRxiv. https://doi.org/10.1101/128348

Wenger et al. 2018. Highly-accurate long-read sequencing improves variant
detection and assembly of a human genome. https://www.biorxiv.org/content/10.1101/519025v2

Jiang et al. 2018. A Hybrid de novo Assembly of the Sea Pansy (Renilla muelleri)
Genome. https://doi.org/10.1101/424614

Vertebrate Genomes Project - https://vgp.github.io/genomeark/




### Genome annotation

Aken et al. 2016. The Ensembl gene annotation system. Database.
https://doi.org/10.1093/database/baw093

Holt and Yandell. 2011. MAKER2: an annotation pipeline and genome-database
management tool for second-generation genome projects. BMC Bioinformatics.
https://doi.org/10.1186/1471-2105-12-491

Mudge and Harrow. 2016 The state of play in higher eukaryote gene annotation.
https://doi.org/10.1038/nrg.2016.119

Yandell and Ence. 2012. A beginner's guide to eukaryotic genome annotation.
Nature Review Genetics. https://doi.org/10.1038/nrg3174





### Genome structure

Li et al. 2011. Chromosome Size in Diploid Eukaryotic Species Centers on the
Average Length with a Conserved Boundary. Molecular Biology and Evolution.
https://doi.org/10.1093/molbev/msr011


Harmston et al. 2017. Topologically associating domains are ancient features
that coincide with Metazoan clusters of extreme noncoding conservation.
Nature Communications. https://doi.org/10.1038/s41467-017-00524-5

Rowley and Victor. 2018. Organizational principles of 3D genome architecture.
Nature Reviews Genetics. https://doi.org/10.1038/s41576-018-0060-8


### Functional genomics


### Reconstructing the history of genome evolution

Kim et al. 2017. Reconstruction and evolutionary history of eutherian chromosomes. PNAS. https://doi.org/10.1073/pnas.1702012114

Demas et al. 2018. Reconstruction of avian ancestral karyotypes reveals differences in the evolutionary history of macro- and microchromosomes. BMC Genome Biology. https://doi.org/10.1186/s13059-018-1544-8

O'Connor et al. 2018. Reconstruction of the diapsid ancestral genome permits chromosome evolution tracing in avian and non-avian dinosaurs. https://www.nature.com/articles/s41467-018-04267-9

### Other sequencing applications and perspectives

Cheifet 2019. Where is genomics going next? Genome Biology.
https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-019-1626-2

Pomerantz et al. 2018. Real-time DNA barcoding in a rainforest using nanopore sequencing: opportunities for rapid biodiversity assessments and local capacity building GigaScience. https://academic.oup.com/gigascience/article/7/4/giy033/4958980


### Making associations between genomes and phenotypes with comparative methods

### Phylogenetic comparative methods

### General computational skills and methods

## Links and other general resources

MarkDown guide - https://guides.github.com/features/mastering-markdown/ . Useful for writing text that
will be pushed to GitHub, among many other things.
