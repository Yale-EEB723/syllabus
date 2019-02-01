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

## Projects

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

### Week 5 - Genome annotations

### Week 6 - Genome annotations

### Week 7 - Phylogenetic methods

### Week 8 - Phylogenetic comparative methods

### Week 9 - Genome evolution

### Week 10 - Functional genomics

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

Mohr et al. 2017. Improved de novo Genome Assembly: Linked-Read Sequencing
Combined with Optical Mapping Produce a High Quality Mammalian Genome at
Relatively Low Cost. BioRxiv. https://doi.org/10.1101/128348

Wenger et al. 2018. Highly-accurate long-read sequencing improves variant
detection and assembly of a human genome. https://www.biorxiv.org/content/10.1101/519025v2

Jiang et al. 2018. A Hybrid de novo Assembly of the Sea Pansy (Renilla muelleri)
Genome. https://doi.org/10.1101/424614

Vertebrate Genomes Project - https://vgp.github.io/genomeark/




### Genome annotation

### Genome structure

Li et al. 2011. Chromosome Size in Diploid Eukaryotic Species Centers on the
Average Length with a Conserved Boundary. Molecular Biology and Evolution.
https://doi.org/10.1093/molbev/msr011

### Functional genomics

### Reconstructing the history of genome evolution

Kim et al. 2017. Reconstruction and evolutionary history of eutherian chromosomes. PNAS. https://doi.org/10.1073/pnas.1702012114

Demas et al. 2018. Reconstruction of avian ancestral karyotypes reveals differences in the evolutionary history of macro- and microchromosomes. BMC Genome Biology. https://doi.org/10.1186/s13059-018-1544-8

O'Connor et al. 2018. Reconstruction of the diapsid ancestral genome permits chromosome evolution tracing in avian and non-avian dinosaurs. https://www.nature.com/articles/s41467-018-04267-9


### Making associations between genomes and phenotypes with comparative methods



### Phylogenetic comparative methods

### General computational skills and methods


## Links and other general resources
