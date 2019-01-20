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

#### Agenda and notes

- Sequencing technology and instruments

- Genome sequencing
	- Challenging factors
		- Large size
		- Repeats
		- Heterozygosity

#### Exercises

We will walk through regular expressions in the exercises at
https://github.com/Yale-EEB723/syllabus/blob/master/regular_expressions.txt .

### Week 3, January 31 - Practical computing skills

#### Reading

*Practical Computing for Biologists* Chapters XX

git chapter draft (to be provided as hard copy in previous week)


#### Agenda and notes

- Quick (<1 minute) description of final project plan

- Working on your laptop
- Working on the cluster via your account
- Getting started with git

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

### Genome sequencing, assembly, and annotation

Alkan et al. 2011. Genome structural variation discovery and genotyping. Nature Reviews Genetics. https://www.nature.com/articles/nrg2958

Goodwin et al. 2016. Coming of age: ten years of next-generation sequencing
technologies. Nature Reviews Genetics.  https://doi.org/10.1038/nrg.2016.49

Heather et al. 2015. The sequence of sequencers: The history of sequencing DNA.
Molecular Cell. https://doi.org/10.1016/j.molcel.2015.05.004

Reuter et al. 2015. High-Throughput Sequencing Technologies. Molecular Cell.
https://doi.org/10.1016/j.molcel.2015.05.004

Shendure et al. 2017. DNA sequencing at 40: past, present and future. Nature.
https://doi.org/10.1038/nature24286

### Functional genomics

### Reconstructing the history of genome evolution

Kim et al. 2017. Reconstruction and evolutionary history of eutherian chromosomes. PNAS. https://doi.org/10.1073/pnas.1702012114

Demas et al. 2018. Reconstruction of avian ancestral karyotypes reveals differences in the evolutionary history of macro- and microchromosomes. BMC Genome Biology. https://doi.org/10.1186/s13059-018-1544-8

O'Connor et al. 2018. Reconstruction of the diapsid ancestral genome permits chromosome evolution tracing in avian and non-avian dinosaurs. https://www.nature.com/articles/s41467-018-04267-9


### Making associations between genomes and phenotypes with comparative methods



### Phylogenetic comparative methods

### General computational skills and methods


## Links and other general resources
