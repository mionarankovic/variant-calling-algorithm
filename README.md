# Variant calling algorithm based on deciding on binomial distribution

> The goal of this project was implementation of variant calling algorithm based on deciding on binomial distribution in python and its results analysis. 

## What is Variant Calling?

<p> Variant calling is a very important workflow in genomic variation studies. It includes identification of variants associated with a specific trait or population. Simply put - it is a process of finding differences between reference genome and sequence data. There are three main steps in the analyses that aims at variant calling. These steps are shown and briefly explained in the following figure. It is important to emphasize that these are not the only steps, but are the essential and the most important ones. </p>


<p align="center">
  <img src="images/intro_diag.png">
</p>

The input of variant calling algorithm is a pileup file. Pileup file is a text file that summarizes the base calls of aligned reads to a reference sequence. A couple of lines in pileup file is shown below. Each line summarize one particular position in the genome. (Find out more 


<p align="center">
  <img src="images/pileup_lines.png">
</p>

There will be one of two possible cases at each position in the pileup file:
  1. All bases will be the same.
  2. In the simplest case there will be two different bases, as there can be only two alleles at a site. 

The first one represents the scenario in which we, as a diploid organism, got the same base both from mother and from father. That base can be the same as the one in the reference genome, when the reference homozygous is called (0/0), which is the most common scenario. This base can differ from the reference base, when genotype is 1/1. In both scenarios it is clear what base should be called - the one that occurred. 


The other case needs some further calculation in order to make a decision. For this purpose we will use formula for binomial distribution 

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=f=\binom{n}{k}p^{k}(1-p)^{n-k}" \Large>
</p>

where n represents total number of bases sequenced at the particular position, p is the probability of success (the probability that the pileup shows just what really is located at the particular position in the DNA) and k represents number of times we got the exact base we want to calculate f for.
After computing probabilities of all possible explanations of pileup, we pick the most probable one.

The output of the variant calling algorithm is VCF 4.2 file.

## Data

As mentioned previously, the algorithm takes pileup file as an input and produces the output VCF 4.2 file. 

For the purposes of this project, pileup file is generated with samtools pileup tool. This tool takes both reference and bam file as an input to produce the output in pileup format. 

The bam file is created from the exome portions of chromosomes 21, 22, Y and MT and is realigned with the latest version of human genome compatible with b37 by GATK. The reference fasta file is used in the 1000 genomes phase 3 project. It contains chromosomes (1-22, X, Y and MT), unlocalized and unplaced sequences (GL contigs), Epstein–Barr virus (NC_007605), decoy sequences (hs37d5, consisting of BAC+Huref+NA12878+Fosmid clones).

The test VCF file is created with bcftools call tool.

## Results

The first choice is p = 0.9. The confusion matrix is shown below. 

<p align="left">
  <img src="images/Confusion_matrix_09.png" height = 250 width = 450>
</p>

According to metrics written next to the confusion matrix, these results seem to be very bad - precision is unacceptably low. What caused this is a large number of false positives, variants that are present in the implemented solution, but not in the bcftools call tool. This choice is good for understanding the decision making though. 

Variants are shown in the figure below (note that not all variants are shown in the figure, for the purpose of making it interpretable).

<p align="left">
  <img src="images/grafikp09.png">
</p>

As it can be seen from the figure, decision making is pretty pessimistic in this case - no matter how large the number of the first allele appearance during the sequencing process is, only one appearance of a different allele will cause the algorithm to call heterozygous variant. This is because the algorithm is pretty much sure there is no error in the previous steps of the analysis, so what is written in the pileupfile is what really occurs in the genome.

It is important to note that how you interprete these metrics strongly depends on what your survey goal is. It is sometimes much more important to detect all positives, so number of false negatives (missed detection) has to be as small as posible. In this case, it is probably less important not to have false positives (false alarm).

The second choice is p = 0.8. It seems to be a pessimistic choice, as the quality strings in the pileup file indicate that base call accuracies go up to 99.9% for each base. Confusion matrix is shown below. 

<p align="left">
  <img src="images/Conf_matrix_p08.png" height = 250 width = 450>
</p>

These results are much better than the previous one. True predictions went up, while false predictions decreased. As false positives decreased twenty-fold, precision is almost ten times larger. 

<p align="left">
  <img src="images/grafikp08.png">
</p>









 

