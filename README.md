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









 

