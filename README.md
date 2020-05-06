# Variant calling algorithm based on deciding on binomial distribution

> The goal of this project was implementation of variant calling algorithm based on deciding on binomial distribution in python and its results analysis. 

## What is Variant Calling?

Variant calling is a very important workflow in genomic variation studies. It includes identification of variants associated with a specific trait or population. Simply put - it is a process of finding differences between reference genome and sequence data. There are three main steps in the analyses that aims at variant calling. These steps are shown and briefly explained in the following figure. It is important to emphasize that these are not the only steps, but are the essential and the most important ones.

<p align="center">
  <img src="images/intro_diag.png">
</p>

The input of variant calling algorithm is a pileup file. Pileup file is a text file that summarizes the base calls of aligned reads to a reference sequence. A couple of lines in pileup file is shown below. 

<p align="center">
  <img src="images/pileup_lines.png">
</p>

### Deciding on binomial distribution

