#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage:
    part1_VCF.py (--pileup=STR) (--p=STR) [-o STR]
    
Description:
    Apply Variant Calling based on binomial distribution on input pileup
Arguments:
    --pileup STR                  Input pileup file
    --p STR                           Probability of success
    
Options:
    -o STR                            Generate VCF output file

Examples:
    part1_VCF.py --pileup=my_pileup_file.pileup --p=0.8 -o
    part1_VCF.py --pileup=my_pileup_file.pileup --p=0.8
"""



from collections import Counter
from docopt import docopt
import pandas as pd
import math
import pysam
import operator
import re
import sys

class PileupRecord:
    
    def __init__(self,line):
        fields = line.split("	")
        self.seq = fields[0]
        self.pos = int(fields[1])
        self.ref = fields[2]
        self.rCount = int(fields[3])
        self.rRes = fields[4]
        self.qual = fields[5][:-1]
        
    def print(self):
        print('{}    {}    {}    {}    {}    {}'.format(self.seq, self.pos, self.ref, self.rCount, self.rRes, self.qual))


def VCFHeader(sample):
    
    '''Creates VariantHeader object for the vcf file containing variant calling results for sample from the input'''
    
    VCFheader = pysam.VariantHeader()
    VCFheader.add_sample(sample)
    faifile = open("human_g1k_v37_decoy.fasta.fai")
    for line in faifile:
            split_line = line.split("\t")
            contig = '##contig=<ID=' + str(split_line[0]) + ', length=' + str(split_line[1]) + '>'
            VCFheader.add_line(contig)
    faifile.close()
    
    VCFheader.add_line("##ALT=<ID=*,Description=Represents allele(s) other than observed.>")
    VCFheader.add_line("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>")
    VCFheader.add_line("##FORMAT=<ID=VAF,Number=1,Type=String,Description=Called genome probability>")
 
    return VCFheader

def VCFBodyLine(rec, alts, genotype, p, VCF_out, sample):
    
    '''input: rec - line in the pileup file at current position
              alts - alternative alleles at current position
              genotype - genotype called at current position
              p - probability of the called genotype
              VCF_out - VCF output file
              sample - sample
      output: a line ready to be written in the vcf file '''
    
    record = VCF_out.header.new_record()
    record.contig = rec.seq
    record.ref = rec.ref
    record.alts = alts
    record.pos = rec.pos
    record.samples[sample]['GT'] = genotype
    record.samples[sample]['VAF'] = str(p)
    
    return record

def DetermineIndelString(readResult):
    
    ind_num = [ ind.end() for ind in (re.finditer(r'[A-Za-z\.,][0-9]*[^A-Za-z]', readResult))]
    if len(ind_num) > 0:
        string = readResult[2:ind_num[0]]
        if len(ind_num) != 1:
            for i in range(0,len(ind_num)-1):
                string += int(float(readResult[ind_num[i]]))*readResult[ind_num[i]+1:ind_num[i+1]]
        string += int(float(readResult[ind_num[-1]]))*readResult[ind_num[-1]+1:]
    else:
        string = readResult[2:]
        
    return string


def DetectPolymorphicSite(readResult):
    
    '''Detects variants at particular location, returns counts and types of them.
       input: string, info about particular location alignment results from all reads
       output: dataframe, columns: variant, count, type'''
    
    readResult = readResult.upper().replace(',','.')

    irrelevant = list(set(re.findall(r'\^[^\.]', readResult)))
    for s in irrelevant:
        readResult = readResult.replace(s,'')
    
    occ = re.findall(r'[\.][+-][ACGT]*[0-9]*[ACGT]*[0-9]*[ACGT]*', readResult)
    var = []
    variants = list(set(occ))
    variants1 = variants
    variants1 = [DetermineIndelString(variant) for variant in variants]
    varCounts = [occ.count(indel) for indel in variants]
    varTypes = ['indel']*len(variants1)
    for i in range(0,len(variants1)):
        var.append([variants1[i], varCounts[i], varTypes[i]])
    
    variants.sort(key = len, reverse = True)
    for s in variants:
        readResult = readResult.replace(s,'')
    
    occ = re.findall(r'[AGCT]', readResult)
    SNVs = list(set(occ))
    SNVCounts = [occ.count(SNV) for SNV in SNVs]
    SNVTypes = ['SNV']*len(SNVs)
    for i in range(0,len(SNVs)):
        var.append([SNVs[i], SNVCounts[i], SNVTypes[i]])

    matchCount = len(re.findall(r'[\.]', readResult))
    var.append(['.', matchCount, 'match'])
    
    var.sort(key = lambda x: x[1], reverse = True)

    if len(var) > 2:
        return var[:2]
    else:
        return var 


def Genotyping(var, p):
    
    '''Determines genotype.'''
    
    if len(var) == 1:
        if var[0][2] == 'match':
            genotype = (0,0)
        else:
            genotype = (1,1)
        P = [1]
    if len(var) == 2:
        pr = [0]*3
        k1 = var[0][1]
        k2 = var[1][1]
        # k1 = 40
        # k2 = 1
        if var[0][2] == 'match' or var[0][2] == 'SNV':
            p0 = 0.8
        else:
            p0 = 0.6
            
        if var[1][2] == 'match' or var[1][2] == 'SNV':
            p1 = 0.8
        else:
            p1 = 0.6
        p2 = (p0 + p1)/2

        pr[0] = math.factorial(k1+k2)//math.factorial(k1)//math.factorial(k2)*(p0**k1)*(1-p0)**(k2) # a1a1
        pr[1] = math.factorial(k1+k2)//math.factorial(k1+k2)//math.factorial(0)*(p2**(k1+k2))*(1-p2)**0 # a1a2
        pr[2] = math.factorial(k1+k2)//math.factorial(k2)//math.factorial(k1)*p1**k2*(1-p1)**(k1) # a2a2
    
        P = [pr[0]/(pr[0]+pr[1]+pr[2]), pr[1]/(pr[0]+pr[1]+pr[2]), pr[2]/(pr[0]+pr[1]+pr[2])]
        index, value = max(enumerate(P), key=operator.itemgetter(1))
    
        if var[0][2] == 'match': 
            if index == 0:
                genotype = (0,0)
            elif index == 1:
                genotype = (0,1)
            else:
                genotype =(1,1)
        elif var[1][2] == 'match':
            if index == 0:
                genotype = (1,1)
            elif index == 1:
                genotype = (0,1)
            else:
                genotype = (0,0)
            P.reverse()
        else:
            if index == 0:
                genotype = (1,1)
            elif index == 1:
                genotype = (1,2)
            else:
                genotype = (2,2)
  
    return genotype, P


def DetermineAltsField(polymorphic_site):
    
    polymorphic_site = [elem for elem in polymorphic_site if(elem[0] != '.')]
    if len(polymorphic_site) == 0:
        alts = ['.']
    else:
        alts = [p[0] for p in polymorphic_site]
    
    return alts


def AvgQual(qual):
    qual_prob = []
    for char in qual:
        qual_prob.append(1 - 10**-((ord(char)-33)/10.0))
    
    return sum(qual_prob)/len(qual_prob)

#%%
    
args = docopt(__doc__)
pileupp = args['--pileup']
pr = args['--p']
out = args['-o']

pileup_file = open(pileupp)
for line in pileup_file:
    print(line)
    break
sampleID = 'sample1'

VCF_header = VCFHeader(sampleID)
if out:
    VCF_out = pysam.VariantFile(out, "w", header = VCF_header)
else:
    VCF_out = pysam.VariantFile('-', 'w', header = VCF_header)


i = 0
for line in pileup_file:
    i+=1
    rec = PileupRecord(line)
    
    polymorphic_site = DetectPolymorphicSite(rec.rRes)
    genotype, probs = Genotyping(polymorphic_site, float(pr))
    
    if genotype != (0,0):
        alts = DetermineAltsField(polymorphic_site)
        VCF_line = VCFBodyLine(rec, alts, genotype, max(probs), VCF_out, sampleID)
        VCF_out.write(VCF_line)

      
VCF_out.close()
if out == False:
    sys.stdout.write(VCF_out)







