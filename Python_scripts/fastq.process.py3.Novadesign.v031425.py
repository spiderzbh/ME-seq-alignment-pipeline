#!/usr/bin/env python

# Author: Sai Ma
# The following program will process the fastq
# handle QC vs full & modify fastq header & split fastqs & add index & trim & split project
# example python3 /mnt/users/sai/Script/Split-seq_Sai/fastq.process.py3.py -a Undetermined_S0_R1_001.fastq.gz -b Undetermined_S0_R2_001.fastq.gz --qc -y /mnt/users/sai/Script/Split-seq_Sai/config_test.yaml

# to do list
## need to add N6 compatability

##### IMPORT MODULES #####
# import necessary for python
import os
import re
import sys
import bz2
import gzip
import string
import Levenshtein
import json
import yaml
from Bio import SeqIO
from Bio import AlignIO
from optparse import OptionParser
import time
import io

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = str.maketrans('ATCGN', 'TAGCN')

def reverse_complement(sequence):
    return sequence.decode("utf-8").upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

def barcodeSet(barcode):
    bases = "ATCGN"
    barcodeSet = set()
    barcodeSet.add(barcode)
    for i, c in enumerate(barcode):
        if c in bases:
            for base in bases:
                if c != base:
                    barcodeSet.add((barcode[:i] + base + barcode[i + 1:]))
                    ## allow 1 mismatch or 1 bp shift
                    # barcodeSet.add((barcode[1:] + base))
                    # barcodeSet.add((base + barcode[:-1]))
    return barcodeSet
        
## next-seq, nova-seq 1.5 P5 barcode
P5revcomp = {'GCGATCTA': 'P1.01',
      'ATAGAGAG': 'P1.02',
      'AGAGGATA': 'P1.03',
      'TCTACTCT': 'P1.04',
      'CTCCTTAC': 'P1.05',
      'TATGCAGT': 'P1.06',
      'TACTCCTT': 'P1.07',
      'AGGCTTAG': 'P1.08',
      'GATTTCCA': 'P1.09',
      'ATCATGTT': 'P1.10',
      'TTTCATCA': 'P1.11',
      'AGTCCGAC': 'P1.12',
      'GCTAGAAA': 'P1.13',
      'CTTGGTTA': 'P1.14',
      'CGATACAC': 'P1.15',
      'TTGATGGA': 'P1.16',
      'TGCACGAA': 'P1.17',
      'GGCAACCT': 'P1.18',
      'ACATAAGG': 'P1.19',
      'CGTTGCTG': 'P1.20',
      'ATTGAACC': 'P1.21',
      'ACGAATGT': 'P1.22',
      'TGGGAATC': 'P1.23',
      'GCAGTCCG': 'P1.24',
      'GAACGGCT': 'P1.25',
      'GACCCAAT': 'P1.26',
      'AGTATGCA': 'P1.27',
      'CCAAGCCC': 'P1.28',
      'GCCACGTC': 'P1.29',
      'AAATTTGC': 'P1.30',
      'GAGGCTGC': 'P1.31',
      'AACTCGGA': 'P1.32',
      'CTTAATGC': 'P1.33',
      'GTTATCGT': 'P1.34',
      'CCCGCAGG': 'P1.35',
      'AACAATCA': 'P1.36',
      'TCCGTGCC': 'P1.37',
      'GAATGATC': 'P1.38',
      'ATGACCAT': 'P1.39',
      'TTGGTACG': 'P1.40',
      'TAAACTGG': 'P1.41',
      'GGGCCGGT': 'P1.42',
      'ACTTCTAG': 'P1.43',
      'ATCTGGCG': 'P1.44',
      'CCATGTGA': 'P1.45',
      'TCGAGTTC': 'P1.46',
      'AACGGTGG': 'P1.47',
      'GTAACTTA': 'P1.48',
      'CACGTCTC': 'P1.49',
      'TTAGGCAA': 'P1.50',
      'CAAGTTAA': 'P1.51',
      'TGTTAAAG': 'P1.52',
      'GGTCTACG': 'P1.53',
      'CGCAAATA': 'P1.54',
      'TCCTGGAT': 'P1.55',
      'CAGGAACA': 'P1.56',
      'CTGCGCGT': 'P1.57',
      'TCGCCAGA': 'P1.58',
      'TGTAGATT': 'P1.59',
      'GGTCAGTA': 'P1.60',
      'CCCTATCG': 'P1.61',
      'TTCTAAGT': 'P1.62',
      'AGATCTCT': 'P1.63',
      'CCTTCACC': 'P1.64',
      'CATTCGAT': 'P1.65',
      'GCTCTTGA': 'P1.66',
      'ACGTGGGC': 'P1.67',
      'ACCGCCCA': 'P1.68',
      'TCCAAGGG': 'P1.69',
      'ACGGTAAT': 'P1.70',
      'CTCGGACT': 'P1.71',
      'CAACAAGT': 'P1.72',
      'TGTATTAC': 'P1.73',
      'TAGACGCC': 'P1.74',
      'AGCAGCGC': 'P1.75',
      'AATGGCAC': 'P1.76',
      'CATACCTA': 'P1.77',
      'TAGGTGTT': 'P1.78',
      'GTTCGGAG': 'P1.79',
      'TGCCGTTG': 'P1.80',
      'CTACATTG': 'P1.81',
      'GGGTAGCC': 'P1.82',
      'CGGACTTT': 'P1.83',
      'CCGCGGAA': 'P1.84',
      'AAGTGCCT': 'P1.85',
      'CACTGAAG': 'P1.86',
      'CTACCGGC': 'P1.87',
      'GGATTGAA': 'P1.88',
      'GTGTGTGG': 'P1.89',
      'GATAATAT': 'P1.90',
      'TGCTTCGG': 'P1.91',
      'ACCGATAC': 'P1.92'}

## mi-seq, nova-seq 1.0 P5 barcode
P5fwdstr = {'TAGATCGC': 'P1.01',
      'CTCTCTAT': 'P1.02',
      'TATCCTCT': 'P1.03',
      'AGAGTAGA': 'P1.04',
      'GTAAGGAG': 'P1.05',
      'ACTGCATA': 'P1.06',
      'AAGGAGTA': 'P1.07',
      'CTAAGCCT': 'P1.08',
      'TGGAAATC': 'P1.09',
      'AACATGAT': 'P1.10',
      'TGATGAAA': 'P1.11',
      'GTCGGACT': 'P1.12',
      'TTTCTAGC': 'P1.13',
      'TAACCAAG': 'P1.14',
      'GTGTATCG': 'P1.15',
      'TCCATCAA': 'P1.16',
      'TTCGTGCA': 'P1.17',
      'AGGTTGCC': 'P1.18',
      'CCTTATGT': 'P1.19',
      'CAGCAACG': 'P1.20',
      'GGTTCAAT': 'P1.21',
      'ACATTCGT': 'P1.22',
      'GATTCCCA': 'P1.23',
      'CGGACTGC': 'P1.24',
      'AGCCGTTC': 'P1.25',
      'ATTGGGTC': 'P1.26',
      'TGCATACT': 'P1.27',
      'GGGCTTGG': 'P1.28',
      'GACGTGGC': 'P1.29',
      'GCAAATTT': 'P1.30',
      'GCAGCCTC': 'P1.31',
      'TCCGAGTT': 'P1.32',
      'GCATTAAG': 'P1.33',
      'ACGATAAC': 'P1.34',
      'CCTGCGGG': 'P1.35',
      'TGATTGTT': 'P1.36',
      'GGCACGGA': 'P1.37',
      'GATCATTC': 'P1.38',
      'ATGGTCAT': 'P1.39',
      'CGTACCAA': 'P1.40',
      'CCAGTTTA': 'P1.41',
      'ACCGGCCC': 'P1.42',
      'CTAGAAGT': 'P1.43',
      'CGCCAGAT': 'P1.44',
      'TCACATGG': 'P1.45',
      'GAACTCGA': 'P1.46',
      'CCACCGTT': 'P1.47',
      'TAAGTTAC': 'P1.48',
      'GAGACGTG': 'P1.49',
      'TTGCCTAA': 'P1.50',
      'TTAACTTG': 'P1.51',
      'CTTTAACA': 'P1.52',
      'CGTAGACC': 'P1.53',
      'TATTTGCG': 'P1.54',
      'ATCCAGGA': 'P1.55',
      'TGTTCCTG': 'P1.56',
      'ACGCGCAG': 'P1.57',
      'TCTGGCGA': 'P1.58',
      'AATCTACA': 'P1.59',
      'TACTGACC': 'P1.60',
      'CGATAGGG': 'P1.61',
      'ACTTAGAA': 'P1.62',
      'AGAGATCT': 'P1.63',
      'GGTGAAGG': 'P1.64',
      'ATCGAATG': 'P1.65',
      'TCAAGAGC': 'P1.66',
      'GCCCACGT': 'P1.67',
      'TGGGCGGT': 'P1.68',
      'CCCTTGGA': 'P1.69',
      'ATTACCGT': 'P1.70',
      'AGTCCGAG': 'P1.71',
      'ACTTGTTG': 'P1.72',
      'GTAATACA': 'P1.73',
      'GGCGTCTA': 'P1.74',
      'GCGCTGCT': 'P1.75',
      'GTGCCATT': 'P1.76',
      'TAGGTATG': 'P1.77',
      'AACACCTA': 'P1.78',
      'CTCCGAAC': 'P1.79',
      'CAACGGCA': 'P1.80',
      'CAATGTAG': 'P1.81',
      'GGCTACCC': 'P1.82',
      'AAAGTCCG': 'P1.83',
      'TTCCGCGG': 'P1.84',
      'AGGCACTT': 'P1.85',
      'CTTCAGTG': 'P1.86',
      'GCCGGTAG': 'P1.87',
      'TTCAATCC': 'P1.88',
      'CCACACAC': 'P1.89',
      'ATATTATC': 'P1.90',
      'CCGAAGCA': 'P1.91',
      'GTATCGGT': 'P1.92'}
P7={'TAAGGCGA': 'P2.01',
'CGTACTAG': 'P2.02',
'AGGCAGAA': 'P2.03',
'TCCTGAGC': 'P2.04',
'GGACTCCT': 'P2.05',
'TAGGCATG': 'P2.06',
'CTCTCTAC': 'P2.07',
'CAGAGAGG': 'P2.08',
'GCTACGCT': 'P2.09',
'CGAGGCTG': 'P2.10',
'AAGAGGCA': 'P2.11',
'GTAGAGGA': 'P2.12',
'TGGATCTG': 'P2.13',
'CCGTTTGT': 'P2.14',
'TGCTGGGT': 'P2.15',
'AGGTTGGG': 'P2.16',
'GTGTGGTG': 'P2.17',
'TGGGTTTC': 'P2.18',
'TGGTCACA': 'P2.19',
'TTGACCCT': 'P2.20',
'CGCGGACA': 'P2.21',
'TTCCATAT': 'P2.22',
'AATTCGTT': 'P2.23',
'GGCGTCGA': 'P2.24',
'ACAAAGTG': 'P2.25',
'TACTTGAA': 'P2.26',
'GTGATAGC': 'P2.27',
'AGTAGATT': 'P2.28',
'ATTGCCGG': 'P2.29',
'TTGCTAAG': 'P2.30',
'ATAAGTTA': 'P2.31',
'ATCACTCG': 'P2.32',
'GTTAACAG': 'P2.33',
'AATGGTAG': 'P2.34',
'GAGCACGT': 'P2.35',
'TTTCGTCA': 'P2.36',
'CAAGAATT': 'P2.37',
'GAAATGCC': 'P2.38',
'AACGCCAT': 'P2.39',
'CCTCGCAG': 'P2.40',
'TACACCTC': 'P2.41',
'GGTCATTT': 'P2.42',
'CAATCTTA': 'P2.43',
'TGTGCCTT': 'P2.44',
'TCTTATTA': 'P2.45',
'GACTTAGT': 'P2.46',
'AGACCAGC': 'P2.47',
'AAATACAG': 'P2.48',
'TTATGAAA': 'P2.49',
'CTTGGGTC': 'P2.50',
'CCAAATAA': 'P2.51',
'GCGTTAAA': 'P2.52',
'CATCCTGT': 'P2.53',
'GGAGTAAG': 'P2.54',
'GACGCTCC': 'P2.55',
'TTCGCGGC': 'P2.56',
'CGGTTCCC': 'P2.57',
'ACCGGCTA': 'P2.58',
'CTCATGGG': 'P2.59',
'TTTAATGC': 'P2.60',
'AAACGGTC': 'P2.61',
'GATCCAAA': 'P2.62',
'ATGATGAT': 'P2.63',
'CCAACACG': 'P2.64',
'TAACAACA': 'P2.65',
'GGTAAACC': 'P2.66',
'CATCGACC': 'P2.67',
'ATGGGAAC': 'P2.68',
'CGGCCAAT': 'P2.69',
'GGGAATGA': 'P2.70',
'GTATTCGG': 'P2.71',
'TCAGCTAT': 'P2.72',
'ATTTATCT': 'P2.73',
'ACAGTTGC': 'P2.74',
'CCCGAGAT': 'P2.75',
'TAATGTCT': 'P2.76',
'GCCAATTC': 'P2.77',
'CGCCGTGC': 'P2.78',
'CTGACCGA': 'P2.79',
'CATTTCGA': 'P2.80',
'GCTTGCCA': 'P2.81',
'TTCTACCA': 'P2.82',
'ACGTGACG': 'P2.83',
'TGTCCGCG': 'P2.84',
'TTAAACTT': 'P2.85',
'ACCACAAC': 'P2.86',
'GCCTCTGG': 'P2.87',
'TCGCCCAC': 'P2.88',
'CACTAGGC': 'P2.89',
'TCGAAGCC': 'P2.90',
'GCATGTAC': 'P2.91',
'GTTCGAGT': 'P2.92',
'CCGGGCGC': 'P2.93',
'AGATTTAA': 'P2.94',
'CACCATTG': 'P2.95',
'AATAAGAC': 'P2.96',
'AAGTAGAG': 'P2.97'}

R1 = {'AACGTGAT': 'R1.001',
'AAACATCG': 'R1.002',
'ATGCCTAA': 'R1.003',
'AGTGGTCA': 'R1.004',
'ACCACTGT': 'R1.005',
'ACATTGGC': 'R1.006',
'CAGATCTG': 'R1.007',
'CATCAAGT': 'R1.008',
'CGCTGATC': 'R1.009',
'ACAAGCTA': 'R1.010',
'CTGTAGCC': 'R1.011',
'AGTACAAG': 'R1.012',
'AACAACCA': 'R1.013',
'AACCGAGA': 'R1.014',
'AACGCTTA': 'R1.015',
'AAGACGGA': 'R1.016',
'AAGGTACA': 'R1.017',
'ACACAGAA': 'R1.018',
'ACAGCAGA': 'R1.019',
'ACCTCCAA': 'R1.020',
'ACGCTCGA': 'R1.021',
'ACGTATCA': 'R1.022',
'ACTATGCA': 'R1.023',
'AGAGTCAA': 'R1.024',
'AGATCGCA': 'R1.025',
'AGCAGGAA': 'R1.026',
'AGTCACTA': 'R1.027',
'ATCCTGTA': 'R1.028',
'ATTGAGGA': 'R1.029',
'CAACCACA': 'R1.030',
'GACTAGTA': 'R1.031',
'CAATGGAA': 'R1.032',
'CACTTCGA': 'R1.033',
'CAGCGTTA': 'R1.034',
'CATACCAA': 'R1.035',
'CCAGTTCA': 'R1.036',
'CCGAAGTA': 'R1.037',
'CCGTGAGA': 'R1.038',
'CCTCCTGA': 'R1.039',
'CGAACTTA': 'R1.040',
'CGACTGGA': 'R1.041',
'CGCATACA': 'R1.042',
'CTCAATGA': 'R1.043',
'CTGAGCCA': 'R1.044',
'CTGGCATA': 'R1.045',
'GAATCTGA': 'R1.046',
'CAAGACTA': 'R1.047',
'GAGCTGAA': 'R1.048',
'GATAGACA': 'R1.049',
'GCCACATA': 'R1.050',
'GCGAGTAA': 'R1.051',
'GCTAACGA': 'R1.052',
'GCTCGGTA': 'R1.053',
'GGAGAACA': 'R1.054',
'GGTGCGAA': 'R1.055',
'GTACGCAA': 'R1.056',
'GTCGTAGA': 'R1.057',
'GTCTGTCA': 'R1.058',
'GTGTTCTA': 'R1.059',
'TAGGATGA': 'R1.060',
'TATCAGCA': 'R1.061',
'TCCGTCTA': 'R1.062',
'TCTTCACA': 'R1.063',
'TGAAGAGA': 'R1.064',
'TGGAACAA': 'R1.065',
'TGGCTTCA': 'R1.066',
'TGGTGGTA': 'R1.067',
'TTCACGCA': 'R1.068',
'AACTCACC': 'R1.069',
'AAGAGATC': 'R1.070',
'AAGGACAC': 'R1.071',
'AATCCGTC': 'R1.072',
'AATGTTGC': 'R1.073',
'ACACGACC': 'R1.074',
'ACAGATTC': 'R1.075',
'AGATGTAC': 'R1.076',
'AGCACCTC': 'R1.077',
'AGCCATGC': 'R1.078',
'AGGCTAAC': 'R1.079',
'ATAGCGAC': 'R1.080',
'ATCATTCC': 'R1.081',
'ATTGGCTC': 'R1.082',
'CAAGGAGC': 'R1.083',
'CACCTTAC': 'R1.084',
'CCATCCTC': 'R1.085',
'CCGACAAC': 'R1.086',
'CCTAATCC': 'R1.087',
'CCTCTATC': 'R1.088',
'CGACACAC': 'R1.089',
'CGGATTGC': 'R1.090',
'CTAAGGTC': 'R1.091',
'GAACAGGC': 'R1.092',
'GACAGTGC': 'R1.093',
'GAGTTAGC': 'R1.094',
'GATGAATC': 'R1.095',
'GCCAAGAC': 'R1.096'}
R2 = {'AACGTGAT': 'R2.001',
'AAACATCG': 'R2.002',
'ATGCCTAA': 'R2.003',
'AGTGGTCA': 'R2.004',
'ACCACTGT': 'R2.005',
'ACATTGGC': 'R2.006',
'CAGATCTG': 'R2.007',
'CATCAAGT': 'R2.008',
'CGCTGATC': 'R2.009',
'ACAAGCTA': 'R2.010',
'CTGTAGCC': 'R2.011',
'AGTACAAG': 'R2.012',
'AACAACCA': 'R2.013',
'AACCGAGA': 'R2.014',
'AACGCTTA': 'R2.015',
'AAGACGGA': 'R2.016',
'AAGGTACA': 'R2.017',
'ACACAGAA': 'R2.018',
'ACAGCAGA': 'R2.019',
'ACCTCCAA': 'R2.020',
'ACGCTCGA': 'R2.021',
'ACGTATCA': 'R2.022',
'ACTATGCA': 'R2.023',
'AGAGTCAA': 'R2.024',
'AGATCGCA': 'R2.025',
'AGCAGGAA': 'R2.026',
'AGTCACTA': 'R2.027',
'ATCCTGTA': 'R2.028',
'ATTGAGGA': 'R2.029',
'CAACCACA': 'R2.030',
'GACTAGTA': 'R2.031',
'CAATGGAA': 'R2.032',
'CACTTCGA': 'R2.033',
'CAGCGTTA': 'R2.034',
'CATACCAA': 'R2.035',
'CCAGTTCA': 'R2.036',
'CCGAAGTA': 'R2.037',
'CCGTGAGA': 'R2.038',
'CCTCCTGA': 'R2.039',
'CGAACTTA': 'R2.040',
'CGACTGGA': 'R2.041',
'CGCATACA': 'R2.042',
'CTCAATGA': 'R2.043',
'CTGAGCCA': 'R2.044',
'CTGGCATA': 'R2.045',
'GAATCTGA': 'R2.046',
'CAAGACTA': 'R2.047',
'GAGCTGAA': 'R2.048',
'GATAGACA': 'R2.049',
'GCCACATA': 'R2.050',
'GCGAGTAA': 'R2.051',
'GCTAACGA': 'R2.052',
'GCTCGGTA': 'R2.053',
'GGAGAACA': 'R2.054',
'GGTGCGAA': 'R2.055',
'GTACGCAA': 'R2.056',
'GTCGTAGA': 'R2.057',
'GTCTGTCA': 'R2.058',
'GTGTTCTA': 'R2.059',
'TAGGATGA': 'R2.060',
'TATCAGCA': 'R2.061',
'TCCGTCTA': 'R2.062',
'TCTTCACA': 'R2.063',
'TGAAGAGA': 'R2.064',
'TGGAACAA': 'R2.065',
'TGGCTTCA': 'R2.066',
'TGGTGGTA': 'R2.067',
'TTCACGCA': 'R2.068',
'AACTCACC': 'R2.069',
'AAGAGATC': 'R2.070',
'AAGGACAC': 'R2.071',
'AATCCGTC': 'R2.072',
'AATGTTGC': 'R2.073',
'ACACGACC': 'R2.074',
'ACAGATTC': 'R2.075',
'AGATGTAC': 'R2.076',
'AGCACCTC': 'R2.077',
'AGCCATGC': 'R2.078',
'AGGCTAAC': 'R2.079',
'ATAGCGAC': 'R2.080',
'ATCATTCC': 'R2.081',
'ATTGGCTC': 'R2.082',
'CAAGGAGC': 'R2.083',
'CACCTTAC': 'R2.084',
'CCATCCTC': 'R2.085',
'CCGACAAC': 'R2.086',
'CCTAATCC': 'R2.087',
'CCTCTATC': 'R2.088',
'CGACACAC': 'R2.089',
'CGGATTGC': 'R2.090',
'CTAAGGTC': 'R2.091',
'GAACAGGC': 'R2.092',
'GACAGTGC': 'R2.093',
'GAGTTAGC': 'R2.094',
'GATGAATC': 'R2.095',
'GCCAAGAC': 'R2.096'}
R3 = {'AACGTGAT': 'R3.001',
'AAACATCG': 'R3.002',
'ATGCCTAA': 'R3.003',
'AGTGGTCA': 'R3.004',
'ACCACTGT': 'R3.005',
'ACATTGGC': 'R3.006',
'CAGATCTG': 'R3.007',
'CATCAAGT': 'R3.008',
'CGCTGATC': 'R3.009',
'ACAAGCTA': 'R3.010',
'CTGTAGCC': 'R3.011',
'AGTACAAG': 'R3.012',
'AACAACCA': 'R3.013',
'AACCGAGA': 'R3.014',
'AACGCTTA': 'R3.015',
'AAGACGGA': 'R3.016',
'AAGGTACA': 'R3.017',
'ACACAGAA': 'R3.018',
'ACAGCAGA': 'R3.019',
'ACCTCCAA': 'R3.020',
'ACGCTCGA': 'R3.021',
'ACGTATCA': 'R3.022',
'ACTATGCA': 'R3.023',
'AGAGTCAA': 'R3.024',
'AGATCGCA': 'R3.025',
'AGCAGGAA': 'R3.026',
'AGTCACTA': 'R3.027',
'ATCCTGTA': 'R3.028',
'ATTGAGGA': 'R3.029',
'CAACCACA': 'R3.030',
'GACTAGTA': 'R3.031',
'CAATGGAA': 'R3.032',
'CACTTCGA': 'R3.033',
'CAGCGTTA': 'R3.034',
'CATACCAA': 'R3.035',
'CCAGTTCA': 'R3.036',
'CCGAAGTA': 'R3.037',
'CCGTGAGA': 'R3.038',
'CCTCCTGA': 'R3.039',
'CGAACTTA': 'R3.040',
'CGACTGGA': 'R3.041',
'CGCATACA': 'R3.042',
'CTCAATGA': 'R3.043',
'CTGAGCCA': 'R3.044',
'CTGGCATA': 'R3.045',
'GAATCTGA': 'R3.046',
'CAAGACTA': 'R3.047',
'GAGCTGAA': 'R3.048',
'GATAGACA': 'R3.049',
'GCCACATA': 'R3.050',
'GCGAGTAA': 'R3.051',
'GCTAACGA': 'R3.052',
'GCTCGGTA': 'R3.053',
'GGAGAACA': 'R3.054',
'GGTGCGAA': 'R3.055',
'GTACGCAA': 'R3.056',
'GTCGTAGA': 'R3.057',
'GTCTGTCA': 'R3.058',
'GTGTTCTA': 'R3.059',
'TAGGATGA': 'R3.060',
'TATCAGCA': 'R3.061',
'TCCGTCTA': 'R3.062',
'TCTTCACA': 'R3.063',
'TGAAGAGA': 'R3.064',
'TGGAACAA': 'R3.065',
'TGGCTTCA': 'R3.066',
'TGGTGGTA': 'R3.067',
'TTCACGCA': 'R3.068',
'AACTCACC': 'R3.069',
'AAGAGATC': 'R3.070',
'AAGGACAC': 'R3.071',
'AATCCGTC': 'R3.072',
'AATGTTGC': 'R3.073',
'ACACGACC': 'R3.074',
'ACAGATTC': 'R3.075',
'AGATGTAC': 'R3.076',
'AGCACCTC': 'R3.077',
'AGCCATGC': 'R3.078',
'AGGCTAAC': 'R3.079',
'ATAGCGAC': 'R3.080',
'ATCATTCC': 'R3.081',
'ATTGGCTC': 'R3.082',
'CAAGGAGC': 'R3.083',
'CACCTTAC': 'R3.084',
'CCATCCTC': 'R3.085',
'CCGACAAC': 'R3.086',
'CCTAATCC': 'R3.087',
'CCTCTATC': 'R3.088',
'CGACACAC': 'R3.089',
'CGGATTGC': 'R3.090',
'CTAAGGTC': 'R3.091',
'GAACAGGC': 'R3.092',
'GACAGTGC': 'R3.093',
'GAGTTAGC': 'R3.094',
'GATGAATC': 'R3.095',
'GCCAAGAC': 'R3.096'}

#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-y", help="<Yaml> yaml file that specifies P5 barcode for each project")
opts.add_option("-a", help="<Read1> Accepts xxx_S1_R1_001.fastq.gz")
opts.add_option("-b", help="<Read2> Accepts xxx_S1_R2_001.fastq.gz")
opts.add_option("--c", help="<Index1> optional xxx_S1_I1_001.fastq.gz", default="NA")
opts.add_option("--d", help="<Index2> optional xxx_S1_I2_001.fastq.gz", default="NA")
#opts.add_option("--outdir", help="output dir")
opts.add_option("--qc", action="store_true", help="QC run with first 3M reads")
opts.add_option("--out", help="Path to the output fastq files")
opts.add_option("-t", help="P5 chemistry accepts rev and fwd")

options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# name input and outputs
p1_in = options.a
p2_in = options.b
i1_in = options.c
i2_in = options.d
yaml_in = options.y
chemistry = options.t
prefix = options.out
qcreads=99999 # number of reads for QC analysis


#if i1_in == "NA":
#    print('2 fqs are supplied')
#else:
#    print('4 fqs are supplied')

#if options.qc == True:
#    print('Running QC on ' + str(qcreads+1)  + " reads")
    
# name outputs and print to working dir
p1_file = p1_in.split('/')[-1]
p2_file = p2_in.split('/')[-1]

#check for file type and open input file
append = p1_in.split('.')[-1]
if append == "gz":
    p1_rds = io.BufferedReader(gzip.open(p1_in,'rb'))
    p2_rds = io.BufferedReader(gzip.open(p2_in,'rb'))
    # set up files for undetermined reads
#    dis1 = io.BufferedWriter(open(prefix + "discard.R1.fq", 'wb'))
#    dis2 = io.BufferedWriter(open(prefix + "discard.R2.fq", 'wb'))
    # read in optional index files
    if i1_in != "NA":
        i1_rds = io.BufferedReader(gzip.open(i1_in,'rb'))
        i2_rds = io.BufferedReader(gzip.open(i2_in,'rb'))
else:
    sys.exit("ERROR! The input file2 must be a .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0;count=1
n=20  # match seq
mismatch=1  # only allow 0-1 mismatches for now
good_r=0

# check if reads are indexed
p1_line = p1_rds.readline()
seqhead1 = p1_line.decode()
# print("first read header:")
# print(seqhead1)

if "+" in seqhead1 or "_" in seqhead1:
    indexed = True
#    print("Fastqs are properly indexed")
else:
    indexed = False
#    print("Will update index")
    if i1_in == "NA":
        sys.exit("warning: one pair of the fastqs are not indexed, but index reads are not supplied")
        

p1_rds.close()
p1_rds = io.BufferedReader(gzip.open(p1_in,'rb'))

# load yaml 
inFile = open(yaml_in, 'r')
config = yaml.load(inFile, Loader=yaml.UnsafeLoader)
projectNames = set()

for key in config.keys():
    if ('Project' in key):
        projectNames.add(key)
# print(projectNames)

# open files to write in
project = dict()
sampletype = dict()
N6type = dict()
printPrimerSet = dict()
for proj in projectNames:
    metaData = config[proj]
    outName = metaData['Name']
    outType = metaData['Type']
    outP5 = metaData['P5']
    if outType == 'RNA2':
        realP5 = metaData['realP5']
    else:
        realP5 = outP5
    outP7 = metaData['P7']
    primer = realP5 + "," + outP7
    outprimer = outP5 + "," + outP7
    # remove existing files
    if os.path.exists(outName + ".R1.fq.gz"):
        os.remove(outName + ".R1.fq.gz")
    if os.path.exists(outName + ".R2.fq.gz"):
        os.remove(outName + ".R2.fq.gz")
      
    project[primer] = outName
    printPrimerSet[primer] = outprimer
    sampletype[primer] = outType
    N6type[primer] = "F"
    if sampletype[primer] == 'RNA' or sampletype[primer] == 'RNA2':
        if 'N6' in metaData.keys():
            N6type[primer] = metaData['N6']
    
#print(project)
#print(sampletype)
#print(N6type)

# generate barcode set
r1set = dict()
r2set = dict()
r3set = dict()
p5set = dict()
p7set = dict()
for barcode, name in R3.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        r3set[bc] = name
for barcode, name in R2.items():
    barcodes = barcodeSet(barcode)
    for bc in barcodes:
        r2set[bc] = name
for barcode, name in R1.items():
    barcodes = barcodeSet(barcode)
#    barcodes2 = barcodeSet(barcode[0:7])
#    barcodes3 = barcodeSet(barcode[0:6])
#    barcodes4 = barcodeSet(barcode[0:5])
    for bc in barcodes:
        r1set[bc] = name
#    for bc in barcodes2:
#        r1set[bc] = name
#    for bc in barcodes3:
#        r1set[bc] = name
#    for bc in barcodes4:
#        r1set[bc] = name
for barcode, name in P7.items():
        barcodes = barcodeSet(barcode)
        for bc in barcodes:
            p7set[bc] = name 
if chemistry == 'rev':
    for barcode, name in P5revcomp.items():
        barcodes = barcodeSet(barcode)
        for bc in barcodes:
            p5set[bc] = name
elif chemistry == 'fwd':            
    for barcode, name in P5fwdstr.items():
        barcodes = barcodeSet(barcode)
        for bc in barcodes:
            p5set[bc] = name
else:
    sys.exit()
files_r1 = dict()
files_r2 = dict()
for proj in project.values():
    f1 = io.BufferedWriter(open((prefix + proj + ".R1.fq"), 'ab'))
    f2 = io.BufferedWriter(open((prefix + proj + ".R2.fq"), 'ab'))
    files_r1[proj] = f1
    files_r2[proj] = f2
#print(files_r1)

start = time.process_time()        
while 1:
    # read lines
    p1_line = p1_rds.readline()
    p2_line = p2_rds.readline()

    # add index to biological reads
    if indexed == False:
        i1_line = i1_rds.readline()
        i2_line = i2_rds.readline()
        
    # break if at end of file
    if not p1_line:
        break

    # load fastq into memory
    if count ==1:
        seqhead1 = p1_line.decode()
        seqhead2 = p2_line.decode()
#        print(seqhead1 + str(i))
        # modify read header
        seqhead1 = seqhead1.replace("1:N:0:1", "1:N:0:")
        seqhead1 = seqhead1.replace("1:N:0:2", "1:N:0:")
        seqhead1 = str.encode(seqhead1.replace("1:N:0:0", "1:N:0:"))
        seqhead2 = seqhead2.replace("4:N:0:1", "2:N:0:")
        seqhead2 = seqhead2.replace("4:N:0:2", "2:N:0:")
        seqhead2 = str.encode(seqhead2.replace("2:N:0:0", "2:N:0:"))        
    elif count ==2:
        seq1 = p1_line.rstrip()
        seq2 = p2_line.rstrip()
        # skip reads in empty
        if seq1.decode() == "+":
            print("warning: found empty biological read. better to use Untrimmed fastq")
            i = i + 1
            p1_line = p1_rds.readline()
            p1_line = p1_rds.readline()
            p2_line = p2_rds.readline()
            p2_line = p2_rds.readline()
            seqhead1 = p1_line.decode()
            seqhead2 = p2_line.decode()
            seqhead1 = seqhead1.replace("1:N:0:1", "1:N:0:")
            seqhead1 = seqhead1.replace("1:N:0:2", "1:N:0:")
            seqhead1 = str.encode(seqhead1.replace("1:N:0:0", "1:N:0:"))
            seqhead2 = seqhead2.replace("4:N:0:1", "2:N:0:")
            seqhead2 = seqhead2.replace("4:N:0:2", "2:N:0:")
            seqhead2 = str.encode(seqhead2.replace("2:N:0:0", "2:N:0:"))
            if indexed == False:
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i1_line = i1_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
                i2_line = i2_rds.readline()
#            print(seqhead1)
#            print(i1_line)

        # print(seq1)
        # print(i1_line)
        if indexed == False:
            id1 = i1_line.rstrip()
            id1 = id1[0:8]
            id2 = i2_line.rstrip()
            id2 = id2[0:8]
            # skip lines in index read is empty
            if "F" in id1.decode():
                print("warning: found empty index read. better to use Untrimmed fastq")
                
            seqhead1 = str.encode(seqhead1.decode().replace("\n", ""))
            seqhead2 = str.encode(seqhead2.decode().replace("\n", ""))
            seqhead1 = (seqhead1 + id1 + b"+" + id2 + b"\n")
            seqhead2 = (seqhead2 + id1 + b"+" + id2 + b"\n")
#        print(id1)    
        # update barcode to R1.xx,R2.xx,R3.xx,P5.xx
        barcode = ""
        barcodeMatch = 0
        index = seqhead1.decode().find("1:N:0:")
        # seqheadMod = seqhead1.decode().replace("_", "+")
        index2 = seqhead1.decode().find("+")
        p7 = seqhead1[(index + 6):(index + 14)].decode()
        p5 = seqhead1[(index + 15):(index + 23)].decode()
        r3 = seq2[0:8].decode()
        r2 = seq2[38:46].decode()
        r1 = seq2[76:84].decode()
        # print(p5)
        # print(r1)
        if (r1 in r1set):
            barcode = r1set[r1]
            barcodeMatch += 1
        if (r2 in r2set):
            barcode = barcode + "," + r2set[r2]
            barcodeMatch += 1
        if (r3 in r3set):
            barcode = barcode + "," + r3set[r3]
            barcodeMatch += 1
        if (p5 in p5set):
            barcode = barcode + "," + p5set[p5]
            barcodeMatch += 1
        if (p7 in p7set):
            barcode = barcode + "," + p7set[p7]
            barcodeMatch += 1
        seqheadMod = seqhead1.decode().replace("_", "+")
        index2 = seqheadMod.find("+")
#        print(barcode)
#        print(barcodeMatch)
    elif count ==3:
        qualhead1 = p1_line
        qualhead2 = p2_line
    elif count ==4:
        qual1 = p1_line.rstrip()
        qual2 = p2_line.rstrip()
        needtrim = "F"
        if (barcodeMatch == 5):
            primerset = p5set[p5] + "," + p7set[p7]
#            print(primerset)
            if primerset in project:
#                print(primerset)
#                print(sampletype[p5set[p5]])
                if sampletype[primerset] == "ATAC" or sampletype[primerset] == "TAPS" or sampletype[primerset] == "DipC":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seq2 = seq2[133: ]
                    qual2 = qual2[133: ]
                   #  needtrim = "T"
                elif sampletype[primerset] == "Purturb":
                    if N6type[primerset] == "F":
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    elif Levenshtein.distance(seq2[124:132].decode(),"TTTTTTTT") < 3:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124]
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124]
                    else:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + b"NNNNNNNNNN"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + b"NNNNNNNNNN"
                elif sampletype[primerset] == "RNA":
                    #print(seq2[124:132])
                    if N6type[primerset] == "F":
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    elif Levenshtein.distance(seq2[124:132].decode(),"TTTTTTTT") < 3:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124]
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124]
                    else:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + b"NNNNNNNNNN"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + b"NNNNNNNNNN"
                    if Levenshtein.distance(seq1[0:20].decode(),"AAGCAGTGGTATCAACGCAG") < 3:
                        seqhead1 = seqhead1 + b"_" + b"PCRRead" + b"\n"
                        seqhead2 = seqhead2 + b"_" + b"PCRRead" + b"\n"
                        seq1 = seq1[30:80]
                        qual1 = qual1[30:80]
                    else:
                        seqhead1 = seqhead1 + b"_" + b"TagRead" + b"\n"
                        seqhead2 = seqhead2 + b"_" + b"TagRead" + b"\n"
                        seq1 = seq1[0:50]
                        qual1 = qual1[0:50]
                        # print(seqhead1)
                elif sampletype[primerset] == "RNA2":
                    #print(seq2[124:132])
                    index_new = barcode.find("P1")
                    new_barcode = barcode[ :index_new] + printPrimerSet[primerset]
                    if N6type[primerset] == "F":
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + seq2[114:124] + b"\n"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + seq2[114:124] + b"\n"
                    elif Levenshtein.distance(seq2[124:132].decode(),"TTTTTTTT") < 3:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + seq2[114:124]
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + seq2[114:124]
                    else:
                        seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + b"NNNNNNNNNN"
                        seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(new_barcode) + b"_" + b"NNNNNNNNNN"
                    if Levenshtein.distance(seq1[0:20].decode(),"AAGCAGTGGTATCAACGCAG") < 3:
                        seqhead1 = seqhead1 + b"_" + b"PCRRead" + b"\n"
                        seqhead2 = seqhead2 + b"_" + b"PCRRead" + b"\n"
                        seq1 = seq1[30:80]
                        qual1 = qual1[30:80]
                    else:
                        seqhead1 = seqhead1 + b"_" + b"TagRead" + b"\n"
                        seqhead2 = seqhead2 + b"_" + b"TagRead" + b"\n"
                        seq1 = seq1[0:50]
                        qual1 = qual1[0:50]
                elif sampletype[primerset] == "crop":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    index3 = seq1.decode().find("GTTTTAG")
                    # remove 22 bp on 5' and anything after GTTTTAG
                    if index3 == -1:
                        seq1 = seq1[22: ]
                        qual1 = qual1[22: ]
                    else:
                        seq1 = seq1[22: index3]
                        qual1 = qual1[22: index3]
                elif sampletype[primerset] == "cite":
                    # remove 21 bp on 5' and keep 22-32
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    seq1 = seq1[21: 31]
                    qual1 = qual1[21: 31]
                elif sampletype[primerset] == "cellhash":
                    # remove 20 bp on 5' and keep 21-31
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"_" + seq2[114:124] + b"\n"
                    seq1 = seq1[20: 30]
                    qual1 = qual1[20: 30]
                elif sampletype[primerset] == "notrim":
                    seqhead1 = seqhead1[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index-1] + b"_" + str.encode(barcode) + b"\n"
                else:
                    seqhead1 = seqhead1[ :index] + b"1:N:0:" + str.encode(barcode) + b"\n"
                    seqhead2 = seqhead2[ :index] + b"2:N:0:" + str.encode(barcode) + b"\n"
        # align reads to themselves
        i = i+1  # total reads

        # trim reads
#        if needtrim == "T":
#            rc_seq2 = reverse_complement(seq2[0:n])
#            idx = seq1.rfind(str.encode(rc_seq2)) # look for perfect match
#        
#            if idx > 0:
#                j = j+1  # 0 mismatchs
#            elif mismatch>0:
#                hold = fuzz_align(rc_seq2,seq1.decode(),mismatch)  # else allow for mismatch
#                if hold:
#                    idx,mis=hold
#                    if mis == 1:
#                        k=k+1  # 1 mismatch
#
#            # trim reads if idx exist
#            if idx > 0:
#                # keep track on how much trimming
#                tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
#            
#                # trim data
#                seq1 = seq1[0:idx+n-1]
#                # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
#                seq2 = seq2[0:idx+n-1]
#                qual1 = qual1[0:idx+n-1]
#                qual2 = qual2[0:idx+n-1]
        # print data
        if barcodeMatch == 5:
            if primerset in project:                 
                outName = project[primerset]
                f1 = files_r1[outName]
                f2 = files_r2[outName]               
                f1.write(seqhead1);f1.write(seq1+b"\n")
                f1.write(qualhead1);f1.write(qual1+b"\n")
                f2.write(seqhead2);f2.write(seq2+b"\n")
                f2.write(qualhead2);f2.write(qual2+b"\n")
                good_r = good_r + 1
#        else:
#            dis1.write(seqhead1);dis1.write(seq1+b"\n")
#            dis1.write(qualhead1);dis1.write(qual1+b"\n")
#            dis2.write(seqhead2);dis2.write(seq2+b"\n")
#            dis2.write(qualhead2);dis2.write(qual2+b"\n")
        if options.qc == True:
            if i > qcreads:
                break
    # increment count
    count = count + 1
    if count == 5:
        count = 1
    else:
        count = count

# close files to write the file
for f in files_r1.values():
    f.close()
for f in files_r2.values():
    f.close()
p1_rds.close();p2_rds.close()
#dis1.close();dis2.close()
time = (time.process_time() - start)/60

# print("%.2g" % time + " min comsumed")

## give summary
try:
#    print(str(i)+" sequences total")
#    print(str(j)+" sequences trimmed with 0 mismatches")
#    print(str(k)+" sequences trimmed with 1 mismatch")
#    mean = tot_b/(j+k)
#    print("%.3g" % mean +" mean number of bases trimmed for reads requiring trimming")
    print(str(good_r)+" sucessfully demultiplexed for "+prefix)
    perc = good_r/(i)*100
    print("%.3g" % perc + "% reads demultiplexed")
except ZeroDivisionError: 
    print("Warning too few reads")
