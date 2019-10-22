#!/usr/bin/env python3


import argparse
import sys
import pprint


parser = argparse.ArgumentParser(description="Debruij algo, main program.")

parser.add_argument("-i", "--fastq",
                    help='Path to the fastq file single end',
                    dest="fastq", metavar="fastq--file")
parser.add_argument("-k", "--kmer-size", help="size of the Kmer --DEFAULT = 21",
                    type=int, default=21, dest="sKmer", metavar="size--Kmer")
parser.add_argument("-r", "--reference-genome",
                    help='path to the file containing the reference genome ',
                     dest="geneRef", metavar="reference--genome")
parser.add_argument("-o", "--config-file",
                    help='path to the file containing the config file ',
                     dest="config", metavar="config--file")

args = parser.parse_args()

#Utiliser les librairies networkx, pytest et pyllint de Python:

#1° CREATION DU GRAPHE de DE BRUIJN
##  a) Identification des kmer unique

def read_fastq(fastq):
    fastq_file = open(fastq)
    for line in fastq_file:
        yield fastq_file.readline()
        seq_id = fastq_file.readline()
        seq = fastq_file.readline().strip('\n')

def cut_kmer(fastq, taille_kmer):
    it = read_fastq(fastq)
    seq = next(it)
    for j in range(len(seq) - taille_kmer + 1) :
        k_mer = seq[j:j+taille_kmer]
        #print(k_mer)
        yield(k_mer)

def dico_kmer(fastq , taille_kmer):
    it = cut_kmer(fastq , taille_kmer)
    dico = {}
    for kmer in it:
        if kmer not in dico :
            dico[kmer] = 0
        dico[kmer] += 1
    return dico

# MAIN :     
file = args.fastq
dic = dico_kmer(file , 3)
print(dic)
