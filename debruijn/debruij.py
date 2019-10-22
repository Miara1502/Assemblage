#!/usr/bin/env python3


import argparse
import sys
import pprint


parser = argparse.ArgumentParser(description="Debruij algo, main program.")

parser.add_argument("-i", "--fastq", type=argparse.FileType('r'),
                    help='Path to the fastq file single end',
                    dest="fastq", metavar="fastq_file")
parser.add_argument("-k", "--kmer-size", help="size of the Kmer --DEFAULT = 21",
                    type=int, default=21, dest="sKmer", metavar="size_Kmer")
parser.add_argument("-r", "--reference-genome", type=argparse.FileType('r'),
                    help='path to the file containing the reference genome ',
                     dest="geneRef", metavar="reference_genome")
parser.add_argument("-o", "--config-file", type=argparse.FileType('r'),
                    help='path to the file containing the config file ',
                     dest="config", metavar="config_file")

args = parser.parse_args()

#Utiliser les librairies networkx, pytest et pyllint de Python:
