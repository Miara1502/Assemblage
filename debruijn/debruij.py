#!/usr/bin/env python3

import argparse
import sys
import pprint
import networkx as nx


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
        yield fastq_file.readline().strip('\n')
        seq_id = fastq_file.readline()
        seq = fastq_file.readline().strip('\n')
    fastq_file.close()

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

##  b) Construction de l'arbre :
def build_graph(dico_kmer) :
    G = nx.DiGraph()
    for i , (kmer, nb) in enumerate(dico_kmer.items()):
        node1 = kmer[:-1]
        node2 = kmer[1:]
        #print(node1, node2, nb)
        G.add_edge(node1 , node2 , weight = nb)
    return G


def starting_nodes(graph) :
    list_entre = []
    for node in graph :
        pred = list(graph.predecessors(node))
        if (not pred) :
            #print("Pas de predecesseur\n")
            list_entre.append(node)
    return list_entre

def sink_nodes(graph) :
    list_sortie = []
    for node in graph :
        pred = list(graph.successors(node))
        if (not pred) :
            #print("Pas de successeurs\n")
            list_sortie.append(node)
    return list_sortie

def get_contigs(graph , noeud_entre , noeud_sortie) :
    '''Retourne (contig , taille du contig)
    '''
    for i in noeud_entre :
        liste_path = []
        for j in noeud_sortie :
            path = list(nx.all_simple_paths(graph , i , j))

            print('HELLO\n')
            tuple  = (path , len(path[0]))
            liste_path.append(tuple)
    return liste_path


################################# MAIN ##################################
file = args.fastq
nb = args.sKmer
dic = dico_kmer(file , nb)

g = build_graph(dic)


entree = starting_nodes(g)
sortie = sink_nodes(g)
print('Liste des noeuds d entrée : ')
print(entree)
print('\n')
print('Liste des noeuds de sortie : ')
print(sortie)


#Contigs :

contig = get_contigs(g , entree , sortie)
print(contig)
'''
path = list(nx.all_simple_paths(g , entree[0] , sortie[0]))

print(path)
print(len(path[0]))
'''
