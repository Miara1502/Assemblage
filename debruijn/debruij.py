#!/usr/bin/env python3

import argparse
import sys , os
import pprint
import networkx as nx
from networkx import algorithms


#Utiliser les librairies networkx, pytest et pyllint de Python:

#1° CREATION DU GRAPHE de DE BRUIJN
##  a) Identification des kmer unique

def read_fastq(fastq):
    '''take a fastq file and return an iterator
    of a given sequence
    '''
    fastq_file = open(fastq)
    lines = iter(fastq_file.readlines())
    for line in lines :
        yield next(lines)
        next(lines)
        next(lines)


def cut_kmer(seq, taille_kmer):
    '''take a sequence and a size of kmer and return
    a kmer as an iterator
    '''
    seq = seq.strip('\n')
    for j in range(len(seq) - taille_kmer + 1):
        yield seq[j:j+taille_kmer]


def build_kmer_dic(fastq, taille_kmer):
    '''take a fastq file and a size of kmer and return
    a dictionnary with the kmer as keys() and his number
    of occurence as values()
    '''
    dico = {}
    for seq in read_fastq(fastq):
        for kmer in cut_kmer(seq, taille_kmer):
            print(kmer)
            if kmer not in dico:
                dico[kmer] = 0
            dico[kmer] += 1
    return dico


##  b) Construction de l'arbre :
def build_graph(dico_kmer):
    '''take a dictionnary of kmer and return the the
    suffix / prefix tree
    '''
    G = nx.DiGraph()
    for i , (kmer, nb) in enumerate(dico_kmer.items()):
        node1 = kmer[:-1]
        node2 = kmer[1:]
        G.add_edge(node1 , node2 , weight = nb)
    return G


def get_starting_nodes(graph):
    '''take a graph and return a list of entry
    nodes
    '''
    list_entre = []
    for node in graph :
        pred = list(graph.predecessors(node))
        if (not pred) :
            #print("Pas de predecesseur\n")
            list_entre.append(node)
    return list_entre

def get_sink_nodes(graph):
    '''take a graph and return a list of exit
    nodes
    '''
    list_sortie = []
    for node in graph :
        pred = list(graph.successors(node))
        if (not pred) :
            #print("Pas de successeurs\n")
            list_sortie.append(node)
    return list_sortie



def get_contigs(graph, list_start_node, list_end_node):
    '''takes a graph, a list of entry nodes and a list of exit nodes
    and returns a list of tuple (contig, contig_size)
    '''
    contigs = []
    for source in list_start_node :
        for target in list_end_node :
            if algorithms.has_path(graph, source, target) == True :
                path = algorithms.shortest_path(graph, source, target)
                contig = path[0]
                for i in range(len(path)-1):
                    contig += path[i+1][-1]
                contigs.append((contig, len(contig)))
        return contigs

def fill(text, width=80):
    return os.linesep.join(text[i:i+width] for i in range(0,
            len(text), width))


def save_contigs(Tuple , output_file):
    '''take a tuple and a name of ouotput file
    and give a file with the list of tuple
    '''
    with open (output_file , 'w+') as filout :
        #>contig_Numéro len=longueur du contig
        for i in range(len(Tuple)) :
            filout.write('>contig_' + str(i) + ' len=' + str(Tuple[i][1]) + '\n'+ str(fill(Tuple[i][0])) + '\n')
    filout.close()



parser = argparse.ArgumentParser(description="Debruij algo, main program.")

parser.add_argument("-i", "--fastq", type = str ,
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


################################# MAIN ##################################
if __name__ == '__main__' :
    file = args.fastq
    size = args.sKmer
    dic = build_kmer_dic(file , size)
    print(dic)
    print('\n')
    g = build_graph(dic)
    entree = get_starting_nodes(g)
    sortie = get_sink_nodes(g)
    print('Liste des noeuds d entrée : ')
    print(entree)
    print('\n')
    print('Liste des noeuds de sortie : ')
    print(sortie)
    print('\n')
    #Contigs :
    print('TUPLE CONTIG : ')
    contig = get_contigs(g , entree , sortie)
    print(contig)
    save_contigs(contig , "lol.txt")
