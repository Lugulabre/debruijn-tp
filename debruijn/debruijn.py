#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq):
    '''Renvoie les séquences contenues dans un fichier fastq
    '''
    with open(fastq, "r") as filin:
        for line in filin:
            yield(next(filin).strip())
            next(filin)
            next(filin)


def cut_kmer(seq, k):
    '''Renvoie les k-mer de taille k d'une séquence
    '''
    for base in range(len(seq)-k+1):
        yield(seq[base:base+k])


def build_kmer_dict(fastq, k):
    '''Renvoie dictionnaire de k-mer
    '''
    dict_kmer = {}
    for seq in read_fastq(fastq):
        for kmer in cut_kmer(seq, k):
            if not kmer in dict_kmer:
                dict_kmer[kmer] = 1
            else:
                dict_kmer[kmer] += 1

    return(dict_kmer)


def build_graph(dict_kmer):
    '''Créer graphe à partir des k-mers
    '''
    G = nx.DiGraph()
    for key, value in dict_kmer.items():
        a = key[:-1]
        b = key[1:]
        G.add_node(a)
        G.add_node(b)
        G.add_edge(a, b, weight = value)
    #nx.draw(G, with_labels = True)
    #plt.show()
    return(G)


def get_starting_nodes(G):
    '''Trouver nœuds d'entrée
    '''
    list_start = []
    for node in G.nodes.data():
        if not list(g.predecessors(node[0])):
            list_start.append(node[0])

    return(list_start)


def get_sink_nodes(G):
    '''Trouver nœuds de sortie
    '''
    list_end = []
    for node in G.nodes.data():
        if not list(g.successors(node[0])):
            list_end.append(node[0])

    return(list_end)


def get_contigs(g, list_start, list_end):
    '''Retourner tuple de contigs
    '''
    list_contig = []
    for start in list_start:
        for end in list_end:
            for path in nx.all_simple_paths(g, start, end):
                list_contig.append((path, len(path)))

    return(list_contig)


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(list_tupple, filename):
    '''Sauvegarder contigs dans fichier de sortie
    '''
    with open(filename, "w") as filout:
        for i in range(len(list_tupple)):
            contig = "".join(list_tupple[i][0][0])
            for read in range(1, len(list_tupple[i][0])):
                contig = contig + list_tupple[i][0][read][-1]

            filout.write(">contig_{} len={}\n{}\n".format(i, len(contig), fill(contig)))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':

    g = build_graph(build_kmer_dict("../data/eva71_two_reads.fq", 20))
    start = get_starting_nodes(g)
    end = get_sink_nodes(g)
    contig = get_contigs(g, start, end)
    save_contigs(contig, "output.fasta")

    main()


