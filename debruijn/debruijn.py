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
import argparse
from operator import itemgetter
import os
import random
from random import randint
import statistics
import sys
import networkx as nx
"""Perform assembly based on debruijn graph."""

random.seed(9001)

__author__ = "Pretet Mael"
__copyright__ = "Universite de Paris"
__credits__ = ["Pretet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Pretet"
__email__ = "pretetmael@gmail.com"
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
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
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
            yield next(filin).strip()
            next(filin)
            next(filin)


def cut_kmer(seq, k):
    '''Renvoie les k-mer de taille k d'une séquence
    '''
    for base in range(len(seq)-k+1):
        yield seq[base:base+k]


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

    return dict_kmer


def build_graph(dict_kmer):
    '''Créer graphe à partir des k-mers
    '''
    graph = nx.DiGraph()
    for key, value in dict_kmer.items():
        key_pre = key[:-1]
        key_suf = key[1:]
        graph.add_node(key_pre)
        graph.add_node(key_suf)
        graph.add_edge(key_pre, key_suf, weight=value)

    return graph


def get_starting_nodes(graph):
    '''Trouver noeuds d'entrée
    '''
    list_start = []
    for node in graph.nodes.data():
        if not list(graph.predecessors(node[0])):
            list_start.append(node[0])

    return list_start


def get_sink_nodes(graph):
    '''Trouver noeuds de sortie
    '''
    list_end = []
    for node in graph.nodes.data():
        if not list(graph.successors(node[0])):
            list_end.append(node[0])

    return list_end


def get_contigs(graph, list_start, list_end):
    '''Retourner tuple de contigs
    '''
    list_contig = []
    for start in list_start:
        for end in list_end:
            for path in nx.all_simple_paths(graph, start, end):
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]

                list_contig.append((contig, len(contig)))

    return list_contig


def fill(text, width=80):
    '''Split text with a line return to respect fasta format
    '''
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(list_tupple, filename):
    '''Sauvegarder contigs dans fichier de sortie
    '''
    with open(filename, "w") as filout:
        for i in range(len(list_tupple)):
            contig = list_tupple[i][0]
            l_contig = list_tupple[i][1]
            filout.write(">contig_{} len={}\n{}\n".format(
                i, l_contig, fill(contig)))


def std(list_val):
    '''Calcul écart-type d'une liste de valeurs
    '''
    if len(list_val) == 1:
        return 0
    return statistics.stdev(list_val)


def path_average_weight(graph, path):
    weight = 0
    for i in range(len(path)-1):
        weight += graph.edges[path[i], path[i+1]]["weight"]

    return weight/(len(path)-1)


def remove_paths(graph, list_path, delete_entry_node, delete_sink_node):
    '''Supprimer un chemin donné du graphique
    '''
    start = 1
    end = 1
    if delete_entry_node:
        start = 0
    if delete_sink_node:
        end = 0

    for path in list_path:
        for node in range(start, len(path)-end):
            graph.remove_node(path[node])

    return graph


def select_best_path(graph, list_path, list_len, list_weight, 
    delete_entry_node=False, delete_sink_node=False):
    '''Sélectionner le meilleur chemin d'une liste de chemin
    '''
    max_weight = max(list_weight)
    max_length = max(list_len)
    pos_max = list(range(0, len(list_path)))
    del_path = []
    pos_tmp = list(pos_max)

    if std(list_weight) != 0:
        for i in pos_tmp:
            if list_weight[i] < max_weight:
                del_path.append(list_path[i])
                pos_max.remove(i)

    list_len = [list_len[i] for i in pos_max]
    pos_tmp = list(pos_max)
    if len(pos_max) > 1 and std(list_len) != 0:
        for i in pos_tmp:
            if list_len[i] < max_length:
                del_path.append(list_path[i])
                pos_max.remove(i)

    if len(pos_max) > 1:
        keep_branch = randint(0,len(pos_max))
        pos_max.remove(keep_branch)
        for i in pos_max:
            del_path.append(list_path[i])

    return remove_paths(graph, del_path, delete_entry_node, delete_sink_node)


def solve_bubble(graph, ancetre, descend):
    '''Résoudre une bulle donnée en entrée en conservant le meilleur chemin
    '''
    list_path = []
    list_weight = []
    list_len = []
    for path in nx.all_simple_paths(graph, ancetre, descend):
        list_path.append(path)
        list_weight.append(path_average_weight(graph, path))
        list_len.append(len(path))

    return select_best_path(graph, list_path, list_len, list_weight)


def simplify_bubbles(graph):
    '''Trouver les bulles présentes dans un graphe et les résoudre
    '''
    list_bubble = []
    for node in graph.nodes:
        list_predecessors = list(graph.predecessors(node))
        if len(list_predecessors) > 1:
            lwca = nx.lowest_common_ancestor(graph, list_predecessors[0], list_predecessors[1], default = -1)
            list_bubble.append([lwca, node])

    for i in range(0, len(list_bubble)):
        graph = solve_bubble(graph, list_bubble[i][0], list_bubble[i][1])

    return graph


def solve_entry_tips(graph, list_in):
    '''Nettoie les entrées indésirables d'un graphe
    '''
    if len(list_in) > 1:
        for node_in in list_in:
            list_descendants = nx.descendants(graph, node_in)
            for descendant in list_descendants:
                list_predecessors = list(graph.predecessors(descendant))
                if len(list_predecessors) > 1:
                    list_path = []
                    list_weight = []
                    list_len = []
                    entry_ancestors = []
                    for entry in nx.ancestors(graph, descendant):
                        if entry in list_in:
                            entry_ancestors.append(entry)

                    for entry in entry_ancestors:
                        path_to_test = nx.shortest_path(graph, entry, descendant)
                        list_path.append(path_to_test)
                        list_weight.append(path_average_weight(graph, path_to_test))
                        list_len.append(len(path_to_test))

                    graph = select_best_path(graph, list_path, list_len, list_weight, True, False)
                    list_in_actual = []
                    for entry in list_in:
                        if entry in list(graph.nodes):
                            list_in_actual.append(entry)
                    return solve_entry_tips(graph, list_in_actual)

    return graph


def solve_out_tips(graph, list_out):
    '''Nettoie les sorties indésirables d'un graphe
    '''
    if len(list_out) > 1:
        for node_out in list_out:
            list_ancestors = nx.ancestors(graph, node_out)
            for ancestor in list_ancestors:
                list_successors = list(graph.successors(ancestor))
                if len(list_successors) > 1:
                    list_path = []
                    list_weight = []
                    list_len = []
                    exit_descendants = []
                    for exit in nx.descendants(graph, ancestor):
                        if exit in list_out:
                            exit_descendants.append(exit)

                    for exit in exit_descendants:
                        path_to_test = nx.shortest_path(graph, ancestor, exit)
                        list_path.append(path_to_test)
                        list_weight.append(path_average_weight(graph, path_to_test))
                        list_len.append(len(path_to_test))

                    graph = select_best_path(graph, list_path, list_len, list_weight, False, True)
                    list_out_actual = []
                    for exit in list_out:
                        if exit in list(graph.nodes):
                            list_out_actual.append(exit)
                    return solve_entry_tips(graph, list_out_actual)

    return graph


# ==============================================================
# Main program
# ==============================================================


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    G = build_graph(kmer_dict)
    ls = get_starting_nodes(G)
    le = get_sink_nodes(G)
    G = simplify_bubbles(G)
    G = solve_entry_tips(G, ls)
    G = solve_out_tips(G, le)
    lc = get_contigs(G, ls, le)
    save_contigs(lc, args.output_file)


if __name__ == '__main__':
    main()

