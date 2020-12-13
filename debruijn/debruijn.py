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

from random import randint
import random
import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
random.seed(9001)
import statistics

__author__ = "Julien Bonin"
__copyright__ = "EISTI (CY Tech)"
__credits__ = ["@julienbnn"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Julien Bonin"
__email__ = "julien729@hotmail.fr"
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


def read_fastq(fastq_file):
    with open(fastq_file, "rt") as monfic:
        for line in monfic:
            yield next(monfic).rstrip("\n")
            next(monfic)
            next(monfic)


def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    '''
    read = ""

    for seq in read_fastq(fastq_file):
        read += seq
    '''
    read = next(read_fastq(fastq_file))
    kmer_list = []
    for k in cut_kmer(read, kmer_size):
        kmer_list.append(k)
    kmer_dict = dict((x, kmer_list.count(x)) for x in set(kmer_list))
    return kmer_dict

def build_graph(kmer_dict, display=False):
    H = nx.Graph()
    G = nx.DiGraph(H)
    for k in kmer_dict:
        G.add_edge(k[:-1], k[1:], weight=kmer_dict[k])
    edge_labels = dict([((u, v, ), d['weight']) for u, v, d in G.edges(data=True)])
    pos = nx.spring_layout(G)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    nx.draw(G, pos, with_labels=True, font_weight='bold')
    if display:
        plt.show()
    return G

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        for i in range(1, len(path)-1):
            graph.remove_node(path[i])

    if delete_entry_node:
        for path in path_list:
            graph.remove_node(path[0])

    if delete_sink_node:
        for path in path_list:
            graph.remove_node(path[-1])
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    max_weight = max(weight_avg_list)
    ind_list = [index for index, element in enumerate(weight_avg_list) if element == max_weight]
    if len(ind_list) > 1:
        new_path_list = [path_list[i] for i in ind_list]
        new_path_length = [path_length[i] for i in ind_list]
        max_length = max(new_path_length)
        ind_list = [index for index, element in enumerate(new_path_length) if element == max_length]
    else:
        path = path_list[ind_list[0]]
    if len(ind_list) > 1:
        new_path_list = [new_path_list[i] for i in ind_list]
        path = new_path_list[randint(0, len(new_path_list)-1)]
    else:
        path = path_list[ind_list[0]]

    path_list.pop(path_list.index(path))
    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)



def path_average_weight(graph, path):
    mean = 0
    for i in range(len(path)-1):
        mean += graph[path[i]][path[i+1]]["weight"]
    return mean/(len(path)-1)

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    if len(path_list) == 0:
        return graph
    path_length = []
    path_avg_weight = []
    for p in path_list:
        path_length.append(len(p))
        path_avg_weight.append(path_average_weight(graph, p))

    return select_best_path(graph, path_list, path_length, path_avg_weight)

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    # get the nodes with in_degree equal to 0
    starting_nodes = [u for u, deg in graph.in_degree() if not deg]
    return starting_nodes

def get_sink_nodes(graph):
    # get the nodes with out_degree equal to 0
    output_nodes = [u for u, deg in graph.out_degree() if not deg]
    return output_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contig_list = []
    for i in starting_nodes:
        for j in ending_nodes:
            path = nx.shortest_path(graph, source=i, target=j)
            way = path[0]
            for p in range(1, len(path)):
                way += path[p][-1]
            contig_list.append((way, len(way)))
    return contig_list

def save_contigs(contigs_list, output_file):
    with open(output_file, 'w') as f:
        id = 0
        for path, path_length in contigs_list:
            id += 1
            f.write(">contig_" + str(id) + " len=" + str(path_length) + "\n" + fill(str(path)) + "\n")
    f.close()

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    #graph = build_graph(kmer_dict, True)

    draw_graph(graph, 'graph.png')

    start_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)

    contig_list = get_contigs(graph, start_nodes, sink_nodes)
    save_contigs(contig_list, args.output_file)

    #graph2 = solve_bubble(graph, a, d)

if __name__ == '__main__':
    main()
