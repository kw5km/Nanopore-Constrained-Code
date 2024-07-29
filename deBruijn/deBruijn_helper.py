# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import itertools
import toyplot
import toyplot.pdf
import numpy as np
import re

def get_all_kmers(k, alphabet):
    l = [''.join(i) for i in itertools.product(alphabet, repeat = k)]
    return np.array(l)

def get_kmer_count_from_sequence(sequence, k):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    kmers = {}
    
    # count how many times each occurred in this sequence (treated as cyclic)
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        #skip end incomplete kmers
        if len(kmer) != k:
            continue
        
        # count occurrence of this kmer in sequence
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    
    return kmers

def get_scores_means(edges, kmers, score_file='/u/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt'):
    means = np.loadtxt(score_file, skiprows=1, usecols=[0,1], dtype=str)

    score_vals = {k: v for k, v in zip(means[:, 0], means[:, 1].astype(float))}

    scores = np.zeros_like(edges)
    for node1 in range(edges.shape[0]):
        for node2 in range(edges.shape[1]):

            # if edges[node1, node2] == 1:
            k1, k2 = kmers[node1], kmers[node2]
            val1 = score_vals[k1]
            val2 = score_vals[k2]
            scores[node1, node2] = np.abs(val1-val2)
    return scores
    
def get_scores(edges, kmers, score_vals):
    base4 = {'A':0, 'C':1, 'G':2, 'T':3}
    scores= np.zeros_like(edges)
    for node1 in range(edges.shape[0]):
        for node2 in range(edges.shape[1]):

            ### For regression value
            # if edges[node1, node2] == 1:

            val1, val2 = 0,0
            for pos, (a,b) in enumerate(zip(kmers[node1][:-1], kmers[node2][:-1])):
                val1 += score_vals[base4[a], pos]
                val2 += score_vals[base4[b], pos]
            scores[node1, node2] = np.abs(val1-val2)
            
            ### For mean value
#            if edges[node1, node2] == 1:
##                print('nodes:', kmers[node1][:-1], kmers[node2][:-1])
#                k1, k2 = kmers[node1][:-1], kmers[node2][:-1]
#                val1 = score_vals[k1]
#                val2 = score_vals[k2]
#                scores[node1, node2] = np.abs(val1-val2)
    return scores

def get_colors(scores, delta):
    colors=[]
    for score in scores:
        if score>=delta:
            colors.append(toyplot.marker.create(shape="r3x1", size=10, label="%.2f" % score, mstyle={"fill":"green"}))
        else: 
            colors.append(toyplot.marker.create(shape="r3x1", size=10, label="%.2f" % score, mstyle={"fill":"white"}))
    return colors
            
        


def get_debruijn_edges_from_kmers(kmers, max_run=0):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps 
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = np.zeros((len(kmers), len(kmers)))
    edges_set = set()
    
    # compare each (k-1)mer
    for i, k1 in enumerate(kmers):
        for j, k2 in enumerate(kmers):
            # if k1 != k2: 
                # if they overlap then add to edges
            if k1[1:] == k2[:-1]:
                ks = k1+k2[-1]
                regex = fr'(\w)\1{{{max_run}}}'
                if (max_run == 0) or ((len(re.findall(regex, ks)) == 0)):
                    edges[i, j] = 1
                    edges_set.add((k1, k2))

    return edges, edges_set


def plot_debruijn_graph(edges_set, scores_set, color=None, width=1000, height=1000):

    layout = toyplot.layout.FruchtermanReingold()
#    mmarkers = [toyplot.marker.create(shape="r3x1", size=10, label="%.2f" % i, mstyle={"fill":"white"}) for i in scores]
    estyle={"stroke": "black", "stroke-width": 1}
    if color==None:
        mmarkers = [toyplot.marker.create(shape="r3x1", size=10, label="%.2f" % i, mstyle={"fill":"white"}) for i in scores_set]
        
    else: mmarkers = color
   
    "returns a toyplot graph from an input of edges"    
    canvas, coord, mark = toyplot.graph(
        [i[0] for i in edges_set],
        [i[1] for i in edges_set],
        mmarker=mmarkers,
        width=width,
        height=height,
        tmarker=">", 
        vsize=50,
        vstyle={"stroke": "black", "stroke-width": 2},
        vlstyle={"font-size": "11px"},
        estyle=estyle,
        layout=layout)
    return canvas, coord, mark
