# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import numpy as np
import itertools
from scipy.stats import laplace, norm
import re
import math
import heapq

def get_regression_scores(filename):
    coeffs = np.loadtxt(filename)
    
    score_mat = coeffs.reshape((-1, 4)).T

    return score_mat

def get_scores(edges, kmers, score_mat):
    base4 = {'A':0, 'C':1, 'G':2, 'T':3}
    scores= np.zeros_like(edges)
    for node1 in range(edges.shape[0]):
        for node2 in range(edges.shape[1]):

            ### For regression value
            if edges[node1, node2] == 1:

                val1, val2 = 0,0
                for pos, (a,b) in enumerate(zip(kmers[node1], kmers[node2])):#[:-1]
                    
                    if a not in ['2', '3']:
                        val1 += score_mat[base4[a], pos]
                    if b not in ['2', '3']:
                        val2 += score_mat[base4[b], pos]
                scores[node1, node2] = np.abs(val1-val2)
            
    return scores

def get_scores_m(edges, kmers, score_file=''):
    means = np.loadtxt(score_file, skiprows=1, usecols=[0,1], dtype=str)

    score_vals = {k: v for k, v in zip(means[:, 0], means[:, 1].astype(float))}

    scores = np.zeros_like(edges)
    for node1 in range(edges.shape[0]):
        for node2 in range(edges.shape[1]):

            if edges[node1, node2] == 1:
                k1, k2 = kmers[node1], kmers[node2]
                val1 = score_vals[k1]
                val2 = score_vals[k2]
                scores[node1, node2] = np.abs(val1-val2)
    return scores, score_vals

def get_all_kmers(k, alphabet):
    l = [''.join(i) for i in itertools.product(alphabet, repeat = k)]
    return np.array(l)


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


def remove_nodes(edges, scores, numerator, v=0):
    loop = True  
    if v==0 or v==2:
        last_sums = np.zeros_like((np.sum(edges, axis=1)))
        while loop==True:
            
            edge_sums = np.sum(edges, axis = 1)
            if np.array_equal(edge_sums, last_sums): loop=False
            
            idxs_bad = [edge_sums<1][0]
            edges[:, idxs_bad] = 0
            
            last_sums=edge_sums

    if v==1 or v==2:
        for row in range(edges.shape[0]):
            pmoves = np.nonzero(edges[row, :])[0]
            sorted_idx = np.argsort(scores[row, pmoves])
            pmoves = pmoves[sorted_idx]
            pmoves = pmoves[2**numerator:]

            edges[row, pmoves] = 0

    return edges


class model:
    def __init__(self, k, delta, b, score_file, trans, max_run=0, score='regression'):

        self.trans = trans
        # Estimated max standard deviation
        b = np.array(b)
        self.b_max = b
        self.k = k
        
        add = '6'
        if delta == 0: add = ''
        rate_dict = {0:2, 1:10, 2:10, 3:10, 4:10, 5:8, 6:8, 7:8, 8:8, 9:8, 10:6}

        set1 = ['A', 'C', 'G', 'T']
        # All possible states
        kmers = get_all_kmers(k, set1)
        ks1 = get_all_kmers(k, set1)
        ks = ks1
        
        self.States = np.array(ks)
        
        edges, edges_set = get_debruijn_edges_from_kmers(kmers, max_run=max_run)
        

        scores, scores_vals = get_scores_m(edges, kmers)
        
        mask1 = [scores<delta][0]
        edges[mask1] = 0

        rate_dict = {0:2, 1:10, 2:10, 3:10, 4:10, 5:8, 6:8, 7:8, 8:8, 9:8, 10:6}
        if delta>0:
            edges = remove_nodes(edges, scores, rate_dict[delta], v=0)
        for d in range(edges.shape[0]):
            edges[d,d] = 1
        
        # Transistion probabilities based on de Bruijn constraints
        new_matrix = np.zeros_like(edges)
        edge_counts = np.sum(edges, axis=1)
        
        # print('Transition: ', self.trans)
        for r in range(edges.shape[0]): 

            new_matrix[r, :] = edges[r, :]*((1/self.trans)/edge_counts[r])
            new_matrix[r,r] = 1-(1/self.trans)
            
        new_matrix[new_matrix==np.nan] = 0

        self.M = new_matrix

        emiss = [*scores_vals.values()]
        # emiss = np.tile(emiss, self.dwell)
        emiss = np.array(emiss)
        self.E = emiss


        # Starting distribution for all states -> Uniform
        check = edge_counts
        check[check>=1]=1 
        self.Start = check/np.sum(check)


    def states(self):
        return np.array(self.States)

    def start(self, state):
        return np.log(self.Start[state])
#        return self.Start[state]

    def emit(self, state, symbol):

        mu = self.E[state]
        sig = self.b_max[state]
        return np.log(norm.pdf(symbol, loc=mu, scale=sig)+1e-50) 
    
    def transition(self, state1, state2, qs=None, states=None):
        return np.log(self.M[state1, state2])

    def emit_vec(self, symbols, states=None):
        if states is None: 
            mus = self.E
            sig = self.b_max
        else: 
            mus = self.E[states]
            sig = self.b_max[states]
        return np.log(norm.pdf(symbols, loc=mus, scale=sig)+1e-50) 
    
    
def viterbi(observations, model):
    states = model.states()
    
    d = -np.inf*np.ones((states.shape[0], observations.shape[0]))

    # Step 1: Initialization
    print('Initializing Model...')
    for i, state in enumerate(states):
        d[i, 0] = model.start(i) + model.emit(i, observations[0])

    # Step 2: Recursion
    print('Recursion Step...')
    for t in range(1, observations.shape[0]):
        
        emit = model.emit_vec(observations[t]).reshape(1, -1) 
        
        ds = d[:, t-1].reshape(-1, 1)

        val_vec =  ds + np.log(model.M) + emit
        val_vec = np.max(val_vec, axis=0)
        
        d[:, t] = np.maximum(val_vec, d[:, t])

    # Step 3: Traceback
    print('Traceback Step...')
    qs = []
    qs_idxs = []
    for t in reversed(range(observations.shape[0])):
        best_state = None
        best_state_val = -np.inf
        for i, state in enumerate(states):
            transition_prob = 0
            if t < observations.shape[0] - 1: 
                qs_idx = qs_idxs[0]
                transition_prob = model.transition(i, qs_idx, qs=qs, states=states)

            if d[i, t]+ transition_prob >= best_state_val:
                best_state_val = d[i, t] + transition_prob
                best_state = state
                best_state_idx = i
            
        qs.insert(0, best_state)
        qs_idxs.insert(0, best_state_idx)

    return qs, qs_idxs
        