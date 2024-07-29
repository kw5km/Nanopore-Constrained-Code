# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import sys
sys.path.append('../DeBruijn/')
import deBruijn_helper
import deBruijn_seqs
import state_splitting
import numpy as np
import os
from fractions import Fraction
import argparse


def remove_nodes(edges, scores, numerator=2, v=0):
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

def apply_split(f_o, f_m, f_s, n, bp, d, scores, edges, messages):

       delta, delta_str = d[0], d[1]
       print('Delta {} State Split Starting...:'.format(delta))

       mask = [scores<delta][0]
       edges[mask] = 0
       A = edges

       set1 = ['A', 'C', 'G', 'T']
       kmers = deBruijn_helper.get_all_kmers(6, set1) 
       
       if delta == 0:   
              r = 2
              deBruijn_seqs.produce_save_seqs(n, bp, A, delta, kmers, f_o+'delta{}_rate2.fasta'.format(delta_str), scores, message=messages, rate=Fraction(2, 1, _normalize=False),
                                                 message_file=f_m+'delta{}_rate2_message.txt'.format(delta))
              return

       denom = 6
       cap = deBruijn_seqs.get_rate(A)
       cap = int(cap * (10 ** 2)) / 10 ** 2
       rate = Fraction(int(np.floor(denom*cap)),denom, _normalize=False)

       num, denom = rate.numerator, rate.denominator
       if (num % 2) != 0: num-=1
       rate = Fraction(num, denom, _normalize=False)
       print('Delta {} Rate:'.format(delta), rate.numerator, rate.denominator, rate.numerator/rate.denominator)
       r = ''+str(rate.numerator)+str(rate.denominator)

       A_q, states_split = state_splitting.state_split(A, num, rate.denominator, kmers)
       scores2 = deBruijn_helper.get_scores_means(A_q, states_split, score_file=f_s) 

       # Produce Seqs
       deBruijn_seqs.produce_save_seqs(n, bp, A_q, delta, states_split, f_o+'delta{}_rate{}{}.fasta'.format(delta_str, num, denom), scores2, message=messages, rate=rate,
                                                 message_file=f_m+'delta{}_rate{}{}_message.txt'.format(delta, num, denom))

if __name__=='__main__':
       
       parser = argparse.ArgumentParser()
       parser.add_argument("-d", type=int, default=0, help="Delta val.")
       parser.add_argument("-b", type=int, default=200, help="Read rength")
       parser.add_argument("-n", type=int, default=1000,help="Num. reads")
       parser.add_argument("-k", type=int, default=6, help="k-mer size for de Bruijn Graph")
       parser.add_argument("-f_o", help="Output file for seqs")
       parser.add_argument("-f_m", help="Output file for messages")
       parser.add_argument("-f_s", help="k-mer mean file location")
       parser.add_argument("-max_run", type=int, default=0, help="Max length of homoplymer repeat allowed")
       parser.add_argument("-m", default='rand', help="Message to encode, 'rand' for random and filepath for file")
       args = parser.parse_args() 
       
       np.random.seed(seed=0)
       delta, bp, n, k, f_o, f_m, f_s, max_run, messages = args.d, args.b, args.n, args.k, args.f_o, args.f_m, args.f_s, args.max_run, args.m
       print('score file:', f_s)
       if messages != 'rand': messages = np.load(messages)
       print('Program Running...')
       set1 = ['A', 'C', 'G', 'T']
       kmers = deBruijn_helper.get_all_kmers(k, set1)

       edges, edges_set = deBruijn_helper.get_debruijn_edges_from_kmers(kmers, max_run=max_run)
       scores = deBruijn_helper.get_scores_means(edges, kmers, score_file=f_s)
       
       deltas = [(0, '000'), (1, '100'), (2, '200'), (3, '300'), (4, '400'), (5, '500'), (6, '600'), (7, '700'),(8, '800'), (9, '900'), (10, '1000')]# 
       d = deltas[delta]

       if d[0] == 0:
              edges0, edges_set0 = deBruijn_helper.get_debruijn_edges_from_kmers(kmers, max_run=0)
              scores0 = deBruijn_helper.get_scores_means(edges0, kmers, score_file=f_s)
              apply_split(f_o, f_m, f_s, n, bp, d, scores0, edges0, messages)
       else: 
              edges = remove_nodes(edges, scores, )
              apply_split(f_o, f_m, f_s, n, bp, d, scores, edges, messages)

