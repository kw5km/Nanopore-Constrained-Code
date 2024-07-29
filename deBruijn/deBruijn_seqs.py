# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""

import deBruijn_helper
import numpy as np
from scipy import sparse
from fractions import Fraction
import math
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from hamming import calcRedundantBits, posRedundantBits, calcParityBits


def get_regression_scores(filename):
    coeffs = np.loadtxt(filename)
    
    score_mat = coeffs.reshape((-1, 4)).T

    return score_mat
    
def remove_nodes(edges):
    loop = True
    
    last_sums = np.zeros_like((np.sum(edges, axis=1)))
    while loop==True:
        
        edge_sums = np.sum(edges, axis = 1)
        if np.array_equal(edge_sums, last_sums): loop=False
        
        idxs_bad = [edge_sums<1][0]
        edges[:, idxs_bad] = 0
        
        last_sums=edge_sums
    
    return edges

def get_fpad(delta, edges, kmers, scores, rate, fpad):
    if fpad is None: return kmers[0], ''
    else:
        seq = ''
        idx_nonzero = np.nonzero(edges)
        idx=0
        next_bases = kmers[idx_nonzero[0][idx]]
        seq = seq+str(next_bases)
        for i in range(0, fpad.shape[0], rate.numerator):             
            base_idx = np.where(kmers==next_bases)[0][0]
            possible_moves = np.nonzero(edges[base_idx, :])[0]

            sorted_idx = np.argsort(scores[base_idx, possible_moves])
                
            possible_moves = possible_moves[sorted_idx]

            moves_idx = fpad[i:i+rate.numerator]
            moves_idx = int("".join(str(x) for x in moves_idx), 2)

            next_idx = possible_moves[moves_idx]
            next_bases = kmers[next_idx]

            seq = seq+str(next_bases[-rate.denominator:])
            
        return next_idx, next_bases, seq
            
def deBruijn_sequence(edges, delta, kmers, scores, message=None, rate=Fraction(2,1), padding=18, start_idx=0, fpad=None):
    if message is None:
        mask = [scores<delta][0]
        edges[mask] = 0
        edges = remove_nodes(edges)

    n_idx, next_bases, seq = get_fpad(delta, edges, kmers, scores, rate, fpad)

    for i in range(0, message.shape[0], rate.numerator):
        base_idx = n_idx 
        possible_moves = np.nonzero(edges[base_idx, :])[0]

        if message is not None:
            sorted_idx = np.argsort(scores[base_idx, possible_moves])
                
            possible_moves = possible_moves[sorted_idx]

            moves_idx = message[i:i+rate.numerator]
            moves_idx = int("".join(str(x) for x in moves_idx), 2)
        else:
            moves_idx = np.random.randint(possible_moves.size, size=1)[0]

        # print(moves_idx, possible_moves)
        next_idx = possible_moves[moves_idx]
        next_bases = kmers[next_idx]
        n_idx = next_idx

        seq = seq+str(next_bases[-rate.denominator:])

    pad = ''
    for i in range(padding//rate.denominator):
        base_idx = n_idx 
        possible_moves = np.nonzero(edges[base_idx, :])[0]
        sorted_idx = np.argsort(scores[base_idx, possible_moves])
        possible_moves = possible_moves[sorted_idx][:2**(rate.numerator)]
        moves_idx = np.random.randint(possible_moves.size, size=1)[0]
        next_idx = possible_moves[moves_idx]
        next_bases = kmers[next_idx]
        n_idx = next_idx

        seq = seq+str(next_bases[-rate.denominator:])
        pad = pad+str(next_bases[-rate.denominator:])
    print('shape bpad:', len(pad))

    return seq
       
def divisorGenerator(n):
    large_divisors = []
    for i in range(1, int(math.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor

def produce_save_seqs(n, bp, edges, delta, kmers, filename, score_dict, message='rand', rate=Fraction(2,1), message_file=None):
    print('Job Delta {} Started'.format(delta), message)
    if message == 'rand':
            print('Encoding Random')
            messages = np.random.choice([0, 1], size=(n, math.ceil((bp*rate.numerator)/rate.denominator)), p=[.5, .5])
    elif message is not None:
        print('Encoding File')
        
        messages = message

        # subtracting 8 bits here for ID
        m_len = math.ceil((bp*rate.numerator)/rate.denominator) - 8
        if messages.size%(m_len)!=0:
            size = messages.size
            remainder = m_len - (size%(m_len))
            pad_char = not message.flat[-1]
            messages = np.pad(messages.astype(float), (0, remainder), mode='constant', constant_values=(pad_char,))
        messages = messages.reshape((-1, m_len))
        print('messages shape:', messages.shape)
        messages = messages.astype(int)

    #front padding, only 6 because initialization kmer adds 6, so result is 12
    pad_message = np.random.choice([0, 1], size=(math.ceil((6*rate.numerator)/rate.denominator),), p=[.5, .5])
    print('shape bits fpad:', pad_message.shape)
    
    seqs = []
    for i in range(messages.shape[0]):    
        if message is not None:
            rn = 8
            bs = f'{i:0{rn}b}'
            seq_idx = np.array(list(bs))       
            message2 = np.concatenate((seq_idx, messages[i, :])) 
            print('shape message2:', message2.shape)
        else: message2 = message[i, :]

        seq = deBruijn_sequence(edges, delta, kmers, score_dict, message=message2, rate=rate, fpad=pad_message)
        print('shape seq:', len(seq))
        if i%1000==0: 
            print('Job Delta {}: {} seqs'.format(delta, i))
            print('Seq Length:', len(seq))
        seqs.append(SeqRecord(Seq(seq), str(i), name='', description=''))

    SeqIO.write(seqs, filename, "fasta")

    np.savetxt(message_file, np.array(messages).astype(int), fmt='%i', delimiter=",")
    print('Job Delta {} Ended'.format(delta))
        
def mat_to_set(edges_mat, scores_mat, k, delta=0):

    good_edges=edges_mat
    edges_set = set()
    scores_set = []
    
    set1 = ['A', 'C', 'G', 'T']
    b = deBruijn_helper.get_all_kmers(k, set1)
    
    for row in range(good_edges.shape[0]):
        for col in range(good_edges.shape[1]):
            if good_edges[row, col] == 1:
                if (b[row][:-1], b[col][:-1]) not in edges_set:
                    scores_set.append(scores_mat[row, col])
                edges_set.add((b[row], b[col]))
               
    return edges_set, scores_set

       
def to_adj_mat(edges, delta, score_mat, k, kmers):
    edges = np.array(list(edges))
    
    scores = deBruijn_helper.get_scores(edges, kmers, score_mat)
    scores = np.array(scores)
    mask = [scores>=delta][0]

    
    good_edges = edges[mask]
    
    

    r, c = 4**k, 4**k

    adj_mat = np.zeros((r, c))
    
    set1 = ['A', 'C', 'G', 'T']
    b = list(deBruijn_helper.get_all_kmers(k, set1))
    
    for row in range(r):
        for col in range(c):

            edge = [b[row], b[col]]
            # print(edge)
            if edge in np.ndarray.tolist(good_edges):
                adj_mat[row, col] = 1
    return adj_mat

def get_rate(edges, delta=None, score_mat=None):

    if ((edges is not None) and (score_mat is not None)):
        mask = [score_mat<delta]
        edges[mask] = 0

    w, v = np.linalg.eig(edges)

    lam = np.max(w.real)

    return np.log2(lam)  
