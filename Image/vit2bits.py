# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import numpy as np
import os, sys, re
from fractions import Fraction
import squiggle_helper
sys.path.append('../DeBruijn/')
import deBruijn_helper
import deBruijn_seqs
sys.path.append('../State_Splitting/')
import state_splitting
import argparse
from best_match import get_best_match

def bitstring_to_bytes(s):
    bits = np.array([int(b) for b in s])
    bytes = np.packbits(bits)
    return bytes
    # byte_hold = int(s, 2).to_bytes((len(s) + 7) // 8, byteorder='big')
    # return byte_hold#bytearray(byte_hold)

def get_graph(k, delta, max_run=2, f_s='/u/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt'):
    delt = str(int(delta))

    print('Generating State Split Graph')
    set1 = ['A', 'C', 'G', 'T']
    kmers = deBruijn_helper.get_all_kmers(k, set1)

    edges, edges_set = deBruijn_helper.get_debruijn_edges_from_kmers(kmers, max_run=max_run)
    scores = deBruijn_helper.get_scores_means(edges, kmers, score_file=f_s)

    mask = [scores<delta][0]
    edges[mask] = 0
    A = edges

    if delta == 0:
        return A, kmers, scores
    
    delt += '00'
    denom = k
    cap = deBruijn_seqs.get_rate(A)
    print('Capacity:', cap)
    cap = int(cap * (10 ** 2)) / 10 ** 2
    rate = Fraction(int(np.floor(denom*cap)),denom, _normalize=False)

    num, denom = rate.numerator, rate.denominator
    if (num % 2) != 0: num-=1
    rate = Fraction(num, denom, _normalize=False)
    print('Rate Calculated:', rate.numerator, rate.denominator, rate)

    A_q, states_split = state_splitting.state_split(A, num, rate.denominator, kmers)
    scores2 = deBruijn_helper.get_scores_means(A_q, states_split, score_file=f_s)

    print('Done')
    return A_q, states_split, scores2

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

def seq2bits(seq, edges, states, scores, rate, k=6, sq=0, idxs=None): 
    seq_len = int(186*(rate.numerator/rate.denominator))

    pos_bits = []
    kmers = squiggle_helper.rolling_window(np.array(list(seq)), k, overlap=rate.denominator)
    print('kmer len:', len(kmers))
    for i, kmer in enumerate(kmers[:-1]):
        next_kmer = kmers[i+1]

        if next_kmer == 'N': break
        idx = idxs[i]

        if i+1 < len(idxs):
            next_idx = idxs[i+1]
        else:
            print('idxs ended eary') 
            next_idx = -1
        

        pattern = kmer.replace('N', '.')
        regex = re.compile(pattern)

        if i==0: 
            val = '' 
            print('getting first state:', kmer, idx)
            pos_states = [i for i,x in enumerate(states) if re.match(regex, x)]
            pos_bits.append([(st, states[st], '') for st in pos_states])
            print('Fisrt state ps:', pos_states)
            pos_bits = pos_bits[0]

        adds = []
        for st_pkmer_bits in pos_bits:
            st, pkmer, bits = st_pkmer_bits[0], st_pkmer_bits[1], st_pkmer_bits[2]

            pmoves = np.nonzero(edges[st, :])[0]
            sorted_idx = np.argsort(scores[st, pmoves])
            pmoves = pmoves[sorted_idx]

            pattern = next_kmer.replace('N', '.')
            regex = re.compile(pattern)
            next_sts = [i for i, x in enumerate(states) if re.match(regex, x)]
        
            for next_st in next_sts:

                loc_idx = np.where(pmoves == next_st)[0]
                
                if loc_idx.size == 0:
                    continue

                bs = f'{loc_idx[0]:0{rate.numerator}b}'
                if (len(bs) != rate.numerator):
                    print('len(bits) != 2^p: seq {}, pos  {}, move_idx {}, {}->{}'.format(sq, i, loc_idx, st, next_st))
                    continue
                adds.append((next_st, states[next_st], bits+bs))
        if len(adds) == 0: 
            print('len(adds) = 0', 'i:', i, kmer, '->', next_kmer)
            print(pos_bits)
            if len(pos_bits)>0: val = pos_bits[0][2]
            return val
            
        pos_bits = adds

    val = pos_bits[0][2]
    print('END: val:', len(val), len(val[:seq_len]))
    return val


if __name__=='__main__':
    modif = '_nofpad'

    parser = argparse.ArgumentParser()
    parser.add_argument("-delta", type=int, default=0, help="Delta constraints to look at")
    parser.add_argument("-f_o", type=str, default=0, help="File path out")
    parser.add_argument("-f_p", type=str, default=0, help="File path padding")
    parser.add_argument("-f_c", type=str, default=0, help="File path check")
    args = parser.parse_args()
    d, f_o, f_p, f_c = args.delta, args.f_o, args.f_p, args.f_c

    np.random.seed(0)
    rates = ["2", "106", "106", "106", "106", "86", "86", "86", "86", "86", "66"]
    deltas = ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']
    for delta, r in zip([deltas[d]], [rates[d]]):

        print(d, delta, r, deltas[0])
        edges, states, scores = get_graph(6, int(delta)/100, max_run=2)
        edges = remove_nodes(edges)

        if delta == '000': rate = Fraction(2,1,_normalize=False)
        else: rate = Fraction(int(r[:-1]), int(r[-1]), _normalize=False)
        print('\nDELTA:', delta, 'RATE:', str(rate.numerator)+'/'+str(rate.denominator), '\n')

        path_var_out = f_o
        filename_basecalls = "basecalls"
        path_var_file = path_var_out + f"delta{delta}/" + filename_basecalls + '/'
        
        true_file = f_c + 'delta{}_rate{}.fasta'.format(delta, r)
        true_seqs = []
        with open(true_file, 'r') as f:
            for line in f: 
                if '>' in line: true_seqs.append(line)
        n_true = len(true_seqs)
        print('n_true:', n_true)

        out_folder = path_var_out+f'delta{delta}/basecalls_bits{modif}/'

        padding_file = f_p + 'delta{}_rate{}{}_fpad.txt'.format(str(int(delta)//100),str(rate.numerator),str(rate.denominator))
        with open(padding_file) as pf:
            pad = pf.readline()
        pad = pad.strip('\n')
        for j in range(n_true):

            for fol in range(1, 100+1):
                path_var_file_full = path_var_file+f'{fol}/j{j}'
                if not os.path.exists(path_var_file_full):
                    print('No file:', path_var_file_full)
                    continue
                path_var_file_full_idx = path_var_file+f'{fol}/j{j}_idxs'

                path_var_file_full_out = out_folder+f'{fol}/'
                if not os.path.exists(path_var_file_full_out):
                    os.makedirs(path_var_file_full_out)
                path_var_file_full_out = path_var_file_full_out+f'j{j}'
                        
                with open(path_var_file_full) as viterbi_file:
                    for line in viterbi_file:
                        viterbi_kmers = line
                    viterbi_file.close()   
                ch = ['\'', '[', ']', '2', '3', ' ']
                # Remove multiple characters from the string
                for character in ch:
                    viterbi_kmers = viterbi_kmers.replace(character, '')  
                
                viterbi_kmers = viterbi_kmers.split(',')  
                viterbi_seq = "" + viterbi_kmers[0]
                last_kmer = viterbi_kmers[0]
                viterbi_kmers_proc = [last_kmer]
                for b in viterbi_kmers[1:]:
                    if b == last_kmer: continue
                    else: 
                        viterbi_kmers_proc.append(b)
                        viterbi_seq=viterbi_seq+b[-1]
                        last_kmer=b
                viterbi_seq = viterbi_seq.strip('\n')
                seq = viterbi_seq

                with open(path_var_file_full_idx) as viterbi_file_idx:
                    for line in viterbi_file_idx:
                        viterbi_idx = line
                    viterbi_file_idx.close() 

                viterbi_idx = viterbi_idx.split(',')
                viterbi_idxs = [''.join(c for c in viterbi_idx[0] if c not in ' []')]
                last_idx = viterbi_idxs[0]
                for b in viterbi_idx[1:]:
                    b = ''.join(c for c in b if c not in ' []')
                    if b == last_idx: continue
                    else: 
                        viterbi_idxs.append(b)
                        last_idx=b
                idxs = viterbi_idxs

                ''' Fpad Removal ''' 
                if modif == '_nofpad':
                    begin_idx_vit = get_best_match(pad, seq)
                    seq = seq[begin_idx_vit-6:]
                    idxs = idxs[begin_idx_vit-6:]

                
                b_seq = seq2bits(seq, edges, states, scores, rate, sq=j, idxs=idxs)

                with open(path_var_file_full_out, 'w') as f:
                    f.write(b_seq)




