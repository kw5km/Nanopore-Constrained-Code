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
from PIL import Image
import math
from pathlib2 import Path
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
import difflib
import argparse

def bitstring_to_bytes(s):
    bits = np.array([int(b) for b in s])
    bytes = np.packbits(bits)
    return bytes
    # byte_hold = int(s, 2).to_bytes((len(s) + 7) // 8, byteorder='big')
    # return byte_hold#bytearray(byte_hold)

def get_graph(k, delta, max_run=2, f_s='/u/kw5km/Research/Nanopore_Channel/Regression_Output/6mer_means.txt', r='2'):
    delt = str(int(delta))
    # print('/u/kw5km/Research/Nanopore_Channel/Output/graphs/Aq_delta{}_rate{}_3run.npy'.format(delt, r))
    if os.path.exists('/u/kw5km/Research/Nanopore_Channel/Output/graphs/Aq_delta{}_rate{}_3run.npy'.format(delt, r)):
        print('Existing Graph Loaded')
        A_q = np.load('/u/kw5km/Research/Nanopore_Channel/Output/graphs/Aq_delta{}_rate{}_3run.npy'.format(delt, r))
        states_split = np.load('/u/kw5km/Research/Nanopore_Channel/Output/graphs/states_delta{}_rate{}_3run.npy'.format(delt, r))
        scores2 = np.load('/u/kw5km/Research/Nanopore_Channel/Output/graphs/scores_delta{}_rate{}_3run.npy'.format(delt, r))

    else:

        print('Generating State Split Graph')
        set1 = ['A', 'C', 'G', 'T']
        kmers = deBruijn_helper.get_all_kmers(k, set1)

        edges, edges_set = deBruijn_helper.get_debruijn_edges_from_kmers(kmers, max_run=max_run)
        scores = deBruijn_helper.get_scores_means(edges, kmers, score_file=f_s)

        mask = [scores<delta][0]
        edges[mask] = 0
        A = edges

        if delta == 0:
            np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/Aq_delta{}_rate{}_3run'.format(delt, r), A)
            np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/states_delta{}_rate{}_3run'.format(delt, r), kmers)
            np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/scores_delta{}_rate{}_3run'.format(delt, r), scores)
            return A, kmers, scores
        
        delt += '00'
        denom = k
        cap = deBruijn_seqs.get_rate(A)
        print('Capacity:', cap)
        cap = int(cap * (10 ** 2)) / 10 ** 2
        rate = Fraction(int(np.floor(denom*cap)),denom, _normalize=False)

        num, denom = rate.numerator, rate.denominator
        r = str(rate.numerator)+str(rate.denominator)
        if (num % 2) != 0: num-=1
        rate = Fraction(num, denom, _normalize=False)
        print('Rate Calculated:', rate.numerator, rate.denominator, rate)

        A_q, states_split = state_splitting.state_split(A, num, rate.denominator, kmers)
        scores2 = deBruijn_helper.get_scores_means(A_q, states_split, score_file=f_s)

        
        np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/Aq_delta{}_rate{}_3run_noloop'.format(delt, r), A_q)
        np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/states_delta{}_rate{}_3run_noloop'.format(delt, r), states_split)
        np.save('/u/kw5km/Research/Nanopore_Channel/Output/graphs/scores_delta{}_rate{}_3run_noloop'.format(delt, r), scores2)
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

def seq2bits(seq, edges, states, scores, rate, k=6, sq=0, idxs=None, delta=6): 
    seq_len = int(186*(rate.numerator/rate.denominator))
    test_file = f'/u/kw5km/Research/Nanopore_Channel/Data_img_v2/messages/delta{delta}_rate{rate.numerator}{rate.denominator}_testconsensus.txt'

    with open(test_file, "a+") as f:
        f.write('\nseq {}:\n'.format(sq))

    pos_bits = []
    kmers = squiggle_helper.rolling_window(np.array(list(seq)), k, overlap=rate.denominator)
    # print('kmer len:', len(kmers))
    for i, kmer in enumerate(kmers[:-1]):
        next_kmer = kmers[i+1]

        if next_kmer == 'N': break
        # idx = idxs[i]

        # if i+1 < len(idxs):
        #     next_idx = idxs[i+1]
        # else:
        #     print('idxs ended eary') 
        #     next_idx = -1
        

        pattern = kmer.replace('N', '.')
        regex = re.compile(pattern)

        if i==0: 
            val = ''
            # print('getting first state:', kmer, idx)
            pos_states = [i for i,x in enumerate(states) if re.match(regex, x)]
            pos_bits.append([(st, states[st], '') for st in pos_states])
            # print('Fisrt state ps:', pos_states)
            pos_bits = pos_bits[0]
            # print(pos_bits)

        adds = []
        check = 0
        for st_pkmer_bits in pos_bits:
            check +=1
            print('check:', check)
            st, pkmer, bits = st_pkmer_bits[0], st_pkmer_bits[1], st_pkmer_bits[2]
            if check>=5000:
                if len(adds) == 0: 
                    if len(pos_bits)>0: val = pos_bits[0][2]
                    else: bits = '' 
                    return val 
                return adds[0][2]
            # print('st_pkmer_bits:', st_pkmer_bits, len(bits))
            # pos_bits.remove(st_pkmer_bits)
            pmoves = np.nonzero(edges[st, :])[0]
            sorted_idx = np.argsort(scores[st, pmoves])
            pmoves = pmoves[sorted_idx]
            # print('pmove', pmoves)

            pattern = next_kmer.replace('N', '.')
            regex = re.compile(pattern)
            next_sts = [i for i, x in enumerate(states) if re.match(regex, x)]

            # print('next_sts:', next_sts)
            for next_st in next_sts:
                # print('next st:', next_st)
                loc_idx = np.where(pmoves == next_st)[0]
                # print('loc idx:', loc_idx)
                
                if (loc_idx.size == 0) or (edges[st, next_st].size == 0):
                    # print('SKIPPED TRANSITION:', st, '->', next_st)
                    # print(edges[st, next_st])
                    continue

                bs = f'{loc_idx[0]:0{rate.numerator}b}'
                if len(bs) != rate.numerator:
                    # print('len(bits) != 2^p: seq {}, pos {}, move_idx {}, {}->{}'.format(sq, i, loc_idx, st, next_st))
                    continue
                # print('kmer:', i, idx, '->', next_idx, 'states[next_idx]:', states[int(next_idx)], 'Next State:', next_st, states[next_st]) 
                adds.append((next_st, states[next_st], bits+bs))
        if len(adds) == 0: 
            # print('len(adds) = 0', 'i:', i, kmer, '->', next_kmer)
            # print(pos_bits)
            if len(pos_bits)>0: val = pos_bits[0][2] 
            return val#[:seq_len]
            
        pos_bits = adds

    val = pos_bits[0][2]
    print('END: val:', len(val), len(val[:seq_len]))
    return val#[:seq_len]

def get_best_match(query, corpus, step=2, flex=5, case_sensitive=False):
    """Return best matching substring of corpus.

    Parameters
    ----------
    query : str
    corpus : str
    step : int
        Step size of first match-value scan through corpus. Can be thought of
        as a sort of "scan resolution". Should not exceed length of query.
    flex : int
        Max. left/right substring position adjustment value. Should not
        exceed length of query / 2.

    Outputs
    -------
    output0 : str
        Best matching substring.
    output1 : float
        Match ratio of best matching substring. 1 is perfect match.

    Source : https://stackoverflow.com/questions/36013295/find-best-substring-match
    """

    def _match(a, b):
        """Compact alias for SequenceMatcher."""
        return difflib.SequenceMatcher(None, a, b).ratio()

    def scan_corpus(step):
        """Return list of match values from corpus-wide scan."""
        match_values = []

        m = 0
        while m + qlen - step <= len(corpus):
            match_values.append(_match(query, corpus[m : m-1+qlen]))
            m += step

        return match_values

    def index_max(v):
        """Return index of max value."""
        return max(range(len(v)), key=v.__getitem__)

    def adjust_left_right_positions():
        """Return left/right positions for best string match."""

        p_l, bp_l = [pos] * 2
        p_r, bp_r = [pos + qlen] * 2

        bmv_l = match_values[p_l // step]
        bmv_r = match_values[p_l // step]

        for f in range(flex):
            ll = _match(query, corpus[p_l - f: p_r])
            if ll > bmv_l:
                bmv_l = ll
                bp_l = p_l - f

            lr = _match(query, corpus[p_l + f: p_r])
            if lr > bmv_l:
                bmv_l = lr
                bp_l = p_l + f

            rl = _match(query, corpus[p_l: p_r - f])
            if rl > bmv_r:
                bmv_r = rl
                bp_r = p_r - f

            rr = _match(query, corpus[p_l: p_r + f])
            if rr > bmv_r:
                bmv_r = rr
                bp_r = p_r + f
        return bp_l, bp_r, _match(query, corpus[bp_l : bp_r])

    if not case_sensitive:
        query = query.lower()
        corpus = corpus.lower()
    qlen = len(query)

    if flex >= qlen/2:
        flex = 3

    match_values = scan_corpus(step)
    pos = index_max(match_values) * step

    pos_left, pos_right, match_value = adjust_left_right_positions()
    if pos_right>20: pos_right = 12
    return pos_right

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-delta", type=int, default=0, help="Delta constraints to look at")
    parser.add_argument("-f_o", type=str, default=0, help="File path out")
    parser.add_argument("-f_c", type=str, default=0, help="File path check")
    parser.add_argument("-f_p", type=str, default=0, help="File path pad")
    args = parser.parse_args()
    d, f_o, f_c, f_p = args.delta, args.f_o, args.f_c, args.f_p

    np.random.seed(0)
    rates = ["2", "106", "106", "106", "106", "86", "86", "86", "86", "86", "66"]
    deltas = ["000", '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']
    read_depths = np.concatenate(([3,5], np.arange(10, 50+1, 10)))

    modif = '_constraint_nofpad'

    for delta, r in zip([deltas[d]], [rates[d]]):
        bp = 186

        print('Delta:', delta)
        edges, states, scores = get_graph(6, int(delta)/100, max_run=2, r=r)
        if delta == '000': rate = Fraction(2,1,_normalize=False)
        else: rate = Fraction(int(r[:-1]), int(r[-1]), _normalize=False)

        path_var_out = f_o
        filename_consensus = "consensus{}".format(modif)
        path_var_file = path_var_out + f"delta{delta}/" + filename_consensus + '/'

        file_check = f_c + "delta{}_rate{}.fasta".format(delta, r)
        f_checks = []
        with open(file_check, 'r') as f:
            for line in f: 
                if '>' in line: f_checks.append(line)
        n_true = len(f_checks)

        seq_len = seq_len = int(bp*(rate.numerator/rate.denominator))
        print('n_true:', n_true, 'len:', seq_len)

        for read_depth in read_depths:
            print(read_depth)
            out_folder = path_var_out+f'delta{delta}/consensus_bits{modif}/read_depth{read_depth}/'
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            seqs = []
            b_seqs = []
            for j in range(n_true):

                consensus_file = path_var_file+'read_depth{}/j{}_consensus_bio.fasta'.format(read_depth, j)
                if not os.path.exists(consensus_file): 
                    print('skipping file ', consensus_file)
                    seqs.append('skip')
                    continue

                if r == '2': r ='21'
                padding_file = f_p + '/delta{}_rate{}_fpad.txt'.format(str(int(delta)//100),r)
                with open(padding_file) as pf:
                    pad = pf.readline()
                pad = pad.strip('\n')

####
                fasta_consensus = SeqIO.parse(open(consensus_file),'fasta')
                for fasta_c in fasta_consensus:
                    fasta_cons_str = str(fasta_c.seq.ungap("-"))
                    seq = fasta_cons_str

                '''Remove fpad'''
                if (modif == '_constraint_fpad') or (modif == '_constraint'):
                    begin_idx_vit = get_best_match(pad, seq)
                    seq = seq[begin_idx_vit-6:]

                if int(delta)//100 == 3: print(j, seq)

                print('delta:', delta, 'j:', j, 'seq:', seq)
                seqs.append(seq)

                b_seq = seq2bits(seq, edges, states, scores, rate, sq=j)
                b_seqs.append(b_seq)
                
                bit_file = out_folder+'j{}_nopad.fasta'.format(j)
                with open(bit_file, 'w') as f:
                    f.write('>\n')
                    f.write(str(b_seq)+'\n')
                    
        