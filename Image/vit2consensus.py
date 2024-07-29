# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from itertools import chain
import os, subprocess
import numpy as np
import argparse
from fractions import Fraction
from best_match import get_best_match


def majority_vote(seqs):
    seq = [] 
    max_seq_len = max([len(s) for s in seqs])
    seqs_list = [list(s)+['']*max(0, max_seq_len-len(s)) for s in seqs]
    seq_array = np.array(seqs_list)
    for base in range(seq_array.shape[1]):
        bases = seq_array[:, base]
        vals, counts = np.unique(bases, return_counts=True)
        ind = np.argmax(counts)
        seq.append(vals[ind])
    seq = ''.join(seq)
    return seq

def multi_align(file_in, j):
    file_out = file_in.strip('.fasta')+'_multi_align3.fasta'
    subprocess.call('muscle -align %s -output %s'%(file_in,file_out), shell=True)
    return file_out

def consensus(file_path, j):
    align_file = multi_align(file_path, j)
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(align_file, "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.gap_consensus(0.01, "N")#, require_multiple=True)#.dumb_consensus(0.01, "N")#
    return consensus



if __name__=='__main__':
    modif = '_nofpad'
    if modif == '_nofpad': sub_mult = 5
    else: sub_mult = 3

    parser = argparse.ArgumentParser()
    parser.add_argument("-delta", type=int, default=0, help="Delta constraints to look at")
    parser.add_argument("-f", type=str, default=0, help="File path")
    parser.add_argument("-f_p", type=str, default=0, help="File path")
    args = parser.parse_args()
    d, f, f_p = args.delta, args.f, args.f_p

    np.random.seed(0)
    rates = ["2", "106", "106", "106", "106", "86", "86", "86", "86", "86", "66"]
    deltas = ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']
    bp = 186
    read_depths = np.ndarray.flatten(np.concatenate(([3,5], np.arange(10,  100+1, 10))))
    
    for delta, r in zip([deltas[d]], [rates[d]]):

        for read_depth in read_depths:
            print('DELTA:', delta, 'READ DEPTH:', read_depth)

            if delta == '000': 
                rate = Fraction(2,1,_normalize=False)
                num = 2*6
            else:
                rate = Fraction(int(r[:-1]), int(r[-1]), _normalize=False)
                num = rate.numerator 

            path_var = f
            path_var_out = path_var + "output/"
            filename_basecalls = "basecalls"
            path_var_file = path_var_out + f"delta{delta}/" + filename_basecalls + '/'

            file_check = path_var + "delta{}_rate{}.fasta".format(delta, r)
            f_checks = []
            with open(file_check, 'r') as f:
                for line in f: 
                    if '>' in line: f_checks.append(line)
            n_true = len(f_checks)
            print('n_true:', n_true)

            bit_check = path_var + "messages/delta{}_rate{}_message.txt".format(int(delta)//100, r)
            bits_check = np.loadtxt(bit_check, delimiter=',')
            bits_check_len = bits_check.shape[1]
            print('bits_check_len:', bits_check_len)
            bits_check_len = (186+12+18)*rate
            print('bits_check_len:', bits_check_len)

            out_folder = path_var_out+f'delta{delta}/consensus_constraint{modif}/read_depth{read_depth}/'
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            padding_file = f_p + 'delta{}_rate{}{}_fpad.txt'.format(str(int(delta)//100), str(rate.numerator), str(rate.denominator))
            with open(padding_file) as pf:
                pad = pf.readline()
            pad = pad.strip('\n')

            for j in range(n_true):
                print('delta, j:', delta, j)
                 
                seqs = []
                # seqs_guppy = []
                
                out_file, out_file_guppy = out_folder+f'j{j}_consensus_bio.fasta', out_folder+f'j{j}_consensus_bio_guppy.fasta'
                in_file, in_file_guppy = out_folder+f'j{j}.fasta', out_folder+f'j{j}_guppy.fasta'
                print('in file:', in_file)
                print('out file:', out_file)

                possible_fols = os.listdir(path_var_file)

                possible = []
                all_fols = []
                for pf in possible_fols:
                    file_check = path_var_file+f'/{pf}/j{j}'
                    if os.path.exists(file_check):
                        all_fols.append(pf)
                        if modif == '': possible.append(pf)
                        else:
                            bits_check = path_var_out + f"delta{delta}/" + f'basecalls_bits{modif}/{pf}/j{j}'
                            b_seq = ''
                            if os.path.exists(bits_check):
                                with open(bits_check, 'r') as bits_seqs:
                                    for bs in bits_seqs: b_seq = str(bs)
                            if len(b_seq)>=bits_check_len-(num*sub_mult):possible.append(pf)

                pos_check = True
                print('len(possible):', len(possible))
                print('possible:', possible)
                if len(possible) < read_depth:
                    if len(possible) == 0: 
                        pos_check = False
                        # check = False
                        continue
                    print('Limiting j:', j, 'to read depth:', len(possible))
                    # check = False

                # print('possible fols:', possible)
                selected_fols = np.random.choice(all_fols, read_depth, replace=False)#np.random.choice(possible, read_depth, replace=False)#
                print('selected fol:', selected_fols)
                fols = [f for f in selected_fols if f in possible]
                print('selected and possible fol:', fols)

                if len(fols) == 0: pos_check = False
                if not pos_check: continue

                for fol in fols:
                    # len(guppy_calls))
                    print('folder', fol)

                    path_var_file_full = path_var_file+f'{fol}/j{j}'
                    
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

                    '''Remove fpad'''
                    if modif == '_nofpad':
                        begin_idx_vit = get_best_match(pad, viterbi_seq)
                        viterbi_seq = viterbi_seq[begin_idx_vit-6:]
                    '''           '''

                    seq = Seq(viterbi_seq)
                    seqs.append(viterbi_seq)

                with open(in_file, 'w') as f:
                    for s, fol in zip(seqs, fols):
                        print('write fol:', fol)
                        f.write(f'>{fol}\n')
                        f.write(s+'\n')

                consensus_seq = consensus(in_file, j)
                with open(out_file, 'w') as f:
                    f.write('>\n')
                    f.write(str(consensus_seq)+'\n')
