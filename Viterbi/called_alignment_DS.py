# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import sys
import numpy as np
import os
from pathlib2 import Path
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse


def format_cigar(cigar):

    pattern_ins = re.compile("(\d+)I")
    ins_counts = pattern_ins.findall(cigar)

    pattern_dels = re.compile("(\d+)D")
    dels_counts = pattern_dels.findall(cigar)

    pattern_subs = re.compile("(\d+)X")
    subs_counts = pattern_subs.findall(cigar)

    return sum(map(int, ins_counts)), sum(map(int, dels_counts)), sum(map(int, subs_counts))

if __name__=='__main__':
       
    parser = argparse.ArgumentParser()
    parser.add_argument("-f_o", help="Output file for seqs")
    parser.add_argument("-f_i", help="Input file for seqs")
    args = parser.parse_args() 
    
    np.random.seed(seed=0)
    f_o, f_i = args.f_o, args.f_i

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    n, bp = 1000, 200

    vit_mean, guppy_mean, ev_mean = [], [], []
    vit_std, guppy_std, ev_std = [], [], []
    vit_indelsub, guppy_indelsub, ev_indelsub = [], [], []

    rates = ["2", "106", "106", "106", "106", "106", "86", "86", "86", "86", "66"]
    deltas = ['000', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000']

    for delt, r in zip(deltas, rates):
        
        print('DELTA:', delt)
        cwd = os.getcwd()

        path_var = f_i
        path_var_long = path_var+ 'n{}_bp{}_delta{}_statesplit_rate{}_DeepSimu/'.format(n, bp, delt,r)

        score_f, fig_f, calls_f = 'scores/', '/Figs', 'basecalls/'

        path_var_out = f_o
        path_var_long_out = path_var_out+ 'n{}_bp{}_delta{}_DeepSimu/'.format(n, bp, delt)

        

        guppy_calls = {}
        filename_guppy = path_var_long+'pass.fastq'  
        with open(filename_guppy) as guppy_file:
            for record in SeqIO.parse(guppy_file, "fastq"):
                guppy_calls[record.id] = record.seq
                # print(record.id, record.seq)
            guppy_file.close() 

        mapping = {}
        for filename in os.listdir(path_var_long+'/fast5/'):
            filename_list = filename.split('_')
            if 'signal' not in filename_list:continue
            mapping[filename_list[-2]] = filename_list[-1].replace(".fast5", "")

        cs = guppy_calls.values()

        scores_dwll, scores_ev , scores_guppy , scores_raw = [], [], [], []
        lev_dwll, lev_ev , lev_guppy , lev_raw, lev_vit = [], [], [], [], []
        errs_vit, errs_ev , errs_guppy , errs_raw = [], [], [], []

        FP_all_vit, FN_all_vit = np.array([]), np.array([])
        FP_all_ont, FN_all_ont = np.array([]), np.array([])


        for j in mapping.keys():

            # print(j)
            ''' Guppy '''
            mapping_id = mapping[j]

            if mapping_id in guppy_calls.keys(): 
                seq_guppy = guppy_calls[mapping_id]
            else: 
                seq_guppy = cs[0]  	
            j = int(j)

            ''' True '''
            filename_true = path_var_long+'processed_genome_{}'.format(j)  
            with open(filename_true) as true_file:
                for line in true_file:
                    if '>' not in line: bases_true = line
                true_file.close()
            seq_true = Seq(bases_true)

            ''' Viterbi Dwell '''
            filename_viterbi_dwll = path_var_long_out+calls_f+'j{}'.format(j)
            # if not Path(filename_viterbi_dwll).is_file(): continue
            if Path(filename_viterbi_dwll).is_file():
                with open(filename_viterbi_dwll) as viterbi_file_dwll:
                    for line in viterbi_file_dwll:
                        viterbi_kmers_dwll = line
                    viterbi_file_dwll.close()   
                ch = ['\'', '[', ']', '2', '3', ' ']
                # Remove multiple characters from the string
                for character in ch:
                    viterbi_kmers_dwll = viterbi_kmers_dwll.replace(character, '')  
                
                viterbi_kmers_dwll = viterbi_kmers_dwll.split(',')
                    
                viterbi_seq_dwll = "" + viterbi_kmers_dwll[0]
                last_kmer = viterbi_kmers_dwll[0]
                viterbi_kmers_proc = [last_kmer]
                for b in viterbi_kmers_dwll[1:]:
                    if b == last_kmer: continue
                    else: 
                        viterbi_kmers_proc.append(b)
                        viterbi_seq_dwll=viterbi_seq_dwll+b[-1]
                        last_kmer=b
                    seq_dwll = Seq(viterbi_seq_dwll)

                file_out = path_var_long_out+calls_f+'j{}_processed'.format(j)
                with open(file_out, "w") as txt_file:
                        txt_file.write(str(viterbi_kmers_proc))
                        txt_file.close()
            else: continue 
            
            ''' Events '''
            filename_ev = path_var_long + 'fasta/signal_{}_{}.fasta'.format(j, mapping_id)  
            with open(filename_ev) as ev_file:
                for line in ev_file:
                    if '>' not in line: bases_ev = line
                ev_file.close()
            seq_ev = Seq(bases_ev)

            ''' Raw '''       
            filename_raw = path_var_long + 'fasta_raw/signal_{}_{}.fasta'.format(j, mapping_id)  
            if not Path(filename_raw).is_file(): continue
            with open(filename_raw) as raw_file:
                for line in raw_file:
                    if '>' not in line: bases_raw = line
                raw_file.close()
            seq_raw = Seq(bases_raw)

            
            ''' Align '''
            if Path(filename_viterbi_dwll).is_file():
                vit_align = edlib.align(str(seq_dwll), str(seq_true), task = "path")
                vit_alignment = edlib.getNiceAlignment(vit_align, str(seq_dwll), str(seq_true))
                lev_vit.append(vit_align["editDistance"])
                vit_ins, vit_dels, vit_subs = format_cigar(vit_align["cigar"])
                errs_vit.append((vit_ins, vit_dels, vit_subs))

            guppy_align = edlib.align(str(seq_guppy), str(seq_true), task = "path")
            guppy_alignment = edlib.getNiceAlignment(guppy_align, str(seq_guppy), str(seq_true))
            lev_guppy.append(guppy_align["editDistance"])
            guppy_ins, guppy_dels, guppy_subs = format_cigar(guppy_align['cigar'])
            errs_guppy.append((guppy_ins, guppy_dels, guppy_subs))

            ev_align = edlib.align(str(seq_ev), str(seq_true), task = "path")
            ev_alignment = edlib.getNiceAlignment(ev_align, str(seq_ev), str(seq_true))
            lev_ev.append(ev_align["editDistance"])
            ev_ins, ev_dels, ev_subs = format_cigar(ev_align["cigar"])
            errs_ev.append((ev_ins, ev_dels, ev_subs))

            raw_align = edlib.align(str(seq_raw), str(seq_true), task = "path")
            raw_alignment = edlib.getNiceAlignment(raw_align, str(seq_raw), str(seq_true))
            lev_raw.append(raw_align["editDistance"])
            raw_ins, raw_dels, raw_subs = format_cigar(raw_align["cigar"])
            errs_raw.append((raw_ins, raw_dels, raw_subs))

            filepath_align = path_var_long_out+score_f
            if not os.path.exists(filepath_align):
                os.makedirs(filepath_align)

            filename_align = filepath_align+'align_j{}.txt'.format(j)
            with open(filename_align, "w") as txt_file:
                if Path(filename_viterbi_dwll).is_file():
                    txt_file.write('Viterbi: \n')
                    txt_file.write("\n".join(vit_alignment.values())) 
                    txt_file.write('\n')
                    txt_file.write('Lev Dist: {}'.format(vit_align["editDistance"]))
                    txt_file.write('\n')
                    txt_file.write('Edits: {}'.format(format_cigar(vit_align["cigar"])))
                    txt_file.write('\n')
                    txt_file.write('\n')
                txt_file.write('Guppy: \n')
                txt_file.write("\n".join(guppy_alignment.values())) 
                txt_file.write('\n')
                txt_file.write('Lev Dist: {}'.format(guppy_align["editDistance"]))
                txt_file.write('\n')
                txt_file.write('Edits: {}'.format(guppy_align["cigar"]))
                txt_file.write('\n')
                txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(guppy_ins, guppy_dels, guppy_subs))
                txt_file.write('\n')
                txt_file.write('\n')
                txt_file.write('Raw: \n')
                txt_file.write("\n".join(raw_alignment.values()))  
                txt_file.write('\n')
                txt_file.write('Lev Dist: {}'.format(raw_align["editDistance"]))
                txt_file.write('\n')
                txt_file.write('Edits: {}'.format(raw_align["cigar"]))
                txt_file.write('\n')
                txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(raw_ins, raw_dels, raw_subs))
                txt_file.write('\n')
                txt_file.write('\n')
                txt_file.write('Event: \n')
                txt_file.write("\n".join(ev_alignment.values())) 
                txt_file.write('\n')
                txt_file.write('Lev Dist: {}'.format(ev_align["editDistance"])) 
                txt_file.write('\n')
                txt_file.write('Edits: {}'.format(ev_align["cigar"]))
                txt_file.write('\n')
                txt_file.write('\n')
                txt_file.write('Insertions: {}, Deletions: {}, Substitutions: {}'.format(ev_ins, ev_dels, ev_subs))
                txt_file.close()

        avg_err_vit = np.mean(errs_vit, axis=0)

        avg_err_guppy = np.mean(errs_guppy, axis=0)
        avg_err_raw = np.mean(errs_raw, axis=0)
        avg_err_ev = np.mean(errs_ev, axis=0)

        filename_scores = filepath_align+'/alignment_scores.txt'
        with open(filename_scores, "w") as txt_file:
            txt_file.write('Viterbi: \n')
            txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_vit), np.std(lev_vit)))
            txt_file.write('\n')
            txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_vit[0], avg_err_vit[1], avg_err_vit[2]))
            txt_file.write('\n')
            txt_file.write('\n')
            txt_file.write('Guppy: \n')
            txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_guppy), np.std(lev_guppy)))
            txt_file.write('\n')
            txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_guppy[0], avg_err_guppy[1], avg_err_guppy[2]))
            txt_file.write('\n')
            txt_file.write('\n')
            txt_file.write('Raw: \n')
            txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_raw), np.std(lev_raw)))
            txt_file.write('\n')
            txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_raw[0], avg_err_raw[1], avg_err_raw[2]))
            txt_file.write('\n')
            txt_file.write('\n')
            txt_file.write('Event: \n')
            txt_file.write('Levenshtein Avg. Distance: {}, std: {}'.format(np.mean(lev_ev), np.std(lev_ev)))
            txt_file.write('\n')
            txt_file.write('Avg. Insertions: {}, Deletions: {}, Substitutions: {}'.format(avg_err_ev[0], avg_err_ev[1], avg_err_ev[2]))
            txt_file.write('\n')
            txt_file.close() 

        vit_mean.append(np.mean(lev_vit))
        guppy_mean.append(np.mean(lev_guppy)) 
        ev_mean.append(np.mean(lev_ev))

        vit_std.append(np.std(lev_vit))
        guppy_std.append(np.std(lev_guppy)) 
        ev_std.append(np.std(lev_ev))

        vit_indelsub.append(avg_err_vit)  
        guppy_indelsub.append(avg_err_guppy)  
        ev_indelsub.append(avg_err_ev) 

 

        
    
