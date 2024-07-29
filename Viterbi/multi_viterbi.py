# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import numpy as np
import argparse
import os, os.path
import viterbi as viterbi
import multiprocessing
import datetime

def process_sig_diff(signal, cutoff=4):

    distances = np.abs(np.diff(signal))
    trans_idx = np.where(distances>cutoff)[0]
    delete_idx = []
    for i in trans_idx:
        if (i+1 in trans_idx) and (np.abs(signal[i]-signal[i+2])>cutoff):
            delete_idx.append(i+1) 

    signal_processed = np.delete(signal,delete_idx)

    return signal_processed

def process_sig_seg(signal, alig_file, remove_num=1):
    alig = np.loadtxt(alig_file, dtype=int)
    alig = alig[:, 1] - 1 

    alig_diff = np.diff(alig)
    trans_idx = np.where(alig_diff!=0)[0]
    delete_idx = np.concatenate((trans_idx-1, trans_idx), dtype=int)
    signal_processed = np.delete(signal, delete_idx)

    return signal_processed 

def process_sig_line(signal, cutoff=5, err=15, win=4):
    times = np.arange(signal.size)
    delete_idx = []
    for i in range(len(signal) - win + 1):
        ys = signal[i:i+win]
        xs = times[i:i+win]
        
        vals = np.polyfit(xs, ys, 1, full=True)
        slope = vals[0][0]
        res = vals[1][0]
        
        p = np.poly1d(vals[0])
        
        if res < err and np.abs(slope) > cutoff:
            delete_idx.extend(range(i+1, i+win-1))    
    signal_processed = np.delete(signal,delete_idx)
    return signal_processed

# produces viterbi calls (with multiprocessing) for deep simulation currents
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-maxrun", type=int, default=0, help="Max length of allowable homopolymer runs")
    parser.add_argument("-b", type=int, default=200, help="Read rength")
    parser.add_argument("-dwell", type=int, default=10, help="Self-Transition mean dwell time")
    parser.add_argument("-n", type=int, default=1000,help="Num. reads")
    parser.add_argument("-k", type=int, default=6, help="k-mer size for de Bruijn Graph")
    parser.add_argument("-f_o", help="Output file location for basecalls")
    parser.add_argument("-f_b", help="Output file name for basecalls")
    parser.add_argument("-f_i", help="Input file location for seqs")
    parser.add_argument("-f_s", help="k-mer mean file location")
    parser.add_argument("-d", type=str, nargs='+', help="deltas")
    parser.add_argument("-r", type=str, default=None, nargs='+', help="rates")
    parser.add_argument("-procc", choices=['none', 'seg', 'diff', 'line'], default='none', help="Process input seq by removing trans. points")
    # parser.add_argument("-r_type", choices=['single', 'repeats'], default='single', help="Does data set contain single reads or multiple repeats?")
    args = parser.parse_args()
    maxrun, procc, dwell, bp, n, k, path_var_out, filename_basecalls, path_var, score_file, deltas, rates = args.maxrun, args.procc, args.dwell, args.b, args.n, args.k, args.f_o, args.f_b, args.f_i, args.f_s, [d+'00' for d in args.d], args.r
    state_split=True
    if rates is None: state_split=False
    
    file_check = path_var+f'/delta{deltas[0]}_rate{rates[0]}/1/signal/'
    n_true = len([name for name in os.listdir(file_check) if os.path.isfile(file_check+name)])
    fols_start, fols_end = 25+1, 50 #50+1, 75 #75+1, 100 #1, 25 #

    # kmer means from score file
    bn = np.loadtxt(score_file, skiprows=1, usecols=[2])

    print('max run:', maxrun)
    print('Num. Seqs:', n_true)

    # n_start, n_end = 0, n_true
    # n_1s, n_2s = np.arange(n_start, n_end-batch_size+1, batch_size), np.arange(batch_size+n_start, n_end+1, batch_size)#
    # print(batch_size, n_true, n_1s, n_2s)
    # remainder = (n_end-n_start)%batch_size
    # if remainder>0:
    #     n_1s = np.append(n_1s, n_2s[-1])
    #     n_2s = np.append(n_2s, n_2s[-1]+remainder)
    #     print(n_1s, n_2s)
    
    # for (n1, n2) in zip(n_1s, n_2s): 
    for j in range(n_true):
        file_check = path_var_out+f"delta{deltas[0]}/"+filename_basecalls+f"/{fols_start}/j{j}"
        if not os.path.exists(file_check):
            print("j{} for fols {}-{}".format(j, fols_start, fols_end))

            template_data = lambda *, delta, r, fol:f"delta{delta}_rate{r}/{fol}/signal/signal_{j}.txt" #lambda *, delta, r, j:f"delta{delta}_rate{r}/signal/signal_{j}.txt"#
            filenames_data = [template_data(delta=delta, r=r, fol=fol) for fol in range(fols_start, fols_end+1) for delta,r in zip(deltas, rates)]#[template_data(delta=delta, r=r, j=j) for delta,r in zip(deltas, rates) for j in range(n1, n2)]#
            
            template_alig = lambda *, delta, r, fol:f"delta{delta}_rate{r}/{fol}/align/align_{j}.ali"#lambda *, delta, r, j:f"delta{delta}_rate{r}/align/align_{j}.ali"#
            filenames_alig = [template_alig(delta=delta, r=r, fol=fol) for fol in range(fols_start, fols_end+1) for delta,r in zip(deltas, rates)]#[template_alig(delta=delta, r=r, j=j) for delta,r in zip(deltas, rates) for j in range(n1, n2)]#

            folders_out = [path_var_out+f"delta{delta}/"+filename_basecalls+f"/{fol}/" for fol in range(fols_start, fols_end+1) for delta in deltas]#[path_var_out+ f"delta{delta}/" for delta in deltas]#

            # Vitebri inputs for multiprocessing
            if procc=='seg':
                datas = [process_sig_seg((np.loadtxt(path_var+fb)-14)/5.7, path_var+fa) for fb, fa in zip(filenames_data, filenames_alig)]
            elif procc=='diff':
                datas = [process_sig_diff((np.loadtxt(path_var+fb)-14)/5.7) for fb in filenames_data]
            elif procc=='line':
                datas = [process_sig_line((np.loadtxt(path_var+fb)-14)/5.7) for fb in filenames_data]
            else:
                datas = [(np.loadtxt(path_var+fb)-14)/5.7 for fb in filenames_data]
            models = [viterbi.model(int(k), int(delta)/100, bn, score_file, dwell, max_run=maxrun, score='mean') for fol in range(fols_start, fols_end+1) for delta in deltas]#[viterbi.model(int(k), int(delta)/100, bn, score_file, dwell, score='mean') for delta in deltas for j in range(n1, n2)]#

            with multiprocessing.Pool(processes=len(filenames_data)) as pool:#batch_size) as pool:
                base_calls, qs_idx = zip(*pool.starmap(viterbi.viterbi, zip(datas, models)))
                pool.close()
                pool.join()

            # save results
            for out_path in folders_out:
                if not os.path.exists(out_path):
                    os.makedirs(out_path)

            files_out = [out_p+f"j{j}" for out_p in folders_out]#[out_p+f"j{j}" for out_p in folders_out for j in range(n1, n2)]#
            files_out_idxs = [out_p+f"j{j}_idxs" for out_p in folders_out]
            time_end = datetime.datetime.now()
            
            for i, file_out_list in enumerate(zip(files_out, files_out_idxs)):
                file_out, file_out_idx = file_out_list[0], file_out_list[1]
                print(file_out, filenames_data[i])

                with open(file_out_idx, "w") as txt_file:
                    txt_file.write(str(qs_idx[i]))
                    txt_file.close()

                with open(file_out, "w") as txt_file:
                    txt_file.write(str(base_calls[i]))
                    txt_file.close()
        else: print("ALREADY COMPLETE: j{} for fols {}-{}".format(j, fols_start, fols_end))





    
