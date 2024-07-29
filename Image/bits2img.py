# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import numpy as np
import os
from fractions import Fraction
from PIL import Image
import math
from Bio import SeqIO
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt
import argparse

def bitstring_to_bytes(s):
    bits = np.array([int(b) for b in s])
    bytes = np.packbits(bits)
    return bytes

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f_o", help="Output file for seqs")
    parser.add_argument("-f_c", help="Input file for check")
    args = parser.parse_args()
    f_o, f_c = args.f_o, args.f_c

    modif, modif2 = '_constraint_nofpad', '_nopad'
    np.random.seed(0)
    rates = ["2", "106", "106", "106", "106", "86", "86", "86", "86"]#
    deltas = ["000", '100', '200', '300', '400', '500', '600', '700', '800']#
    bp = 186
    read_depths = np.concatenate(([3],[30]))#np.concatenate(([3,5], np.arange(10, 100+1, 10)))
    rd_3, rd_30 = [], []

    for delta, r in zip(deltas, rates):

        if delta == '000': rate = Fraction(2,1,_normalize=False)
        else: rate = Fraction(int(r[:-1]), int(r[-1]), _normalize=False)
        print('\nDELTA:', delta, 'RATE:', str(rate.numerator)+'/'+str(rate.denominator), '\n')

        message_file = f_c + 'delta{}_rate{}_message.txt'.format(int(delta)//100, r)
        message = np.loadtxt(message_file, delimiter=',')
        bit_len = 186*(rate.numerator/rate.denominator) - 8

        seq_len = int(bp*(rate.numerator/rate.denominator)) - 8
        num_seqs = math.ceil((64*64*8)/seq_len)

        messages = np.reshape(message, (-1, num_seqs))
        print('message size:', message.shape)

        
        for rd in read_depths:
            print('RD:', rd)
            print('\nREAD DEPTH:', rd, '\n')
            out_folder = f_o + 'delta{}/'.format(delta)
            fig_file = out_folder + 'img_read_depth{}{}_noRS_v2.png'.format(rd, modif2)
            fig_file2 = out_folder + 'img_read_depth{}{}_test.png'.format(rd, modif2)
            fig_file3 = out_folder + 'img_read_depth{}{}_diff.png'.format(rd, modif2)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)

            seqs = []
            bits = ''
            ids = []
            error_vec = []
            for j in range(num_seqs):

                bits_file = f_o + 'delta{}/consensus_bits{}/read_depth{}/j{}{}.fasta'.format(delta, modif, rd, j, modif2)
                if not os.path.exists(bits_file): 
                    seqs.append('skip')
                    error_vec.append(1)
                    continue
                else: error_vec.append(0)

                fasta_consensus = SeqIO.parse(open(bits_file),'fasta')
                for fasta_c in fasta_consensus:
                    fasta_cons_str = str(fasta_c.seq.ungap("-"))
                    seq = fasta_cons_str
                    seq = seq.replace('N', '0')
                    seqs.append(seq)

            print('seqs:', len(seqs))
            for sq, seq in enumerate(seqs):
                if (seq == 'skip') or (len(seq)==0): 
                    b_seq = '0'*(seq_len)
                    ids.append[-1] 
                else:
                    seq = seq.strip('\n')
                    ids.append(int(seq[:8]))
                    b_seq = seq[8:seq_len]
                    if len(b_seq)<(seq_len): 
                        b_seq += '0'*(seq_len-len(b_seq))
                bits += b_seq

            print('bits len:', len(bits), type(bits))
            message_test = ''.join(str(int(x)) for x in messages.flatten())
            message_test = message_test[:64*64*8]
            print('message_test len:', len(message_test), type(message_test))            

            img_data = {'data': bits+'0'*(num_seqs*int(bit_len)-len(bits))}
            erasure_data = {'erasure_idxs': error_vec}
            savemat(f_o + 'delta{}/img_data_rd{}_noers.mat'.format(delta, rd), img_data)
            savemat(f_o + 'delta{}/erasure_data_rd{}.mat'.format(delta, rd), erasure_data)

            if len(bits)<64*64*8:
                bits += '0'*(64*64*8 - len(bits))
            else:
                bits = bits[:64*64*8]

            bytesarr = bitstring_to_bytes(bits)
            bytesarr = bytesarr.reshape((64, 64)).T

            im = Image.new('L', (64, 64))
            image = im.load()

            indices = np.array(ids)
            for pixel in range(64*64):
                indxs = np.where(indices == pixel)[0]
                if indxs.size == 0:
                    i = ids.index(-1)
                else:
                    i = indxs[0]
                image[pixel] = int(bytesarr[int(i,2)])

            messageArr = bitstring_to_bytes(message_test)
            messageArr = messageArr.reshape((64, 64)).T
            im2 = Image.new('L', (64, 64))
            image2 = im2.load()
            for row in range(64):
                for col in range(64):
                    image2[row, col] = int(messageArr[row, col])

            im.save(fig_file)
            im2.save(fig_file2)

            count = sum(1 for c1, c2 in zip(bits, message_test) if c1 != c2)
            count_perct = count/(64*64*8)
            print('different indices:', count, count_perct)
            if rd==3: rd_3.append(count_perct)
            elif rd==30: rd_30.append(count_perct)

            errs = []
            for b, m in zip(bits, message_test):
                if b==m: errs.append(0)
                else: errs.append(1)
            errs = np.array(errs)
            errs = errs.reshape((64, 64*8))
            print(errs.shape)
            plt.subplot(111)
            plt.imshow(errs, cmap='Greys',  interpolation='nearest')
            plt.show()
            plt.savefig(fig_file3)
            plt.close()

    # headers = [r"Reads = 3", r"Reads = 30"]
    # data = dict()
    # data[r"$\delta=0$"] = [format(rd_3[0], '.3f'), format(rd_30[0], '.3f')]
    # data[r"$\delta=1$"] = [format(rd_3[1], '.3f'), format(rd_30[1], '.3f')]
    # data[r"$\delta=2$"] = [format(rd_3[2], '.3f'), format(rd_30[2], '.3f')]
    # data[r"$\delta=3$"] = [format(rd_3[3], '.3f'), format(rd_30[3], '.3f')]
    # data[r"$\delta=4$"] = [format(rd_3[4], '.3f'), format(rd_30[4], '.3f')]
    # data[r"$\delta=5$"] = [format(rd_3[5], '.3f'), format(rd_30[5], '.3f')]
    # data[r"$\delta=6$"] = [format(rd_3[6], '.3f'), format(rd_30[6], '.3f')]
    # data[r"$\delta=7$"] = [format(rd_3[7], '.3f'), format(rd_30[7], '.3f')]
    # data[r"$\delta=8$"] = [format(rd_3[8], '.3f'), format(rd_30[8], '.3f')]


    # textabular = f"l|{'r'*len(headers)}"
    # texheader = " & " + " & ".join(headers) + "\\\\"
    # texdata = "\\hline\n"
    # for label in sorted(data):
    #     texdata += f"{label} & {' & '.join(map(str,data[label]))} \\\\\n"

    # print("\\begin{tabular}{"+textabular+"}")
    # print(texheader)
    # print(texdata,end="")
    # print("\\end{tabular}")
            
                
                
