# -*- coding: utf-8 -*-
"""

@author: Kallie Whritenour
"""
import numpy as np
from PIL import Image


def bitstring_to_bytes(s):
    bits = np.array([int(b) for b in s])
    bytes = np.packbits(bits)
    return bytes


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-delta", type=int, default=0, help="Delta constraints to look at")
    parser.add_argument("-f_o", type=str, default=0, help="File path out")
    parser.add_argument("-f_i", type=str, default=0, help="File path input image")
    parser.add_argument("-f_i_gray", type=str, default=0, help="File path input image gray")
    args = parser.parse_args()
    d, f_o, f_i, f_i_gray = args.delta, args.f_o, args.f_i, args.f_i_gray


    img = Image.open(f_i).convert('L')
    img.save(f_i_gray)

    im = Image.open(f_i_gray).tobytes()
    bytes_im = np.frombuffer(im, dtype=np.uint8)
    bits_im = np.unpackbits(bytes_im)

    np.save(f_o + 'Fall_2020_Grounds_scaled_gray', bits_im)