# Constrained Coding for Error Mitigation in Nanopore-Based DNA Data Storage

This repository includes code and data for the implementation of a de Bruijn based constrained code tailord to error correction for nanopore seqeuncing. This error correcting code is intended for use in a DNA data storage pipeline when nanopore sequencing is also used.

## Compatibility
This code was run on `Ubuntu 22.04.4` with `conda 4.12.0`

## Contents Description
* `requirements.txt`:
  run the following in order to create a conda env with the library requirements for this repository.
 ```
 conda create --name <env_name> --file requirements.txt
 ```
* deBruijn: utils and functions related to building a de Bruijn graph and applying it to encoding problems.
* StateSplittng: implementation of State Splitting algorithm to enforce a minimum out degree on a given graph
* Viterbi: threaded viterbi for handling multiple sequences at once
* Image: code for decoding sequence reads back into binary data
  
## Running Instructions
### Producing Sequences
'StateSplitting/apply_splitting.py' applies the state splitting algorithm to de Bruijn graphs of different $\delta$ constraints, then uses your `-m` filepath as the input data for encoding to DNA bases. 
```
python State_Splitting/apply_splitting.py -d delta -max_run 2 -b 186 -k 6 -f_o <file_path> -f_m <file_path> + messages/ -f_s <file_path_score> + 6mer_means.txt -m <file_path_img>+Fall_2020_Grounds_scaled_gray.npy
```
* `-d`: $\delta$ constraint value
* `-max_run`: maximum number of consecutive homopolymers 
* `-b`: sequence length
* `-k`: $k$-mer size
* `-f_o`: file path out
* `-f_m`: file path out for binary messages
* `-f_s`: file path ONT $k$-mer [means file](https://github.com/nanoporetech/kmer_models) 
* `-m`: file path in for message to be encoded

To produce Reed Solomon redundancy bits, use `RS_singleCodeword.m`
In our pipeline, these encoded bases are then used to simulate nanopore sequencing current values with [DeepSimulator](https://github.com/liyu95/DeepSimulator).

### Basecalling with Viterbi
Signal measurement files can then be basecalle dusing our Viterbi Algorithm implementation, `Viterbi/multi_viterbi.py`
```
python Viterbi/multi_viterbi.py -maxrun $max_run -procc diff -k 6 -dwell 8 -d delta -r <rate> -f_i <file_path_in> -f_o <file_path_out> -f_b <file_path_out_basecalls> -f_s <file_path_score> + 6mer_means.txt
```
* `-maxrun`: maximum number of consecutive homopolymers
* `-procc`: ISI mitigation preprocessing step. `diff` is used exclusively in this work.
* `-k`: $k$-mer size
* `-d`: $\delta$ constraint value
* `-f_o`: file path prefix out
* `-f_b`: file path out for basecalls
* `-f_s`: file path ONT $k$-mer [means file](https://github.com/nanoporetech/kmer_models) 

### Decoding Basecalls
Basecalled sequences produced using our Viterbi implementation can be decoded back to binary using the following pipeline.

```
python Image/vit2bits.py -delta delta -f_o <file_path_out> -f_p <file_path_padding> -f_c <file_path_check>
python Image/vit2consensus.py -delta delta -f_o <file_path_out> -f_p <file_path_padding>
python Image/consensus2bits.py -delta delta -f_o <file_path_out> -f_p <file_path_padding> -f_c <file_path_check>
python Image/bits2img.py -f_o <file_path_out> -f_c <file_path_check>
```
* `-delta`: $\delta$ constraint value
* `-f_o`: file path out
* `-f_p`: file path location of saved padding info
* `-f_c`: file path for teh orignal message outouts to check against decoded outputs
* 
  After decoding, RS decoding can be applied to further correct any residual errors or missing sequences using `RS_singleCodeword_decode.m`. Results can be viewed using `paper_PLOT.m`.
