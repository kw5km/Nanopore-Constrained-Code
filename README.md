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
python State_Splitting/apply_splitting.py -d $\delta$ -max_run 2 -b 186 -k 6 -f_o <file_path> -f_m <file_path> + messages/ -f_s <file_path_score> + 6mer_means.txt -m <file_path_img>+Fall_2020_Grounds_scaled_gray.npy
```
* `-d`: $\delta$ constraint value
* `-max_run`: maximum number of consecutive homopolymers 
* `-b`: sequence length
* `-k`: $k$-mer size
* `-f_o`: file path out
* `-f_m`: file path out for binary messages
* `-f_s`: file path ONT $k$-mer [means file](https://github.com/nanoporetech/kmer_models) 
* `-m`: file path in for message to be encoded
In our pipeline, these encoded bases are then used to simulate nanopore sequencing current values with [DeepSimulator](https://github.com/liyu95/DeepSimulator).

### Basecalling with Viterbi
Signal measurement files can then be basecalle dusing our Viterbi Algorithm implementation, `Viterbi/multi_viterbi.py`
```
python Viterbi/multi_viterbi.py -maxrun $max_run -procc 'diff' -k 6 -dwell 8 -d $SLURM_ARRAY_TASK_ID -r ${rateArray[$SLURM_ARRAY_TASK_ID]} -f_i $path_var -f_o $path_var_out -f_b $file_out -f_s $score_file
```

### Decoding Basecalls

