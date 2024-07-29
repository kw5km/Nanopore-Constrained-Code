#Constrained Coding for Error Mitigation in Nanopore-Based DNA Data Storage

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
* Viterbi:
* Image: code for decoding sequence reads back into binary data
