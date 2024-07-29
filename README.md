#Constrained Coding for Error Mitigation in Nanopore-Based DNA Data Storage
## Contents Description
* requirements.txt: run in order to create a conda env with the library requirements for this repository.
 ```
  conda create --name <env_name> --file requirements.txt
 ```
* deBruijn: utils and funcctions related to building a de Bruijn graph and applying it to encoding problems.
* StateSplittng: implementation of State Splitting algorithm to enforce a minimum out degree on a given graph
* Viterbi:
* Image: code for decoding sequence reads back into binary data
