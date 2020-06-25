# HMMGeneFinder
A Hidden Markov Model to identify coding regions and genes in sequences of E.Coli DNA

----

### Description
The HMM takes in an input E.Coli DNA Sequence and identifies Gene Coding regions of DNA, employing the Viterbi algorithm and other probabilistic dynamically programmed techniques. 

### Usage
The script takes in 3 arguments:

f - FASTA file with sequences in FASTA format. 
out - Output filename to display intervals of gene-coding regions. 

Example Usage:
```
python nwlinear.py -f <sequence filename> -out <output filename> 
```
A sample sequences file and output file are in this repository. Running the algorithm with an existing output file will overwrite the current contents of the output file. 
