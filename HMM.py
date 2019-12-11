#!/usr/bin/env python3


#Hidden Markov Model to find the most likely state path through E. coli DNA Sequence


"""
we want the probability of an emission of a base. if its the first in the triplet, we go 
by sum of the table. if its the second we calculate based on the first
if its the third we calculate based on the first again (or the second)
then we go back to the start. 
we should keep a counter of whether the base is the first, second or third in thesequenc 

if this is a non coding region it doesnt matter. 

so for example in coding state: chance of emitting A as first base is the sum prob of all the codons starting A
then chance of emitting C as the second base in triplet is different for each first base. 
We go back and check the previous one. if it is A for example, chance of getting C is chance of an AC_ codon out of all
codons starting with A. 
"""


"""
This program contains a script to find the most likely state path through
sequences of E.coli DNA, thus helping us identify coding and non coding regions. 

Arguments: 
    -out: file to which the intervals should be output to, line by line. 

"""

import argparse 
import numpy as np

def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s



def viterbi(obs, init_probs, noncod_emiss, coding_emiss_1, coding_emiss_2, coding_emiss_3, trans_probs):
    #we will first set up a matrix that two rows such that the first row represents the coding state C and the second row 
    #represents a non coding state N 


    matrix = np.zeros((2, len(obs)))
    traceback_matrix = np.zeros((2, len(obs)))

    #We then code the initial probabilties
    #we want to track the position in the codon - this will help during later calculations of emission probabilties. 
    codon_pos = 1
    matrix[0][0] = np.log(0.9) + coding_emiss_1[obs[0]]
    codon_pos = 2
    matrix[1][0] = np.log(0.1) + noncod_emiss[obs[0]]

    #Having set these probabilities, we want to iterate over the rest of the sequence   
    
    for j in range(1, len(obs)):
        base = obs[j]
        previous_coding = matrix[0][j-1]
        previous_non  = matrix[1][j-1]

        #This helps me update the score and update the traceback matrix at the same time to backtrack later

        #If we are at the first base in the codon
        if codon_pos == 1:
            matrix[0][j] = coding_emiss_1[base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos += 1
        
        #If we are at the second base in the codon 
        elif codon_pos == 2:
            matrix[0][j] = coding_emiss_2[obs[j-1]][base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos = 3

        #If we are at the third base in the codon 
        elif codon_pos == 3:
            matrix[0][j] = coding_emiss_3[obs[j-2]][base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos = 1

        #Set the correct value in the traceback matrix
        if (previous_coding + +trans_probs['C']['C'] > previous_non +trans_probs['N']['C']):
            traceback_matrix[0][j] = 1
        else:
            traceback_matrix[0][j] = 2


        #now code for tracking the non coding regions probabilities in our scoring matrix 

        matrix [1][j] = noncod_emiss[base] + max(previous_non+trans_probs['N']['N'], previous_coding +trans_probs['C']['N'])
        if (previous_non+trans_probs['N']['N'] > previous_coding +trans_probs['C']['N']):
            traceback_matrix[1][j] = 2
        else:
            traceback_matrix[1][j] =  1

    
    #Now we want to loop back through our traceback matrix that we have been building and identify 
    #the path with the maximum likelihood. 

    p = max(matrix[1][len(obs)-1], matrix[0][len(obs)-1])

    traceback_seq = ''

    if p == matrix[0][len(obs)-1]:
        traceback_seq += 'C'
        if traceback_matrix[0][len(obs)-1] == 1:
            traceback_seq += 'C'
        elif traceback_matrix[0][len(obs)-1] == 2:
            traceback_seq += 'N'
    else:
        traceback_seq += 'N'
        if traceback_matrix[1][len(obs)-1] == 1:
            traceback_seq += 'C'
        elif traceback_matrix[1][len(obs)-1] == 2:
            traceback_seq += 'N'

    k = len(obs) - 2
    while k!= 0:

        if traceback_seq[-1] == 'C':
            if traceback_matrix[0][k] == 1:
                traceback_seq += 'C'
            elif traceback_matrix[0][k] == 2:
                traceback_seq += 'N'

        elif traceback_seq[-1] == 'N':
            if traceback_matrix[1][k] == 1:
                traceback_seq += 'C'
            elif traceback_matrix[1][k] == 2:
                traceback_seq += 'N'
        k-=1


    traceback_seq = traceback_seq[::-1]

    l = []

    for element in traceback_seq:
        l.append(element)

    print(l)
    return l, p 





def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into coding and non coding regions')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    # parser.add_argument('-out', action="store", dest="out", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    # intervals_file = args.out

    #We read the FASTA file to obtain our required sequence
    obs_sequence = read_fasta(fasta_file)


    #we set up the initial probabilities
    initial_probabilities = {'n': np.log(0.9), 'c': np.log(0.1)}


    #for this, we simply use the relative frequencies
    non_coding_emission = {
        'A': np.log(0.2366),  'C': np.log(0.2530), 'G': np.log(0.2789), 'T': np.log(0.2315)
    }

    coding_emission_1 = {
        'A': np.log(0.245), 'C': np.log(0.245), 'G': np.log(0.36), 'T': np.log(0.147)
    }

    #for coding_emission_2, what we want is the probability of a second base appearing, GIVEN that we already know the first base
    #so for C, this will be the probabilty of all AC_ codons over the probability of all A _ _ codons (i.e x/0.245)
    coding_emission_2  = {
        'A': {'A': np.log(0.084/0.245), 'C': np.log(0.053/0.245), 'G': np.log(0.025/0.245), 'T': np.log(0.083/0.245)},
        'C': {'A': np.log(0.066/0.245), 'C': np.log(0.044/0.245), 'G': np.log(0.056/0.245), 'T': np.log(0.079/0.245) }, 
        'G': {'A': np.log(0.115/0.36), 'C': np.log(0.097/0.36), 'G': np.log(0.076/0.36), 'T': np.log(0.072/0.36)}, 
        'T': {'A': np.log(0.029/0.147), 'C': np.log(0.032/0.147), 'G': np.log(0.026/0.147), 'T': np.log(0.06/0.147)}
    }

    #for coding_emission_3, we want the probabilty of the 3rd codon given that we know the probabilty of the second
    coding_emission_3 = {
        'A': {'A': np.log(0.044/0.245), 'C': np.log(0.092/0.245), 'G': np.log(0.051/0.245), 'T': (0.058/0.245)}, 
        'C': {'A': np.log(0.027/0.245) , 'C': np.log(0.049/0.245), 'G': np.log(0.117/0.245), 'T': np.log(0.052/0.245)},
        'G': {'A': np.log(0.08/0.36), 'C': (0.094/0.36), 'G': np.log(0.091/0.36), 'T': np.log(0.095/0.36)},
        'T': {'A': np.log(0.017/0.147), 'C': np.log(0.048/0.147), 'G': np.log(0.034/0.147), 'T': np.log(0.048/0.147)}
    }

    transition_probabilities = {'C': {'C':0.9, 'N':0.1}, 'N': {'C':0.1, 'N':0.9}}


    viterbi(obs_sequence,initial_probabilities,non_coding_emission,coding_emission_1, coding_emission_2,coding_emission_3,transition_probabilities)


if __name__ == "__main__":
    main()