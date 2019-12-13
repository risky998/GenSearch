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


'''
A helper function to read and parse a FASTA file to give us appropriate data to work with. 

Arguments: 
    filename: The name of the FASTA file that is to be parsed. 
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s

'''
The viterbi algorithm finds and returns the maximum likelihood path through a given sequence

Arguments: 
    obs: The sequence to be analyzed 
    init_probs: A dictionary of initial probabilties that represents the chance of the HMM starting in either of our two states
    noncod_emiss: A dictionary that represents the probability of a base emission in the intergenic state of the model
    coding_emiss_1: A dictionary that represents the probabilty of a base emission given that the base is the first in the sequence
    coding_emiss_2: A dictionary that represents the probabilty of a base emission given that the base is the second in the sequence
    coding_emiss_2: A dictionary that represents the probabilty of a base emission given that the base is the third in the sequence
    trans_prob: A dictionary that represents the transition probabilties between the coding and non-coding states

'''

def viterbi(obs, init_probs, noncod_emiss, coding_emiss_1, coding_emiss_2, coding_emiss_3, trans_probs):

    #Create a dynamic programming score and traceback matrix using numpy
    matrix = np.zeros((2, len(obs)))
    traceback_matrix = np.zeros((2, len(obs)))

    #Add the initial probabilities for each state to the score matrix 
    codon_pos = 1
    matrix[0][0] = np.log(0.9) + coding_emiss_1[obs[0]]
    codon_pos = 2
    matrix[1][0] = np.log(0.1) + noncod_emiss[obs[0]]

    #Iterate over the rest of the observed sequence    
    for j in range(1, len(obs)):
        base = obs[j]
        previous_coding = matrix[0][j-1]
        previous_non  = matrix[1][j-1]

        #For the first base in the codon 
        if codon_pos == 1:
            matrix[0][j] = coding_emiss_1[base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos += 1
        
        #For the second base in the codon
        elif codon_pos == 2:
            matrix[0][j] = coding_emiss_2[obs[j-1]][base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos = 3

        #For the third base in the codon
        elif codon_pos == 3:
            matrix[0][j] = coding_emiss_3[obs[j-2]][base] + max(previous_coding+trans_probs['C']['C'], previous_non +trans_probs['N']['C'])
            codon_pos = 1

        #Set the correct value in the traceback matrix
        if (previous_coding + +trans_probs['C']['C'] > previous_non +trans_probs['N']['C']):
            traceback_matrix[0][j] = 1
        else:
            traceback_matrix[0][j] = 2

        #Updating score and traceback matrix for intergenic state
        matrix [1][j] = noncod_emiss[base] + max(previous_non+trans_probs['N']['N'], previous_coding +trans_probs['C']['N'])
        if (previous_non+trans_probs['N']['N'] > previous_coding +trans_probs['C']['N']):
            traceback_matrix[1][j] = 2
        else:
            traceback_matrix[1][j] =  1

    
    #Use traceback matrix to find maximum likelihood path
    p = max(matrix[1][len(obs)-1], matrix[0][len(obs)-1])

    #initialize an empty string for us to use
    traceback_seq = ''

    #Add the first traceback to the string based on final observation in score matrix
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

    #Loop back through our observations and find the maximum probability in our traceback matrix
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

    #Reverse the traceback sequence to get the correct ordering
    traceback_seq = traceback_seq[::-1]
    l = []

    #Create a list that is easier to work with. 
    for element in traceback_seq:
        l.append(element)

    print(l)
    return l, p 

''' 
This function takes in the viterbi output and return the intervals that are coding regions 
Arguments:
	sequence: list of hidden states
Returns:
	intervals: list of tuples (i, j), 1 <= i <= j <= len(sequence), that
                   describe coding regions in the input list of hidden states.
'''
def find_intervals(sequence):

    viterbi = ''.join(sequence)

    list = []
    current_start =  viterbi.find('C')
    if current_start == -1:
        return list

    current_end = viterbi.find('N', current_start)
    if current_end == -1:
        current_end = len(viterbi)-1

    list.append((current_start+1, current_end))


    while current_start < len(viterbi) - 1:
        current_start = viterbi.find('C', current_end)
        if current_start == -1:
            return list
        current_end = viterbi.find('N', current_start)
        if current_end == -1:
            current_end = len(sequence)
        list.append((current_start+1, current_end))

    return list

def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into coding and non coding regions')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-out', action="store", dest="out", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    intervals_file = args.out

    obs_sequence = read_fasta(fasta_file)


    #we set up the initial probabilities according the expected value gene percentage in E.Coli 
    initial_probabilities = {'n': np.log(0.9), 'c': np.log(0.1)}

    #For non coding region emissions, we use relative probabilities of the 4 bases
    non_coding_emission = {
        'A': np.log(0.2366),  'C': np.log(0.2530), 'G': np.log(0.2789), 'T': np.log(0.2315)
    }

    #For coding emission on the first base in the codon, we calculate the probabilitities
    #of emitting any codon beginning with that base at each stage in the HMM
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

    #for coding_emission_3, we want the probabilty of the 3rd codon given that we know the probabilty of the first
    coding_emission_3 = {
        'A': {'A': np.log(0.044/0.245), 'C': np.log(0.092/0.245), 'G': np.log(0.051/0.245), 'T': (0.058/0.245)}, 
        'C': {'A': np.log(0.027/0.245) , 'C': np.log(0.049/0.245), 'G': np.log(0.117/0.245), 'T': np.log(0.052/0.245)},
        'G': {'A': np.log(0.08/0.36), 'C': (0.094/0.36), 'G': np.log(0.091/0.36), 'T': np.log(0.095/0.36)},
        'T': {'A': np.log(0.017/0.147), 'C': np.log(0.048/0.147), 'G': np.log(0.034/0.147), 'T': np.log(0.048/0.147)}
    }

    #Transition probabilities from non-coding region to coding region is based on the percentage of E.coli genome that codes for DNA 
    #Transition probabilities from coding to stop codon to non-coding region is based on 1/length of gene
    #Based on the study the average gene is 42 codons long
    transition_probabilities = {'C': {'C':41/42, 'N':1/42}, 'N': {'C':0.1, 'N':0.9}}

    #we call our viterbi sequence
    sequence, p = viterbi(obs_sequence,initial_probabilities,non_coding_emission,coding_emission_1, coding_emission_2,coding_emission_3,transition_probabilities)
    print(sequence)

    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("(%d,%d)" % (start, end)) for (start, end) in intervals]))
        f.write("\n")
    print("Viterbi probability: {:.2f}".format(p))


if __name__ == "__main__":
    main()