from scripts.genome import find_melting_temperature
from scripts.duplex import Sequence

"""
    Miscellaneous functions
"""

def fformat(number):
    """
        Description: Formats negative and position numbers to fit uniformly
    """
    if(number >= 0): 
        return str(number) + " " + '\t'
    else:
        return str(number) + '\t'

def run_sequences(sequences, name, same_file):
    """
        Description: Find energy, enthalphy, entropy and temperature of each sequence in list

        Input: sequences - list of sequence, name - seudo-sequence, smae_file - sentinel to specify
            if output should be in the same file
    """
    # create output file with filename + "output.txt"
    if same_file == True:
        output = open("output/output.txt", 'a')
    else:
        output = open("output/" + name + "_output.txt", 'w')
    for x, sequence in enumerate(sequences):
        if not (x != 0 and sequences[x] == sequences[x-1]):
            prime_3 = sequence.split('/')[0]
            prime_5 = sequence.split('/')[1]
            s = find_melting_temperature(Sequence(prime_3, prime_5, .004))
            output.write('\n'+name+ '\n\n')
            output.write(str(sequence) + '\t' + fformat(s.energy) + fformat(s.enthalpy) + fformat(s.entropy) + fformat(s.temperature) + str(len(prime_3)))
            output.write('\n')

def add_tuples(t1, t2, size):
    """
        Definition: Add tuples of variable size
    """
    for x in range(len(t1)):
        t1[x] += t2[x]


