"""
    Compiler for a DNA duplex to calculate the dissociation and
    association melting temperaeure
"""

import random
import math
import sys
from subprocess import call
from duplex import Sequence

from energy import get_initiation_energy as i_energy
from energy import get_nearest_neighbor_energy as nn_energy
from energy import get_nearest_neighbor_mismatch_energy as nn_mm_energy
from energy import get_symmetry_g as get_sym

# from graphic import visualize_genome as visualize

"""
    Calculate melting temperature of dissociation based on nearest
    neighbor method
"""


def find_melting_temperature(s):
    s.set_complementary(True)
    t = 0
    t1 = 0
    sl = s.length
    for x in range(sl - 1):
        t = s.get_energy()
        # print("E: " + str(round((t - t1), 2)))
        # print(s.get_energy())
        # print(s.get_3_s(x, x + 2) + " " + s.get_5_s(x, x + 2))
        t1 = t
        # check for symmetry
        if(s.get_3(x) + s.get_5(x) != s.get_5(sl - 1 - x) + s.get_3(sl - 1 - x)):
            s.set_complementary(False)
        # find complementary base pair nearest neighbor energy
        if is_complement(s.get_5(x +1), s.get_3(x +1)) and is_complement(s.get_3(x), s.get_5(x)):
            # print(s.three_prime[x]+s.three_prime[x+1] +", "+ s.five_prime[x]+s.five_prime[x+1])
            s.add(nn_energy(s.get_3_s(x,x + 2), s.get_5_s(x, x + 2)))
        # look for base pair inverse mismatch value
        else:
        # print(s.three_prime[x]+s.three_prime[x+1] +", "+ s.five_prime[x]+s.five_prime[x+1])
            s.add(nn_mm_energy( s.get_3_s(x, x + 2), s.get_5_s(x, x + 2)))
            s.set_complementary(False)

# if s is symmetrical, add symmetry value
    s.add(get_sym(s.get_complementary()))

    # print("E: "+str(s.get_energy()))
    # calculate energy values of first pair ignoring nearest neighbor
    s.add(i_energy(s.get_initial(), s.get_terminal(), s.get_complementary()))

    # Account for Other factors that will affect temperature/gibbs energy

    """
        Calculate temperature from Gibbs Energy:
        [A], [B] - molarity of dissociated single strand
        [AB] - molarity of associated double strand
        If no additional nucleic acids are present, then [A], [B], and [AB]
        will be equal, and equal to half the initial concentration of
        double-stranded nucleic acid, [AB]initial. This gives an expression for
        the melting point of a nucleic acid duplex of
    """
    R = 1.987
    s.set_temperature((s.get_enthalpy()*1000) / (s.get_entropy() + (R * math.log(s.get_oligo_molarity() , math.exp(1)))))

    """
    temperature calculation with salt
    a = math.pow(3.92,-5)
    b = math.pow(-9.11,-6)
    c = math.pow(6.26,-5)
    d = math.pow(1.42,-5)
    e = math.pow(-4.82,-4)
    f = math.pow(5.25,-4)
    g = math.pow(8.31,-5)
    temperature_salt = (1/temperature) + a + b*math.log(Mg_mol, math.exp(1)) +/
        (f * (c+d*math.log(Mg_mol,math.exp(1))))+ ((1/(2*(length-1)))*(e +\
                (f*math.log(Mg_mol, math.exp(1))) +\
                (g*math.pow(math.log(Mg_mol, math.exp(1)), 2))))
    """

    return s


"""
   string from text file and creates DNA sequence for 5'->3' and 3'->5' strand
"""


def parse_genome(number):
    prime_5 = ""
    prime_3 = ""
    with open("dna_samples/dna_sample_" + str(number) + ".txt") as f:
        seq = f.read()
    seq = seq.split()
    even = True
    for x in range(len(seq)):
        if seq[x] != '/n' and seq[x] != ' ':
            if even:
                prime_3 += seq[x]
                even = False
            else:
                prime_5 += seq[x]
                even = True
    return Sequence(prime_3, prime_5, .004)


# Show help menu
def help():
    print("Commands:")
    print("")
    print("help: Show help menu")
    print("test: Run program on sequence and compare to expected results")
    print("show: Show sequence and expected result for that sequence")
    print("details: Explains experimental conditions and procedure")
    print("multiple: Enter menu for running and recording multiple sequences\n")


# runs program on dna sequence and compares it to expected value
def test(number):
    sequence = parse_genome(number)
    data = find_melting_temperature(sequence)
    print("Length: " + str(data.length))
    print("Gibb's free energy: " + str(data.energy))
    print("Entropy; " + str(data.entropy))
    print("Enthalpy: " + str(data.enthalpy))
    print("Temperature: " + str(data.temperature))
    print()
    # print correct result
    show_results(number)


# details about the program
def details():
    print("\nParameters:\n")
    print("pH: 7")
    print("oligonucleotide: 1E-4 for self-complementary sequences and 4E-4 for all others")
    print("\nMeasurements:\n")
    print("Gibb's free energy: kcal/mol")
    print("Enthalpy: kcal/mol")
    print("Entropy: eu")
    print("Temperature: K")
    print("\nGoal:\n")
    print("To find the temperature at which half of a DNA duplex")
    print("dissociates 1/2 of its base pairs from the Gibb's energy")
    print("contribution of each base pair bond.\n")


# prints out a dna sequence from a sample file
def show(number):
    print("Sequence\n")
    # print dna sequence
    with open("dna_samples/dna_sample_" + str(number) + ".txt") as f:
        print(f.read())
    # printing expected output for DNA sample
    show_results(number)


# show expected results on file
def show_results(num):
    print("Expected Output:")
    start = ((int(num) - 1) * 7) + 1
    end = start + 5
    s = str(start) + "," + str(end) + "p"
    call(["sed", "-n", s, "dna_samples/catalog.txt"])


# show all sample sequences
def show_all(n):
    for x in range(n):
        show(x + 1)


# test all samples
def test_all(n):
    for x in range(n):
        test(x + 1)
        print()


# test if bases are complements
def is_complement(b1, b2):
    if (b1 == 'T' and b2 == 'A') or (b1 == 'A' and b2 == 'T'):
        return True
    if (b1 == 'C' and b2 == 'G') or (b1 == 'G' and b2 == 'C'):
        return True
    return False


# Formatting negative and position numbers to fit uniformly
def fformat(number):
    if(number >= 0): 
        return str(number) + " " + '\t'
    else:
      return str(number) + '\t'

# interface to run multiple sequence
def multiple():
    print()
    print("Enter numbers to indicate what sequences you would like to test")
    print("Enter a number for the number of adjacent matches then M for a mismatch")
    print("This program will find the average melting temperature for G-C matches and A-T matches with")
    print("All possible mismatches. Putting R before the number will make the matches random")
    print()
    print("Example 1: 3M3 -- three matches , a mismatch and three more matches")
    print("Example 2: 4M2 -- four matches , a mismatch and two more matches")
    print("Example 3: 3 -- three matches")
    print("Example 4: R3MR3 -- three random matches , a mismatch and three more random matches")
    print("Enter 'quit', 'exit', or 'q' to exit")
    print()
    i = 1
    sequence = ""
    while sequence != "exit" and sequence != "quit" and sequence != "q":
        sequence = input(str(i) + ": ")
        if sanitize_m(sequence):
            s = seudo_sequence_converter(sequence)
            run_sequences(add_mismatches(s), sequence)
        i += 1

# sanitize multiple sequence input
def sanitize_m(input_m):
    for x, item in enumerate(input_m):  
        if item != "M" and item != "R" and not item.isdigit():
            return False
        elif item == "M":
            if x == 0 or x == len(input_m)-1:
                print("Mismatch cannot be on the end")
                return False
        elif item == "R":
            if x == len(input_m)-1 or not input_m[x].isdigit():
                return False
    return True

# creates and tests sequences
def seudo_sequence_converter(sequence):
    # [3'-high, 5'-high, 3'-low, 5'-low]
    strands = ["","","",""]
    random = False
    sequences = []
    for x, item in enumerate(sequence):	
        if item == "R":
            random = True
        elif item.isdigit() and random:
            add_stretch_rand(strands, item)
            random = False
        elif item.isdigit():
            add_stretch(strands, item)
        elif item == "M":
            strands[0] += "M"
            strands[1] += "M"
            strands[2] += "M"
            strands[3] += "M"
        else:
            print("seudo-sequence character not recognized")
    sequences.append(strands[0] + '/' + strands[1])
    sequences.append(strands[2] + '/' + strands[3])
    return sequences

# Add stretch of matches/mismatches to a strand
def add_stretch(strands,  number):
    basepairs = ["AT", "GC"]
    number = int(number)
    strands[0] += number * basepairs[1][0]
    strands[1] += number * basepairs[1][1]
    strands[2] += number * basepairs[0][0]
    strands[3] += number * basepairs[0][1]
    return strands

# add random stretch of base pairs
def add_stretch_rand(strands,  number):
    basepairs = ["AT", "GC"]
    for x in range(number):
        integer = random.randint(0,1)
        strands[0] += basepairs[integer][0]
        strands[1] += basepairs[integer][1]
        strands[2] += basepairs[integer][0]
        strands[3] += basepairs[integer][1]
    return strands

# find energy, enthalphy, entropy and temperature of sequences
def run_sequences(sequences, name):
    # create output file with filename + "output.txt"
    output = open(name + "_output.txt", 'w')
    for x, sequence in enumerate(sequences):
        prime_3 = sequence.split('/')[0]
        prime_5 = sequence.split('/')[1]
        s = find_melting_temperature(Sequence(prime_3, prime_5, .004))
        output.write(str(sequence) + '\t' + fformat(s.energy) + fformat(s.enthalpy) + fformat(s.entropy))
        output.write('\n')

# add tuples of variable size
def add_tuples(t1, t2, size):
    for x in range(len(t1)):
        t1[x] += t2[x]

# go through sequences and add possible mismatches is stipulated
# this will only work if there is one mismatch in the sequence
def add_mismatches(sequence):
    bases = ["AC", "AG", "AA", "GA", "GG", "GT", "CA", "CC", "CT", "TC", "TT", "TG"]
    new_sequence = []
    for x, item in enumerate(sequence):
        if "M" in item:
            matched = item.split("M")
            for x in range(len(bases)):
                new_sequence.append(matched[0]+ bases[x][0] + matched[1] + bases[x][1] + matched[2])
        else:
            new_sequence.append(item)
    return new_sequence

if __name__ == '__main__':


    if len(sys.argv) > 1:
        print("Program only accepts two total command line arguments")
        print("usage: python compiler.py")
        sys.exit(0)

    details()
    help()

    # number of samples
    with open("dna_samples/catalog.txt") as f:
        samples = int(sum(1 for _ in f) / 7)

    # Command prompt loop
    while True:
        command = input("[DNA_compiler]: ").strip().split()
        if len(command)  < 1:
            print("Invalid command: enter help for command menu")
        elif command[0] == "help" or command[0] == "h":
            help()
        elif command[0] == "multiple" or command[0] == "m":
            multiple()
        elif command[0] == "test" or command[0] == "t":
            if len(command) > 1:
                if command[1] == "all":
                    test_all(samples)
                elif not command[1].isdigit():
                    print("2nd argument has to be an integer")
                elif not int(command[1]) < 26 and not int(command[1]) > 0:
                    print("2nd argument has to be an integer between 0 and 26")
                else:
                    test(command[1])
            else:
                file_name = input("select a number from 1 and 25: ").strip()
                if file_name.isdigit():
                    test(file_name)
                elif file_name == "all":
                    test_all(samples)
                else:
                    print("Integer between 0 and 26 required")
                    continue
        elif command[0] == "show" or command[0] == "s":
            if len(command) > 1:
                if command[1] == "all":
                    show_all(samples)
                elif not command[1].isdigit():
                    print("2nd argument has to be an integer")
                elif int(command[1]) > 25 or int(command[1]) < 1:
                    print("2nd argument has to be an integer between 0 and 26")
                else:
                    show(command[1])
            else:
                file_name = input("select a number from 1 and 25: ").strip()
                if file_name.isdigit():
                    show(file_name)
                elif file_name == "all":
                    show_all(samples)
                else:
                    print("Integer between 0 and 26 required")
                    continue
        elif command[0] == "exit" or command[0] == "e" or command[0] == "quit" or command[0] == "q":
            sys.exit()
        elif command[0] == "details" or command[0] == "d":
            details()
        else:
            print("Invalid command: enter help for command menu")
