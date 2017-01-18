"""
    Compiler for a DNA duplex to calculate the dissociation and
    association melting temperature with one or two base pair
    mismatches
"""

import random
import numpy
import math
import sys
from duplex import Sequence
from termcolor import colored
from subprocess import call

from energy import get_initiation_energy as i_energy
from energy import get_nearest_neighbor_energy as nn_energy
from energy import get_nearest_neighbor_mismatch_energy as nn_mm_energy
from energy import get_symmetry_g as get_sym

"""
    Global variables
"""
bases_pairs = ["AC", "AG", "AA", "GA", "GG", "GT", "CA", "CC", "CT", "TC", "TT", "TG"]
bases = ['A', 'C', 'G', 'T']


def find_melting_temperature(s):
    """
        Description: Calculate melting temperature of dissociation based on nearest
        neighbor method

        Returns: Sequence Objects with energy, enthalpy, entropy, and temperature values

        Input:Sequence object with strand information
    """
    s.set_complementary(True)
    sl = s.length
    for x in range(sl - 1):
        # check for symmetry
        if(s.get_3(x) + s.get_5(x) != s.get_5(sl - 1 - x) + s.get_3(sl - 1 - x)):
            s.set_complementary(False)
        # find complementary base pair nearest neighbor energy
        if is_complement(s.get_5(x +1), s.get_3(x +1)) and is_complement(s.get_3(x), s.get_5(x)):
            s.add(nn_energy(s.get_3_s(x,x + 2), s.get_5_s(x, x + 2)))
        # look for base pair inverse mismatch value
        else:
            s.add(nn_mm_energy( s.get_3_s(x, x + 2), s.get_5_s(x, x + 2)))
            s.set_complementary(False)

    # if s is symmetrical, add symmetry value
    s.add(get_sym(s.get_complementary()))

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

    return s


def parse_genome(number):
    """
       Description: Creates DNA sequence for 5'->3' and 3'->5' strand

       Returns: Sequence object with 3' and 5' sequences

       Input: Number of DNA sample file to read from
    """
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

    return Sequence(prime_3, prime_5, .0004)


def help():
    """
        Description: Prints help menu
    """
    print("Commands:")
    print("")
    print("help: Show help menu")
    print("test: Run program on sequence and compare to expected results")
    print("show: Show sequence and expected result for that sequence")
    print("details: Explains experimental conditions and procedure")
    print("multiple: Enter menu for running and recording multiple sequences")
    print("difference: Enters menu for dining average difference of mismatch substitution")
    print("manual: Enters mode where user can manually enter sequence")
    print()

def title():
    """
        Description: Prints introduction title
    """
    width = 125
    red = 'red'
    green = 'green'
    print('\n')
    print(colored('#' * 127, green))
    print('\n')
    print(colored(" /=\              ______  __   __     ___                    /=\  ", red).center(width, ' '))
    print(colored("/===\             |  _  \ | \ | |    / _ \                  /===\ ", red).center(width, ' '))
    print(colored("|===|             | | | | |  \| |   / /_\ \                 |===| ", red).center(width, ' '))
    print(colored("\===/             | |_| | | |\  |  / /___\ \                \===/ ", red).center(width, ' '))
    print(colored(" \=/              |_____/ |_| \_| /_/     \_\                \=/  ", red).center(width, ' '))
    print(colored("  X     ___   ___   _      _  ____   _   _     ___   ___      X   ", red).center(width, ' '))
    print(colored(" /=\   /  _\ /   \ | \    / | |   | | | | |   | __| | . |    /=\  ", red).center(width, ' '))
    print(colored("/===\  | |   | | | |  \  /  | | | | | | | |   | |_  | __|   /===\  ", red).center(width, ' '))
    print(colored("|===|  | |   | | | |   \/   | | __| | | | |   |  _| |  \    |===| ", red).center(width, ' '))
    print(colored("\===/  | |_  | | | | |\__/| | | |   | | | |_  | |_  | | \   \===/ ", red).center(width, ' '))
    print(colored(" \=/   \___/ \___/ |_|    |_| |_|   |_| |___| |___| |_||_\   \=/  ", red).center(width, ' '))
    print('\n')
    print(colored('#' * 127, green))


def test(number):
    """
        Description: runs program on dna sequence and compares it to expected value

        Input: Number of DNA sample file to test

        Output: Sequence object credentials (energy, enthalpy, etc.)
    """
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


def details():
    """
        Description: Prints details about the program
    """
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


def show(number):
    """
        Description: Prints DNA sequence from sample file

        Input: Number of DNA sample file to read from

        Output: DNA Sequence
    """
    print("Sequence\n")
    # print dna sequence
    with open("dna_samples/dna_sample_" + str(number) + ".txt") as f:
        print(f.read())
    # printing expected output for DNA sample
    show_results(number)


def show_results(num):
    """
        Description: Prints out a dna sequence from a sample file

        Input: Number of DNA sample file to read from

        Output: Expected values for the sample DNA sequence
    """
    print("Expected Output:")
    start = ((int(num) - 1) * 7) + 1
    end = start + 5
    s = str(start) + "," + str(end) + "p"
    call(["sed", "-n", s, "dna_samples/catalog.txt"])


def show_all(n):
    """
        Description: Shows all sample DNA file sequences and expected values
    """
    for x in range(n):
        show(x + 1)


def test_all(n):
    """
        Description: Tests all sample DNA file sequences and expected values
    """
    for x in range(n):
        test(x + 1)
        print()


def is_complement(b1, b2):
    """ 
        Description: Tests if bases are complements
    """
    if (b1 == 'T' and b2 == 'A') or (b1 == 'A' and b2 == 'T'):
        return True
    if (b1 == 'C' and b2 == 'G') or (b1 == 'G' and b2 == 'C'):
        return True
    return False


def fformat(number):
    """
        Description: Formats negative and position numbers to fit uniformly
    """
    if(number >= 0): 
        return str(number) + " " + '\t'
    else:
        return str(number) + '\t'

def multiple():
    """
        Description: Interface to run multiple sequence

        Output(file): Sequence credentials
    """
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
    print()
    print("Enter 'quit', 'exit', or 'q' to exit")
    print("Would you like to store all reslts in the same file (or separate files)?")
    same_file = False
    choice = input("(y or n): ")
    print()
    if choice == "y" or choice == "yes":
        open("output.txt",'w')
        same_file = True

    i = 1
    sequence = input(str(i) + ": ")
    while sequence != "exit" and sequence != "quit" and sequence != "q":
        if sanitize_m(sequence):
            s = seudo_sequence_converter(sequence)
            run_sequences(add_mismatches(s), sequence, same_file)
        i += 1
        sequence = input(str(i) + ": ")

def sanitize_m(input_m):
    """
        Description: Sanitize multiple sequence input

        Returns: True if the input is valid, False otherwise

        Input: Suedo-sequence for multiple sequence generation
    """
    for x, item in enumerate(input_m):  
        if item != "M" and item != "R" and not item.isdigit():
            print("item: "+ str(item))
            print("Input must be either a n 'M', 'R', or a number")
            return False
        elif item == "M":
            if x == 0 or x == len(input_m)-1:
                print("Mismatch cannot be on the end")
                return False
        elif item == "R":
            if x == len(input_m)-1 or not input_m[x+1].isdigit():
                print("R cannot be the last character and must be followed by a digit")
                return False
    return True

def seudo_sequence_converter(sequence):
    """
        Descripton: Creates sequences from seudo-sequences and adds them to a list

        Returns: List of generated sequences

        Input: Seudo-sequence, ex: R3MR3 --three random basepairs, mismatch, three random base pairs
    """
    # [3'-high, 5'-high, 3'-low, 5'-low]
    strands = ["","","",""]
    random = False
    sequences = []
    x = 0
    while x < len(sequence):	
        if sequence[x] == "R":
            random = True
        elif sequence[x].isdigit() and random:
            temp = sequence[x]
            i = 1
            for i in range(1, len(sequence)-x+2):
                if not sequence[x:x+i].isdigit():
                    break
                else:
                    temp = sequence[x:x+i]

            add_stretch_rand(strands, temp)
            x+= (i-2)
            random = False
        elif sequence[x].isdigit():
            temp = sequence[x]
            i = 1
            for i in range(1, len(sequence)-x+2):
                if not sequence[x:x+i].isdigit():
                    break
                else:
                    temp = sequence[x:x+i]

            add_stretch(strands, temp)
            x+= (i-2)
        elif sequence[x] == "M":
            strands[0] += "M"
            strands[1] += "M"
            strands[2] += "M"
            strands[3] += "M"
        else:
            print("seudo-sequence character not recognized")
        x+=1
    sequences.append(strands[0] + '/' + strands[1])
    sequences.append(strands[2] + '/' + strands[3])
    return sequences

def add_stretch(strands,  number):
    """
        Description Add stretch of matches/mismatches to a strand

        Returns: Strands with base pair matches added

        Input: List - with strands to appended to, Number - number of additions of a base pair 
    """
    basepairs = ["AT", "GC"]
    number = int(number)
    strands[0] += number * basepairs[1][0]
    strands[1] += number * basepairs[1][1]
    strands[2] += number * basepairs[0][0]
    strands[3] += number * basepairs[0][1]
    return strands

def add_stretch_rand(strands,  number):
    """
        Description Add stretch of random matches/mismatches to a strand

        Returns: Strands with random base pair matches added

        Input: strands - with strands to appended to, number - number of additions of a base pair 
    """
    basepairs = ["AT", "GC", "TA", "CG"]
    for x in range(int(number)):
        integer = random.randint(0,3)
        strands[0] += basepairs[integer][0]
        strands[1] += basepairs[integer][1]
        strands[2] += basepairs[integer][0]
        strands[3] += basepairs[integer][1]
    return strands

def add_stretch_rand_double(number):
    """
        Description Add stretch of random matches/mismatches to a strand

        Returns: Strands with random base pair matches added to two strands

        Input: strands - two strands to appended to, number - number of additions of a base pair 
    """

    three = ""
    five = ""
    basepairs = ["AT", "GC", "TA", "CG"]

    for x in range(int(number)):
        integer = random.randint(0,3)
        three += basepairs[integer][0]
        five += basepairs[integer][1]

    return three + '/' + five



def run_sequences(sequences, name, same_file):
    """
        Description: Find energy, enthalphy, entropy and temperature of each sequence in list

        Input: sequences - list of sequence, name - seudo-sequence, smae_file - sentinel to specify
            if output should be in the same file
    """
    # create output file with filename + "output.txt"
    if same_file == True:
        output = open("output.txt", 'a')
    else:
        output = open(name + "_output.txt", 'w')
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

def add_mismatches(sequence):
    """
        Description: Goes through sequences and add possible mismatches as stipulated

        Rule: This will only work if there is one mismatch in the sequence
    """
    new_sequence = []
    for x, item in enumerate(sequence):
        if "M" in item:
            matched = item.split("M")
            for x in range(len(bases)):
                new_sequence.append(matched[0]+ base_pairs[x][0] + matched[1] + base_pairs[x][1] + matched[2])
        else:
            new_sequence.append(item)
    return new_sequence


def difference():
    """
        Description: Calculates average difference in melting temperature after a mismatch substitution

        Output(database): Average difference of sequence melting temperature and single mismatch sequence
            sequence melting temperature
    """

    print("Calculate the average difference between a DNA sequence of length k and that sequence with one")
    print("mismatch substituted in")
    print("Enter d or default for the default lengths ([5, 10, 20, 50, 100, 130]) or enter manual (or m)")
    print("to enter the numbers manually")

    """
    lengths = []
    data = input(("d or m: "))
    while data != "d" and data != "default" and data != "manual" and data != "m":
        print("Command not recognized")
        data = input(("d or m: "))
    if data == "m" or data == "manual":
        i = 1
        while data != "" and data != "quit" and data != "q":
            data = intput(str(i) + ": ")
            lengths.append(data)
            i += 1
    else:
        lengths = [5, 10, 20, 50, 100, 130]
    """

    lengths = [5, 10, 20, 50, 100, 130]

    for x in lengths:
        val = calculate_differences(x)
        print("Length: " + str(x))
        print("Mean : " + str(val[0]) + " Kelvin")
        print("Median : " + str(val[1]) + " Kelvin")
        print("90% Confidence : " + val[2] + " Kelvin")

def calculate_differences(length):
    """
        Description: Calculate 100 matched sequence of size length, from those 100 sequences replace one match 
            with a mismatch (10 times for each sample) and calculate the difference in temperature after the 
            mismatch substitute. from these values, find a mean, median and 90% confidence

        Returns: tuple with mean, median, 90% confidence
    """

    sequences = []
    mismatched_sequences = []
    
    # add Sequence strings
    for x in range(100):
        sequences.append(add_stretch_rand_double(length))
        mismatched_sequences.append([])
        for i in range(10):
            mismatched_sequences[x].append(substitute_mismatch(sequences[x]))

    # replace strings with temperatures
    for x in range(100):
        s1 = find_melting_temperature(Sequence(sequences[x].split('/')[0], sequences[x].split('/')[1], .004))
        for i in range(10):
            s2 = find_melting_temperature(Sequence(mismatched_sequences[x][i].split('/')[0], mismatched_sequences[x][i].split('/')[1], .004))
            mismatched_sequences[x][i] = math.fabs(s1.temperature - s2.temperature)

    mismatched_sequences = list(sum(mismatched_sequences, []))
    mismatched_sequences.sort()
    mean = sum(mismatched_sequences)/1000
    median = (mismatched_sequences[500] + mismatched_sequences[499])/2
    ninety_percent_conf = str(mean) + " ± " + str(get_confidence_level(90, mismatched_sequences))
    return [mean, median, ninety_percent_conf]

def get_confidence_level(percentage, mismatched_sequences):
    """
            Descripton: find confidence level from data -- formula( x̅ ± Za/2 * σ/√(n)) 
    """
    std_dev = numpy.std(mismatched_sequences)
    confidence_level = (lookup_t_distribution(percentage/200)) * (std_dev/math.sqrt(1000))
    return confidence_level

def lookup_t_distribution(value):

    tenths = [0.0,0.00  ,0.01  ,0.02  ,0.03  ,0.04  ,0.05  ,0.06  ,0.07  ,0.08  ,0.0900]

    ones = [0.0,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7 , 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 
    2.4, 2.5, 2.6, 2.7, 2.8,2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7 ,3.8, 3.9, 4.0] 
    
    Z_list = [
        0.000 ,0.0040,0.0080,0.0120,0.0160,0.0199,0.0239,0.0279,0.0319,0.0359,
        0.0398,0.0438,0.0478,0.0517,0.0557,0.0596,0.0636,0.0675,0.0714,0.0753,
        0.0793,0.0832,0.0871,0.0910,0.0948,0.0987,0.1026,0.1064,0.1103,0.1141,
        0.1179,0.1217,0.1255,0.1293,0.1331,0.1368,0.1406,0.1443,0.1480,0.1517,
        0.1554,0.1591,0.1628,0.1664,0.1700,0.1736,0.1772,0.1808,0.1844,0.1879,
        0.1915,0.1950,0.1985,0.2019,0.2054,0.2088,0.2123,0.2157,0.2190,0.2224,
        0.2257,0.2291,0.2324,0.2357,0.2389,0.2422,0.2454,0.2486,0.2517,0.2549,
        0.2580,0.2611,0.2642,0.2673,0.2704,0.2734,0.2764,0.2794,0.2823,0.2852,
        0.2881,0.2910,0.2939,0.2967,0.2995,0.3023,0.3051,0.3078,0.3106,0.3133,
        0.3159,0.3186,0.3212,0.3238,0.3264,0.3289,0.3315,0.3340,0.3365,0.3389,
        0.3413,0.3438,0.3461,0.3485,0.3508,0.3531,0.3554,0.3577,0.3599,0.3621,
        0.3643,0.3665,0.3686,0.3708,0.3729,0.3749,0.3770,0.3790,0.3810,0.3830,
        0.3849,0.3869,0.3888,0.3907,0.3925,0.3944,0.3962,0.3980,0.3997,0.4015,
        0.4032,0.4049,0.4066,0.4082,0.4099,0.4115,0.4131,0.4147,0.4162,0.4177,
        0.4192,0.4207,0.4222,0.4236,0.4251,0.4265,0.4279,0.4292,0.4306,0.4319,
        0.4332,0.4345,0.4357,0.4370,0.4382,0.4394,0.4406,0.4418,0.4429,0.4441,
        0.4452,0.4463,0.4474,0.4484,0.4495,0.4505,0.4515,0.4525,0.4535,0.4545,
        0.4554,0.4564,0.4573,0.4582,0.4591,0.4599,0.4608,0.4616,0.4625,0.4633,
        0.4641,0.4649,0.4656,0.4664,0.4671,0.4678,0.4686,0.4693,0.4699,0.4706,
        0.4713,0.4719,0.4726,0.4732,0.4738,0.4744,0.4750,0.4756,0.4761,0.4767,
        0.4772,0.4778,0.4783,0.4788,0.4793,0.4798,0.4803,0.4808,0.4812,0.4817,
        0.4821,0.4826,0.4830,0.4834,0.4838,0.4842,0.4846,0.4850,0.4854,0.4857,
        0.4861,0.4864,0.4868,0.4871,0.4875,0.4878,0.4881,0.4884,0.4887,0.4890,
        0.4893,0.4896,0.4898,0.4901,0.4904,0.4906,0.4909,0.4911,0.4913,0.4916,
        0.4918,0.4920,0.4922,0.4925,0.4927,0.4929,0.4931,0.4932,0.4934,0.4936,
        0.4938,0.4940,0.4941,0.4943,0.4945,0.4946,0.4948,0.4949,0.4951,0.4952,
        0.4953,0.4955,0.4956,0.4957,0.4959,0.4960,0.4961,0.4962,0.4963,0.4964,
        0.4965,0.4966,0.4967,0.4968,0.4969,0.4970,0.4971,0.4972,0.4973,0.4974,
        0.4974,0.4975,0.4976,0.4977,0.4977,0.4978,0.4979,0.4979,0.4980,0.4981,
        0.4981,0.4982,0.4982,0.4983,0.4984,0.4984,0.4985,0.4985,0.4986,0.4986,
        0.4987,0.4987,0.4987,0.4988,0.4988,0.4989,0.4989,0.4989,0.4990,0.4990,
        0.4990,0.4991,0.4991,0.4991,0.4992,0.4992,0.4992,0.4992,0.4993,0.4993,
        0.4993,0.4993,0.4994,0.4994,0.4994,0.4994,0.4994,0.4995,0.4995,0.4995,
        0.4995,0.4995,0.4995,0.4996,0.4996,0.4996,0.4996,0.4996,0.4996,0.4997,
        0.4997,0.4997,0.4997,0.4997,0.4997,0.4997,0.4997,0.4997,0.4997,0.4998,
        0.4998,0.4998,0.4998,0.4998,0.4998,0.4998,0.4998,0.4998,0.4998,0.4998,
        0.4998,0.4998,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,
        0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,
        0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999,0.4999
    ]
    
    min_index = min(range(len(Z_list)), key=lambda i: abs(Z_list[i]-value))
    col =  int(min_index % 10)
    row = int(min_index / 10)
    print(Z_list[min_index])
    return tenths[col] + ones[row]

def substitute_mismatch(strand):
    sub_index = random.randint(1, len(strand.split('/')[0]) - 2)

    # keep getting random base until it is different from present one
    base = bases[random.randint(0,3)]
    while strand[sub_index] == base:
        base = bases[random.randint(0,3)]

    return strand.split('/')[0][:sub_index] + base + strand.split('/')[0][sub_index+1:] + '/' + strand.split('/')[1]

def manual():
    """
        Description: Manually enter DNA sequence (with / separating single strands)
            Ex: >> ACT/TGA
        Returns: energy, enthalpy, entropy and temperature values for duplex
    """
    print("Enter DNA sequence separated by '/', enter nothing, 'quit' or q to finish")
    sequence = input(">> ")
    while sequence != "" and sequence != "quit" and sequence != "q":
        if sanitize_sequence(sequence):
            seq = find_melting_temperature(Sequence(sequence.split('/')[0], sequence.split('/')[1], .0004))
            print("Energy: " + str(seq.energy))
            print("Enthalpy: " + str(seq.enthalpy))
            print("Entropy: " + str(seq.entropy))
            print("Temperature: " + str(seq.temperature))
        else:
            print("Invalid Input")
        sequence = input(">> ")

def sanitize_sequence(iput):
    if len(iput.split('/')) != 2:
        return False
    if len(iput.split('/')[0]) != len(iput.split('/')[1]):
        return False
    for x in iput:
        if x != 'C' and x != 'G' and x != 'A' and x != 'T' and x != '/':
            return False
    return True

def store():
    """
        Description: use SQL queries to store in database
    """
    pass

if __name__ == '__main__':

    if len(sys.argv) > 1:
        print("Program only accepts two total command line arguments")
        print("usage: python compiler.py")
        sys.exit(0)

    call(['clear'])
    title()
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
        elif command[0] == "difference" or command[0] == "d":
            difference()
        elif command[0] == "manual" or command[0] == "man":
            manual()
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
