"""
    Compiler for a DNA duplex to calculate the dissociation and
    association melting temperature
"""

import math
import sys
from subprocess import call

from energy import get_end_energy as e_energy
from energy import get_initiation_energy as i_energy
from energy import get_nearest_neighbor_energy as nn_energy
from energy import get_nearest_neighbor_mismatch_energy as nn_mm_energy
from energy import get_symmetry_g as get_sym

# from graphic import visualize_genome as visualize

# string that will hold the dna sequence
prime_5 = ""
prime_3 = ""
# molarity of duplex and each strand alone
molarity_duplex = 0
molarity_p5 = 0
molarity_p3 = 0

"""
    Calculate melting temperature of dissociation based on nearest i
    neighbor method
"""


def find_melting_temperature():
    global prime_5, prime_3
    # iterates over shortest strand
    length = len(prime_5) if len(prime_5) < len(prime_3) else len(prime_3)
    # initial and terminal base pairs
    initial = prime_3[0]+prime_5[0]
    terminal = prime_3[length-1]+prime_5[length-1]
    # create tuple that will hold [enthalpy, entropy, gibb's free energy]
    sequence_data = [0.0, 0.0, 0.0]
    # calculate energy values of first pair ignoring nearest neighbor
    add(sequence_data, i_energy(initial, terminal))

    symmetrical = True

    for x in range(length-1):
        # check for symmetry
        if(prime_3[x]+prime_5[x] != prime_5[length-1-x]+prime_3[length-1-x]):
            symmetrical = False
        # calculates base-pair energy based on nearest_neighbor
        # find complementary base pair nearest neighbor energy
        if (not is_complement(prime_5[x], prime_3[x]) and x == 0) or \
                (not is_complement(prime_5[x+1], prime_3[x+1])
                    and x+1 == length):
            # We cannot evaluate a mismatch on the end of a sequence yet
            continue
        elif is_complement(prime_5[x+1], prime_3[x+1]) and \
                is_complement(prime_3[x], prime_5[x]):
            # print(prime_3[x]+prime_3[x+1] +", "+ prime_5[x]+prime_5[x+1])
            add(sequence_data, nn_energy(prime_3[x:x+2], prime_5[x:x+2]))
        # look for base pair inverse mismatch value
        else:
            # print(prime_3[x]+prime_3[x+1] +", "+ prime_5[x]+prime_5[x+1])
            add(sequence_data, nn_mm_energy(prime_3[x:x+2], prime_5[x:x+2]))

    # checks for terminal 5'-> 3' T-A base pair
    add(sequence_data, e_energy(initial, terminal))
    # if sequence is symmetrical, add symmetry value
    add(sequence_data, get_sym(symmetrical))

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
    AB_mol = 1
    # Mg_mol = .5
    R = 8.3144598
    temperature = -(sequence_data[2])/(R*math.log(AB_mol/2, math.exp(1)))

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

    return sequence_data+[length]+[temperature]


"""
   string from text file and creates DNA sequence for 5'->3' and 3'->5' strand
"""


def parse_genome(number):
    global prime_5, prime_3
    prime_5 = ""
    prime_3 = ""
    with open("dna_samples/dna_sample_"+str(number)+".txt") as f:
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


# Show help menu
def help():
    print("Commands:")
    print("")
    print("help: Show help menu")
    print("test: Run program on sequence and compare to expected results")
    print("show: Show sequence and expected result for that sequence")
    print("details: Explains experimental conditions and procedure\n")


# runs program on dna sequence and compares it to expected value
def test(number):
    parse_genome(number)
    data = find_melting_temperature()
    print("Length: " + str(round(data[3], 2)))
    print("Gibb's free energy: " + str(round(data[2], 2)))
    print("Entropy; " + str(round(data[1], 2)))
    print("Enthalpy: " + str(round(data[0], 2)))
    print("Temperature: " + str(round(data[4], 2)))
    # print correct result
    show_results(number)


# details about the program
def details():
    print("Parameters:")
    print("pH: 7")
    print("Mg+ molarity: .5")
    print("Oligonucleotide probe molarity: 1\n")
    print("Measurements:")
    print("Gibb's free energy: kcal/mol")
    print("Enthalpy: kcal/mol")
    print("Entropy: eu")
    print("Temperature: K\n")
    print("Goal:")
    print("To find the temperature at which half of a DNA duplex")
    print("dissociates 1/2 of its base pairs from the Gibb's energy")
    print("contribution of each base pair bond.\n")


# prints out a dna sequence from a sample file
def show(number):
    print("Sequence\n")
    # print dna sequence
    with open("dna_samples/dna_sample_"+str(number)+".txt") as f:
        print(f.read())
    # printing expected output for DNA sample
    show_results(number)


# show expected results on file
def show_results(num):
    print("Expected Output:")
    start = ((int(num)-1)*7)+1
    end = start+5
    s = str(start)+","+str(end)+"p"
    call(["sed", "-n", s, "dna_samples/catalog.txt"])


# show all sample sequences
def show_all(n):
    for x in range(n):
        show(x+1)


# test all samples
def test_all(n):
    for x in range(n):
            test(x+1)


# add tuples of size 3
def add(t1, t2):
    t1[0] += t2[0]
    t1[1] += t2[1]
    t1[2] += t2[2]


# test if bases are complements
def is_complement(b1, b2):
    if (b1 == 'T' and b2 == 'A') or (b1 == 'A' and b2 == 'T'):
        return True
    if (b1 == 'C' and b2 == 'G') or (b1 == 'G' and b2 == 'C'):
        return True
    return False


if __name__ == '__main__':

    details()

    # number of samples
    with open("dna_samples/catalog.txt") as f:
        samples = int(sum(1 for _ in f)/7)

    # Command prompt loop
    while True:
        command = (input("[DNA_compiler]: ")).strip()
        if command == "help":
            help()
        elif command == "test":
            file_name = (input("select a number from 1 and 25: ")).strip()
            if file_name.isdigit():
                test(file_name)
            elif file_name == "all":
                test_all(samples)
            else:
                print("Invalid input: an integer between 1 and 25 required")
                continue
        elif command == "show":
            file_name = (input("select a number from 1 and 25: ")).strip()
            if file_name.isdigit():
                show(file_name)
            elif file_name == "all":
                show_all(samples)
            else:
                print("Invalid input: an integer between 1 and 25 required")
                continue
        elif command == "exit":
            sys.exit()
        elif command == "details":
            details()
        else:
            print("Invalid command: enter help for command menu")
