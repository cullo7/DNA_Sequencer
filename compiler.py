"""
    Compiler for a DNA duplex to calculate the dissociation and
    association melting temperature
"""

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
        t = s.get_energy()
        print(s.get_3_s(x, x + 2) + " " + s.get_5_s(x, x + 2))
        print("E: " + str(round((t - t1), 2)))
        print(s.get_energy())
        t1 = t
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

    print("========")
    print(s.get_enthalpy())
    print(s.get_entropy())
    print(s.get_oligo_molarity())
    print("========")

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
    print("details: Explains experimental conditions and procedure\n")


# runs program on dna sequence and compares it to expected value
def test(number):
    sequence = parse_genome(number)
    data = find_melting_temperature(sequence)
    print("Length: " + str(data.length))
    print("Gibb's free energy: " + str(data.energy))
    print("Entropy; " + str(data.entropy))
    print("Enthalpy: " + str(data.enthalpy))
    print("Temperature: " + str(data.temperature))
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


# test if bases are complements
def is_complement(b1, b2):
    if (b1 == 'T' and b2 == 'A') or (b1 == 'A' and b2 == 'T'):
        return True
    if (b1 == 'C' and b2 == 'G') or (b1 == 'G' and b2 == 'C'):
        return True
    return False


if __name__ == '__main__':

    if len(sys.argv) > 2:
        print("Program only accepts two total command line arguments")
        print("usage: python compiler.py")
        sys.exit(0)

    details()

    # number of samples
    with open("dna_samples/catalog.txt") as f:
        samples = int(sum(1 for _ in f) / 7)

    # Command prompt loop
    while True:
        command = input("[DNA_compiler]: ").strip().split()
        if command[0] == "help" or command[0] == "h":
            help()
        elif command[0] == "test" or command[0] == "t":
            if len(command) > 1:
                if command[1] == "all":
                    test_all(samples)
                if not command[1].isdigit():
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
                if not command[1].isdigit():
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
        elif command[0] == "exit" or command[0] == "e" or command[0] == "quit":
            sys.exit()
        elif command[0] == "details" or command[0] == "d":
            details()
        else:
            print("Invalid command: enter help for command menu")
