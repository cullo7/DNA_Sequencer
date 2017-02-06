from subprocess import call
from scripts.duplex import Sequence
from scripts.genome import parse_genome
from scripts.genome import find_melting_temperature

"""
    Tests algorithm on sample cases
"""


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


