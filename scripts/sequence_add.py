import random

"""
    Functions to add to DNA strands
"""

bases = ['A', 'C', 'G', 'T']
base_pairs = ["AT", "GC", "TA", "CG"]

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

    for x in range(int(number)):
        integer = random.randint(0,3)
        three += basepairs[integer][0]
        five += basepairs[integer][1]

    return three + '/' + five

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


def substitute_mismatch(strand):
    sub_index = random.randint(1, len(strand.split('/')[0]) - 2)

    # keep getting random base until it is different from present one
    base = bases[random.randint(0,3)]
    while strand[sub_index] == base:
        base = bases[random.randint(0,3)]

    return strand.split('/')[0][:sub_index] + base + strand.split('/')[0][sub_index+1:] + '/' + strand.split('/')[1]



