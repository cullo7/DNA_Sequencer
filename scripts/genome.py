from scripts.misc_func import is_complement
from scripts.duplex import Sequence
from scripts.energy import get_initiation_energy as i_energy
from scripts.energy import get_nearest_neighbor_energy as nn_energy
from scripts.energy import get_nearest_neighbor_mismatch_energy as nn_mm_energy
from scripts.energy import get_symmetry_g as get_sym

import math

"""
    Parses and finds temperature for dna strand
"""

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

    if int(number) > 15 or int(number) < 0:
        raise ValueError("Number has to be less than 15 and greater than 0")

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


