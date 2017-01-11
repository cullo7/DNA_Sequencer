"""
    All possible base pairs delta G energy from
    https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
"""
base_stability_energy = {
    "TA": -0.12,
    "TG": -0.78,
    "CA": -0.78,
    "CG": -1.44,
    "AG": -1.29,
    "CT": -1.29,
    "AA": -1.04,
    "TT": -1.04,
    "AT": -1.27,
    "GA": -1.66,
    "TC": -1.66,
    "CC": -1.97,
    "GG": -1.97,
    "AC": -2.04,
    "GT": -2.04,
    "GC": -2.70
}

"""
    Gibb's free energy from each adjacent base pair using the nearest neighbor
    approach from https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
"""
backup_nearest_neighbor_energy = {
    "AATT": [-8.4, -23.6, -1.02],
    "ATTA": [-6.5, -18.8, -0.73],
    "TAAT": [-6.3, -18.5, -0.60],
    "CAGT": [-7.4, -19.3, -1.38],
    "GTCA": [-8.6, -23.0, -1.43],
    "CTGA": [-6.1, -16.1, -1.16],
    "GACT": [-7.7, -20.3, -1.46],
    "CGGC": [-10.1, -25.5, -2.09],
    "GCCG": [-11.1, -28.4, -2.28],
    "GGCC": [-6.7, -15.6, -1.77]
}

nearest_neighbor_energy = {
    "AATT": [-7.9, -22.2, -1.00],
    "ATTA": [-7.2, -20.4, -0.88],
    "TAAT": [-7.2, -21.3, -0.58],
    "CAGT": [-8.5, -22.7, -1.45],
    "GTCA": [-8.4, -22.4, -1.44],
    "CTGA": [-7.8, -21.0, -1.28],
    "GACT": [-8.2, -22.2, -1.30],
    "CGGC": [-10.6, -27.2, -2.17],
    "GCCG": [-9.8, -24.4, -2.24],
    "GGCC": [-8.0, -19.9,-1.84]
}

"""
    Gibb's free energy from each adjacent mismatch base pair using the nearest
    neighbor approach with trimer surrounding mismatch.
    TODO: C-A, C-C, G-G, A-A, T-T, G-A, G-T
"""
nearest_neighbor_mismatch_energy = {
    # C-T
    "ACTT": [0.7, 0.2, 0.64],
    "ATTC": [-1.2, -6.2, 0.73],
    "CCGT": [-0.8, -4.5, 0.62],
    "CTGC": [-1.5, -6.1, 0.40],
    "GCCT": [2.3, 5.4, 0.62],
    "GTCC": [5.2, 13.5, 0.98],
    "TCAT": [1.2, 0.7, 0.97],
    "TTAC": [1.0, 0.7, 0.75],

    # C-A
    "AATC": [2.3, 4.6, 0.88],
    "ACTA": [5.3, 14.6, 0.77],
    "CAGC": [1.9, 3.7, 0.75],
    "CCGA": [0.6, -0.6, 0.79],
    "GACC": [5.2, 14.2, 0.81],
    "GCCA": [-0.7, -3.8, 0.47],
    "TAAC": [3.4, 8.0, 0.92],
    "TCAA": [7.6, 20.2, 1.33],

    # G-A
    "AATG": [-0.6, -2.3, 0.14],
    "AGTA": [-0.7, -2.3, 0.02],
    "CAGG": [-0.7, -2.3, 0.03],
    "CGGA": [-4.0, -13.2, 0.11],
    "GACG": [-0.6, -1.0, -0.25],
    "GGCA": [0.5, 3.2, -0.52],
    "TAAG": [0.7, 0.7, 0.42],
    "TGAA": [3.0, 7.4, 0.74],

    # G-T
    "AGTT": [1.0, 0.9, 0.71],
    "ATTG": [-2.5, -8.3, 0.07],
    "CGGT": [-4.1, -11.7, -0.47],
    "CTGG": [-2.8, -8.0, -0.32],
    "GGCT": [3.3, 10.4, 0.08],
    "GTCG": [-4.4, -12.3, -0.59],
    "TGAT": [-0.1, -1.7, 0.43],
    "TTAG": [-1.3, -5.3, 0.34],

    # C-C
    "ACTC": [0.0, -4.4, 1.33],
    "CCGC": [-1.5, -7.2, 0.7],
    "GCCC": [3.6, 8.9, 0.79],
    "TCAC": [6.1, 16.4, 1.05],

    # A-A
    "AATA": [1.2, 1.7, 0.61],
    "CAGA": [-0.9, -4.2, 0.43],
    "GACA": [-2.9, -9.8, 0.17],
    "TAAA": [4.7, 12.9, 0.69],

    # G-G
    "AGTG": [-3.1, -9.5, -0.13],
    "CGGG": [-4.9, -15.3, -0.11],
    "GGCG": [-6.0, -15.8, -1.11],
    "TGAG": [1.6, 3.6, 0.44],

    # T-T
    "ATTT": [-2.7, -10.8, 0.69],
    "CTGT": [-5.0, -15.8, -0.12],
    "GTCT": [-2.2, -8.4, 0.45],
    "TTAT": [0.2, -1.5, 0.68]
}


"""
    Initial and terminal base pair contributions for
    mismatch DNA duplexes
"""
initial_energy = {
   "CG": [0.1, -2.8, 0.98],
   "AT": [2.3, 4.1, 1.03],
   "TERM_AT": [0.4, 0, 0.4],
   "SYM": [0, -1.4, 0.4]
}

# Find energy contribution of head and tail base pairs in DNA sequence
def get_initiation_energy(init_pair, term_pair, comp):
    total = [0.0, 0.0, 0.0]
    print("init: "+init_pair+" term pair: "+term_pair)
    """
        Check if initial is A-T or G-C base pair
    """
    if init_pair == "AT" or init_pair == "TA":
        print("Adding AT val")
        add(total, initial_energy["AT"])
    elif init_pair == "GC" or init_pair == "CG":
        print("Adding CG val")
        add(total, initial_energy["CG"])
    """
        Check if terminal pair is A-T or G-C base pair 
    """
    if term_pair == "AT" or term_pair == "TA":
        print("Adding AT val")
        add(total, initial_energy["AT"])
    elif term_pair == "CG" or term_pair == "GC":
        print("Adding CG val")
        add(total, initial_energy["CG"])
    """
        Check for terminal 3'->5' A-T
    """
    if term_pair == "AT":
        print("Adding A-T energy")
        add(total, initial_energy["TERM_AT"]) 
    return total


# get value for symmetry when sequence is symmetrical
def get_symmetry_g(sym):
    if(sym):
        print("Adding Symmetry value")
        return initial_energy["SYM"]
    else:
        return [0, 0, 0]


"""
    Find energy contribution of inner adjacent dimer (AT/TA)
    from nearest neigbor model
"""


def get_nearest_neighbor_energy(pair_one, pair_two):
    #print("nn")
    regular = pair_one + pair_two
    inverse = pair_two[1] + pair_two[0] + pair_one[1] + pair_one[0]
    if regular in nearest_neighbor_energy:
        # print(nearest_neighbor_energy[pair_one+pair_two][2])
        return nearest_neighbor_energy[regular]
    # if sequence isn't found we use its inverse (e. g. TC/AG == GA/CT)
    elif inverse in nearest_neighbor_energy:
        # print(nearest_neighbor_energy[inverse[2])
        return nearest_neighbor_energy[inverse]
    else:
        print("Nearest neighbor sequence: " + str(regular) + "not found")


# Energy contribution from mismatch trimer(e.g. ATA/TTT)
def get_nearest_neighbor_mismatch_energy(pair_one, pair_two):
    #print("mm")
    regular = pair_one + pair_two
    inverse = pair_two[1] + pair_two[0] + pair_one[1] + pair_one[0]
    if regular in nearest_neighbor_mismatch_energy:
        # print(nearest_neighbor_mismatch_energy[pair_one+pair_two][2])
        return nearest_neighbor_mismatch_energy[regular]
    # look for equivalent inverse sequence (e.g. TCT/AGC == CGA/TCT)
    elif inverse in nearest_neighbor_mismatch_energy:
        # print(nearest_neighbor_mismatch_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]][2])
        return nearest_neighbor_mismatch_energy[inverse]
    else:
        print("Mismatch sequence: " + str(regular) + "not found")


# adding trimers -- tuples
def add(t1, t2):
    t1[0] = round(t1[0] + t2[0], 2)
    t1[1] = round(t1[1] + t2[1], 2)
    t1[2] = round(t1[2] + t2[2], 2)
