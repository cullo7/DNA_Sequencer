"""
    Functions for energy additions such as nearest neighbor contribution.
    symmetry, and initial base pair
"""

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


def get_initiation_energy(init_pair, term_pair, comp):
    """
        Description: Find energy contribution of head and tail base pairs in DNA sequence

        Returns: Total energy, enthalpy, and entropy contribution

        Input: init_pair - initial base pair, term_pair - terminal base pair, 
            comp - sentinal to specifcy complementarity
    """
    total = [0.0, 0.0, 0.0]
    """
        Check if initial is A-T or G-C base pair
    """
    if init_pair == "AT" or init_pair == "TA":
        add(total, initial_energy["AT"])
    elif init_pair == "GC" or init_pair == "CG":
        add(total, initial_energy["CG"])

    """
        Check if terminal pair is A-T or G-C base pair 
    """
    if term_pair == "AT" or term_pair == "TA":
        add(total, initial_energy["AT"])
    elif term_pair == "CG" or term_pair == "GC":
        add(total, initial_energy["CG"])

    """
        Check for terminal 3'->5' A-T
    """
    if term_pair == "AT":
        add(total, initial_energy["TERM_AT"]) 
    return total


def get_symmetry_g(sym):
    """ 
        Description: Retrieves the contribution of a symmetrical duplex

        Returns: Value contributions of a symmetrical duplex
    """
    if(sym):
        return initial_energy["SYM"]
    else:
        return [0, 0, 0]


def get_nearest_neighbor_energy(pair_one, pair_two):
    """
        Description: Find energy contribution of inner adjacent dimer (AT/TA)
        from nearest neigbor model

        Return: Nearest neighbor contribution of adjacent base pairs

        Input: pair_one - first base pair, pair_two - second, adjacent base pair
    """
    regular = pair_one + pair_two
    inverse = pair_two[1] + pair_two[0] + pair_one[1] + pair_one[0]

    if regular in nearest_neighbor_energy:
        return nearest_neighbor_energy[regular]
    # if sequence isn't found we use its inverse (e. g. TC/AG == GA/CT)
    elif inverse in nearest_neighbor_energy:
        return nearest_neighbor_energy[inverse]
    else:
        print("Nearest neighbor sequence: " + str(regular) + "not found")


def get_nearest_neighbor_mismatch_energy(pair_one, pair_two):
    """
        Description: Energy contribution from mismatch trimer(e.g. ATA/TTT)
        
        Return: Nearest neighbor contribution of adjacent mismatch base pairs

        Input: pair_one - first base pair, pair_two - second, adjacent base pair
    """
    regular = pair_one + pair_two
    inverse = pair_two[1] + pair_two[0] + pair_one[1] + pair_one[0]
    if regular in nearest_neighbor_mismatch_energy:
        return nearest_neighbor_mismatch_energy[regular]
    # look for equivalent inverse sequence (e.g. TCT/AGC == CGA/TCT)
    elif inverse in nearest_neighbor_mismatch_energy:
        return nearest_neighbor_mismatch_energy[inverse]
    else:
        print("Mismatch sequence: " + str(regular) + "not found")


def add(t1, t2):
    """ 
        Description: Adding trimers -- tuples
    """
    t1[0] = round(t1[0] + t2[0], 2)
    t1[1] = round(t1[1] + t2[1], 2)
    t1[2] = round(t1[2] + t2[2], 2)
