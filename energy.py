"""
    All possible base pairs delta G energy from https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
"""
base_stability_energy = {
    "TA" : -0.12,
    "TG" : -0.78,
    "CA" : -0.78,
    "CG" : -1.44,
    "AG" : -1.29,
    "CT" : -1.29,
    "AA" : -1.04,
    "TT" : -1.04,
    "AT" : -1.27,
    "GA" : -1.66,
    "TC" : -1.66,
    "CC" : -1.97,
    "GG" : -1.97,
    "AC" : -2.04,
    "GT" : -2.04,
    "GC" : -2.70
}

"""
    Gibb's free energy from each adjacent base pair using the nearest neighbor 
    approach from https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics
"""
nearest_neighbor_energy = {
    "AATT" : [-8.4, -23.6, -1.02],
    "ATTA" : [-6.5, -18.8, -0.73],
    "TAAT" : [-6.3, -18.5, -0.60],
    "CAGT" : [-7.4, -19.3, -1.38],
    "GTCA" : [-8.6, -23.0, -1.43],
    "CTGA" : [-6.1, -16.1, -1.16],
    "GACT" : [-7.7, -20.3, -1.46],
    "CGGC" : [-10.1, -25.5, -2.09],
    "GCCG" : [-11.1, -28.4, -2.28],
    "GGCC" : [-6.7, -15.6, -1.77],
    "AT" : [0.0, -9.0, 2.8],
    "CG" : [0.0, -5.9, 1.82],
    "SYM" : [0.0, -1.4, 0.4],
    "TA" : [0.4, 0, 0.4]
}

"""
    Gibb's free energy from each adjacent mismatch base pair using the nearest neighbor 
    approach with trimer surrounding mismatch.
    TODO: C-A, C-C, G-G, A-A, T-T, G-A, G-T
"""
nearest_neighbor_mismatch_energy = {
    #C-T 
    "ACTT" : [0.7, 0.2, 0.64],
    "ATTC" : [-1.2, -6.2, 0.73],
    "CCGT" : [-0.8, -4.5, 0.62],
    "CTGC" : [-1.5, -6.1, 0.40],
    "GCCT" : [2.3, 5.4, 0.62],
    "GTCC" : [5.2, 13.5, 0.98],
    "TCAT" : [1.2, -20.3, 0.97],
    "TTAC" : [1.0, 0.7, 0.75],

    #C-A
    "AATC" : [2.3, 4.6, 0.88],
    "ACTA" : [5.3, 14.6, 0.77],
    "CAGC" : [1.9, 3.7, 0.75],
    "CCGA" : [0.6, -0.6, 0.79],
    "GACC" : [5.2, 14.2, 0.81],
    "GCCA" : [-0.7, -3.8, 0.47],
    "TAAC" : [3.4, 8.0, 0.92],
    "TCAA" : [7.6, 20.2, 1.33],

    #G-A
    "AATG" : [-0.6, -2.3, 0.14],
    "AGTA" : [-0.7, -2.3, 0.02],
    "CAGG" : [-0.7, -2.3, 0.03],
    "CGGA" : [-4.0, -13.2, 0.11],
    "GACG" : [-0.6, -1.0, -0.25],
    "GGCA" : [0.5, 3.2, -0.52],
    "TAAG" : [0.7, 0.7, 0.42],
    "TGAA" : [3.0, 7.4, 0.74],

    #G-T
    "AGTT" : [1.0, 0.9, 0.71],
    "ATTG" : [-2.5, -8.3, 0.07],
    "CGGT" : [-4.1, -11.7, -0.47],
    "CTGG" : [-2.8, -8.0, -0.32],
    "GGCT" : [3.3, 10.4, 0.08],
    "GTCG" : [-4.4, -12.3, -0.59],
    "TGAT" : [-0.1, -1.7, 0.43],
    "TTAG" : [-1.3, -5.3, 0.34],

    #C-C
    "ACTC" : [0.0, -4.4, 1.33],
    "CCGC" : [-1.5, -7.2, 0.7],
    "GCCC" : [3.6, 8.9, 0.79],
    "TCAC" : [6.1, 16.4, 1.05],

    #A-A
    "AATA" : [1.2, 1.7, 0.61],
    "CAGA" : [-0.9, -4.2, 0.43],
    "GACA" : [-2.9, -9.8, 0.17],
    "TAAA" : [4.7, 12.9, 0.69],

    #G-G
    "AGTG" : [-3.1, -9.5, -0.13],
    "CGGG" : [-4.9, -15.3, -0.11],
    "GGCG" : [-6.0, -15.8, -1.11],
    "TGAG" : [1.6, 3.6, 0.44],

    #T-T
    "ATTT" : [-2.7, -10.8, 0.69],
    "CTGT" : [-5.0, -15.8, -0.12],
    "GTCT" : [-2.2, -8.4, 0.45],
    "TTAT" : [0.2, -1.5, 0.68],
}

# Find energy contribution of head and tail base pairs in DNA sequence
def get_initiation_energy(init_pair, term_pair):
    total = [0.0, 0.0, 0.0]
    # if either the initial or terminal base pair is C-G or G-C we add thes energy 
    print("initial energy: "+init_pair+" "+term_pair)
    if (init_pair == "GC" or term_pair == "CG") or (init_pair == "CG" or term_pair == "GC"):
        print("Adding GC val")
        total = add(total, nearest_neighbor_energy["CG"])
    # If intial and terminal base pair are A-T or T-A then we add a certain enrgy value
    if (init_pair == "AT" or term_pair == "TA") and (init_pair == "TA" or term_pair == "AT"):
        print("Adding AT val")
        total = add(total, nearest_neighbor_energy["AT"])
    print(total)
    return total

# get value for symmetry when sequence is symmetrical
def get_symmetry_g():
    return nearest_neighbor_energy["SYM"]

# evaluate energy contribution of the terminal pair of a sequence
def get_end_energy(init_pair, term_pair):
    total = [0.0,0.0,0.0];
    # If terminal pair is 3' -> 5' == A -> T
    if term_pair == "AT":
        print("Adding A-T energy")
        total=add(total,nearest_neighbor_energy["AT"])
    return total

# Find energy contribution of inner adjacent dimer (AT/TA) from nearest neigbor model
def get_nearest_neighbor_energy(pair_one, pair_two):
    print("nn")
    if pair_one+pair_two in nearest_neighbor_energy:
        #print(nearest_neighbor_energy[pair_one+pair_two][2])
        return nearest_neighbor_energy[pair_one+pair_two]
    # if base pair sequence isn't found we can look for its reverse which will have the same energy (TC/AG == GA/CT)
    elif pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0] in nearest_neighbor_energy:
        #print(nearest_neighbor_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]][2])
        return nearest_neighbor_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]]
    else:
        print("Nearest neighbor sequence: "+str(pair_one+pair_two)+ "not found")

# # Find energy contribution of inner base pair mismatch from trimer(ATA/TAT) with the mismatch as its middle pair
def get_nearest_neighbor_mismatch_energy(pair_one, pair_two):
    print("mm")
    if pair_one+pair_two in nearest_neighbor_mismatch_energy:
        #print(nearest_neighbor_mismatch_energy[pair_one+pair_two][2])
        return nearest_neighbor_mismatch_energy[pair_one+pair_two]
    # if base pair sequence isn't found we can look for its reverse which will have the same energy (TCT/AGC == CGA/TCT)
    elif pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0] in nearest_neighbor_mismatch_energy:
        #print(nearest_neighbor_mismatch_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]][2])
        return nearest_neighbor_mismatch_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]]
    else:
        print("Mismatch sequence: "+str(pair_one+pair_two)+ "not found")

# adding trimers -- tuples
def add(t1, t2):
    return [t1[0]+t2[0], t1[1]+t2[1], t1[2]+t2[2]]
