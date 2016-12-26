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
    approach with trimer surrounding mismatch. Mismatch pair was restricted to C-T
"""
nearest_neighbor_mismatch_energy = {
# seven linearly independent trimers
"ACCTTG" : [5.9, 13.8, 1.62], 
"CCCGTG" : [4.4, 9.0, 1.60], 
"GCACTT" : [3.3, 6.2, 1.37], 
"GCCCTG" : [7.5, 19.0, 1.60], 
"GCGCTC" : [0.8, -0.7, 1.02], 
"GCTCTA" : [1.1, -0.8, 1.35],
"TCCATG" : [6.4, 14.3, 1.95],

# nine other trimer contexts
"ACATTT" : [1.7, 1.0, 1.39], 
"ACGTTC" : [-0.8, -5.9, 1.04],
"ACTTTA" : [-0.5, -6.0, 1.37],
"CCAGTT" : [0.2, -3.8, 1.37], 
"CCGGTC" : [-2.3, -10.7, 1.02], 
"CCTGTA" : [-2.0, -10.8, 1.35], 
"TCAATT" : [2.2, 1.5, 1.72], 
"TCGATC" : [-0.3, -5.4, 1.37], 
"TCTATA" : [0.0, -5.5, 1.70],
}

# Find energy contribution of head and tail base pairs in DNA sequence
def get_initiation_energy(pair_one, pair_two):
    # if either the initial or terminal base pair is C-G or G-C we add thes energy 
    if (pair_one == "GC" or pair_one == "CG") or (pair_two == "GC" or pair_two == "CG"):
        return nearest_neighbor_energy["CG"]
    # If intial and terminal base pair are A-T or T-A then we add a certain enrgy value
    elif (pair_one == "AT" or pair_one == "TA") and (pair_two == "TA" or pair_two == "AT"):
        return nearest_neighbor_energy["AT"]
    else:
        print("Initiation sequence: "+str(pair_one+pair_two)+ "not found")
        return [0.0,0.0,0.0]

def get_end_energy(init_pair, term_pair):
    total = [0.0,0.0,0.0];
    if (init_pair == "AT" and term_pair == "TA") or (init_pair == "GC" and term_pair == "CG") or (init_pair == "TA" and term_pair == "AT") or (init_pair == "CG" and term_pair == "GC"):
        total=add(total,nearest_neighbor_energy["SYM"])
    if term_pair == "AT":
        total=add(total,nearest_neighbor_energy["TA"])
    return total

# Find energy contribution of inner adjacent dimer (AT/TA) from nearest neigbor model
def get_nearest_neighbor_energy(pair_one, pair_two):
    if pair_one+pair_two in nearest_neighbor_energy:
        return nearest_neighbor_energy[pair_one+pair_two]
    # if base pair sequence isn't found we can look for its reverse which will have the same energy (TC/AG == GA/CT)
    elif pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0] in nearest_neighbor_energy:
        return nearest_neighbor_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]]
    else:
        print("Nearest neighbor sequence: "+str(pair_one+pair_two)+ "not found")

# # Find energy contribution of inner base pair mismatch from trimer(ATA/TAT) with the mismatch as its middle pair
def get_nearest_neighbor_mismatch_energy(pair_one, pair_two):
    if pair_one+pair_two in nearest_neighbor_mismatch_energy:
        return nearest_neighbor_mismatch_energy[pair_one+pair_two]
    # if base pair sequence isn't found we can look for its reverse which will have the same energy (TCT/AGC == CGA/TCT)
    elif pair_two[2]+pair_two[1]+pair_two[0]+pair_one[2]+pair_one[1]+pair_one[0] in nearest_neighbor_mismatch_energy:
        return nearest_neighbor_mismatch_energy[pair_two[2]+pair_two[1]+pair_two[0]+pair_one[2]+pair_one[1]+pair_one[0]]
    else:
        print("Mismatch sequence: "+str(pair_one+pair_two)+ "not found")

def add(t1, t2):
    return [t1[0]+t2[0], t1[1]+t2[1], t1[2]+t2[2]]
