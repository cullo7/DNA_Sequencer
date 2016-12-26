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
    "AATT" : [-33.1, -92.9, -4.26],
    "ATTA" : [-30.1, -85.4, -3.67],
    "TAAT" : [-30.1, -89.1, -2.50],
    "CAGT" : [-35.6, -95.0, -6.12],
    "GTCA" : [-35.1, -93.7, -6.09],
    "CTGA" : [-32.6, -87.9, -5.40],
    "GACT" : [-34.3, -92.9, -5.51],
    "CGGC" : [-44.4, -113.8, -9.07],
    "GCCG" : [-41.0, -102.1, -9.36],
    "GGCC" : [-33.5, -83.3, -7.66],
    "AT" : [9.6, 17.2, 4.31],
    "GC" : [0.4, -11.7, 4.05]
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
def get_initiation_energy(base_one, base_two):
    if pair_one+pair_two in nearest_neighbor_energy:
        return nearest_neighbor_energy[pair_one+pair_two]
    # if base pair isn't found we can look for its reverse which will have the same energy (T/A == A/T)
    elif pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0] in nearest_neighbor_energy:
        return nearest_neighbor_energy[pair_two[1]+pair_two[0]+pair_one[1]+pair_one[0]]
    else:
        print("Initiation sequence: "+str(pair_one+pair_two)+ "not found")

# Find energy contribution of inner adjacent dimer (AT/TA) from nearest neigbor model
def get_nearest_neighbor_energy(pair_one, pair_two):
    print("nn")
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

