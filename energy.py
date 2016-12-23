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

def get_initiation_energy(base_one, base_two):
    return nearest_neighbor_energy[base_one+base_two][2]

def get_nearest_neighbor_energy(pair_one, pair_two):
    return nearest_neighbor_energy[pair_one+pair_two][2]
