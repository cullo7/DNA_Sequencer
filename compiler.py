"""
    Compiler for a DNA duplex to calculate the dissociation and association melting temperature 
"""

from sys import argv
from energy import get_initiation_energy as i_energy
from energy import get_nearest_neighbor_energy as nn_energy
#from graphic import visualize_genome as visualize

# string that will hold the dna sequence 
prime_5 = ""
prime_3 = ""
# molarity of duplex and each strand alone
molarity_duplex = 0
molarity_p5 = 0
molarity_p3 = 0

"""
  Calculate melting temperature of dissociation based on nearest neighbor method
"""
def find_melting_temperature():
    global prime_5, prime_3
    #iterates over shortest strand
    length = len(prime_5) if len(prime_5) < len(prime_3) else len(prime_3)
    delta_g = 0
    # calculate free energy of first pair ignoring nearest neighbor
    delta_g += i_energy(prime_5[0],prime_3[0])

    for x in range(length-1):
        # For right now, we will exit is we get a mismatch
        # In the future mismatch will correspond to different affect on Gibbs enrgy/temperature
        if(prime_5[x] != prime_3[x]):
            return 0
        # calculates base-pair energy based on nearest_neighbor
        delta_g += nn_energy(prime_5[x:x+1],prime_3[x:x+1])

    # calculate free energy of last pair ignoring nearest neighbor
    delta_g += i_energy(prime_5[length-1],prime_3[length-1])
    
    # Account for Other factors that will affect temperature/gibbs energy
    # Calculate temperature from Gibbs Energy

    return delta_g

"""
   string from text file and creates DNA sequence for 5'->3' and 3'->5' strand
"""
def parse_genome(dna_string):
    global prime_5, prime_3
    


if __name__ == '__main__':
    # Read genome from file
    filename = argv[1]

    txt = open(filename)

    print ("Here's your file %r:" % (filename))
    print (txt.read())

    print ("Type the filename again:")
    file_again = input("> ")

    txt_again = open(file_again)

    print (txt_again.read())

    parse_genome(txt_again)

    m_temp = find_melting_temperature()

    if m_temp == 0:
        print("Mismatch")
    else:   
        print("melting temperature: "+ m_temp)
