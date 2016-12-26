"""
    Compiler for a DNA duplex to calculate the dissociation and association melting temperature 
"""

import sys
from subprocess import call
from energy import get_initiation_energy as i_energy
from energy import get_end_energy as e_energy
from energy import get_nearest_neighbor_energy as nn_energy
from energy import get_nearest_neighbor_mismatch_energy as nn_mm_energy
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
    #create tuple that will hold [enthalpy, entropy, gibb's free energy]
    sequence_data = [0.0,0.0,0.0] 
    # calculate energy values of first pair ignoring nearest neighbor
    sequence_data = add(sequence_data, i_energy(prime_5[0]+prime_3[0],prime_5[length-1]+prime_3[length-1]))

    for x in range(length-1):
        # calculates base-pair energy based on nearest_neighbor
        # if base pair is complementary then we find its nearest neighbor energy
        if not is_complement(prime_5[x], prime_3[x]) and x != 0:
            continue  
        elif is_complement(prime_5[x+1], prime_3[x+1]):
            print(prime_5[x+1] +", "+ prime_3[x+1])
            sequence_data = add(sequence_data, nn_energy(prime_5[x:x+2],prime_3[x:x+2]))
        # if base pair isn't a complement we look for its mismatch value as long as it is not the 
        # end of the sequence
        elif x+2 < length:
            sequence_data = add(sequence_data, nn_mm_energy(prime_5[x:x+3],prime_3[x:x+3]))
        # We cannot evaluate a mismatch on the end of a sequence
        else:
            print("Error: unable to evaluate mismatch on the end of a sequence")
            return -1


    # calculate free energy, entropy and enthalpy of last pair ignoring nearest neighbor
    sequence_data = add(sequence_data, e_energy(prime_5[0]+prime_3[0],prime_5[length-1]+prime_3[length-1]))
    
    # Account for Other factors that will affect temperature/gibbs energy
    # Calculate temperature from Gibbs Energy

    return sequence_data+[length]

"""
   string from text file and creates DNA sequence for 5'->3' and 3'->5' strand
"""
def parse_genome(number):
    global prime_5, prime_3
    with open("dna_samples/dna_sample_"+str(number)+".txt") as f:
        seq = f.read()
    seq = seq.split()
    even = True
    for x in range(len(seq)):
        if seq[x] != '/n' and seq[x] != ' ':
            if even:
                prime_3+=seq[x]
                even = False
            else:
                prime_5+=seq[x]
                even = True

# Show help menu
def help():
    print("Commands:")
    print("")
    print("help: show help menu")
    print("test: run program on sequence and compare to expected results")
    print("show: show sequence and expected result for that sequence")

# runs program on dna sequence and compares it to expected value from catalog.txt file
def test(number):
    parse_genome(number)
    data = find_melting_temperature()
    print("Length: "+ str(data[3]))
    print("Gibb's free energy: "+ str(data[2]))
    print("Entropy; "+ str(data[1]))
    print("Enthalpy: "+ str(data[0]))
    # print correct result
    show_results(number)


# prints out a dna sequence from a sample file
def show(number):
    print("Sequence\n")
    # print dna sequence
    with open("dna_samples/dna_sample_"+str(number)+".txt") as f:
        print(f.read())
    # printing expected output for DNA sample
    show_results(number)
    
# show expected results on file
def show_results(num):
    print("Expected Output:")
    start = ((int(num)-1)*5)+1
    end =  start+4
    s = str(start)+","+str(end)+"p"
    call(["sed", "-n",s , "dna_samples/catalog.txt"])

# show all sample sequences
def show_all(n):
    for x in range(n):
        show(x+1)

# test all samples
def test_all(n):
    for x in range(n):
            test(x+1)
    
# add tuples of size 3
def add(t1, t2):
    return [t1[0]+t2[0], t1[1]+t2[1], t1[2]+t2[2]]

# test if bases are complements
def is_complement(b1, b2):
    if (b1 == 'T' and b2 == 'A') or (b1 == 'A' and b2 == 'T'):
        return True
    if (b1 == 'C' and b2 == 'G') or (b1 == 'G' and b2 == 'C'):
        return True
    return False

if __name__ == '__main__':

    print("Units:")
    print("Gibb's free energy: kcal/mol")
    print("Enthalpy: kcal/mol")
    print("Entropy: eu\n")

    # nnumber of samples
    samples = 5

    # Command prompt loop
    while True:
        command = input("[DNA_compiler]: ") 
        if  command == "help":
            help()
        elif command == "test":
            file_name = input("select a number between 1 and 25, inclusive: ")
            if file_name.isdigit():
                test(file_name)
            elif file_name == "all":
                test_all(samples)
            else:
                print("Invalid input: an integer between 1 and 25 required")
                continue
        elif command == "show":
            file_name = input("select a number between 1 and 25, inclusive: ")
            if file_name.isdigit():
                show(file_name)
            elif file_name == "all":
               show_all(samples)
            else:
                print("Invalid input: an integer between 1 and 25 required")
                continue
        elif command == "exit":
            sys.exit()
        else:
            print("Invalid command: enter help() for command menu")
