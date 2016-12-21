from sys import argv

prime_5 = ""
prime_3 = ""

def find_melting_temperature(prime_5, prime_3):
    length = len(prime_5) if len(prime_5) < len(prime_3) else len(prime_3)
    for x in range(length):
        if(prime_5[x] != prime_3[x]):
            return False
    return True

"""
   string from text file and creates DNA sequence for 5'->3' and 3'->5' strand
"""
def parse_genome(dna_string):
    global prime_5, prime_3


if __name__ == '__main__':
    # Read genome from file
    script, filename = argv

    txt = open(filename)

    print ("Here's your file %r:" % (filename))
    print (txt.read())

    print ("Type the filename again:")
    file_again = input("> ")

    txt_again = open(file_again)

    print (txt_again.read())

    parse_genome(txt_again)
