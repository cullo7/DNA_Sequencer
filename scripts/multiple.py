from scripts.sequence_add import add_stretch
from scripts.misc_func import run_sequences
from scripts.sequence_add import add_mismatches

"""
    Enter multiple different sequences
"""

def multiple():
    """
        Description: Interface to run multiple sequence

        Output(file): Sequence credentials
    """
    print()
    print("Enter numbers to indicate what sequences you would like to test")
    print("Enter a number for the number of adjacent matches then M for a mismatch")
    print("This program will find the average melting temperature for G-C matches and A-T matches with")
    print("All possible mismatches. Putting R before the number will make the matches random")
    print()
    print("Example 1: 3M3 -- three matches , a mismatch and three more matches")
    print("Example 2: 4M2 -- four matches , a mismatch and two more matches")
    print("Example 3: 3 -- three matches")
    print("Example 4: R3MR3 -- three random matches , a mismatch and three more random matches")
    print()
    print("Enter 'quit', 'exit', or 'q' to exit")
    print("Would you like to store all results in the same file (or separate files)?")
    same_file = False
    choice = input("(y or n): ")
    print()
    if choice == "y" or choice == "yes":
        open("output/output.txt",'w+')
        same_file = True

    i = 1
    sequence = input(str(i) + ": ")
    while sequence != "exit" and sequence != "quit" and sequence != "q":
        if sanitize_m(sequence):
            s = seudo_sequence_converter(sequence)
            run_sequences(add_mismatches(s), sequence, same_file)
        i += 1
        sequence = input(str(i) + ": ")

def sanitize_m(input_m):
    """
        Description: Sanitize multiple sequence input

        Returns: True if the input is valid, False otherwise

        Input: Suedo-sequence for multiple sequence generation
    """
    for x, item in enumerate(input_m):  
        if item != "M" and item != "R" and not item.isdigit():
            print("item: "+ str(item))
            print("Input must be either a n 'M', 'R', or a number")
            return False
        elif item == "M":
            if x == 0 or x == len(input_m)-1:
                print("Mismatch cannot be on the end")
                return False
        elif item == "R":
            if x == len(input_m)-1 or not input_m[x+1].isdigit():
                print("R cannot be the last character and must be followed by a digit")
                return False
    return True

def seudo_sequence_converter(sequence):
    """
        Descripton: Creates sequences from seudo-sequences and adds them to a list

        Returns: List of generated sequences

        Input: Seudo-sequence, ex: R3MR3 --three random basepairs, mismatch, three random base pairs
    """
    # [3'-high, 5'-high, 3'-low, 5'-low]
    strands = ["","","",""]
    random = False
    sequences = []
    x = 0
    while x < len(sequence):	
        if sequence[x] == "R":
            random = True
        elif sequence[x].isdigit() and random:
            temp = sequence[x]
            i = 1
            for i in range(1, len(sequence)-x+2):
                if not sequence[x:x+i].isdigit():
                    break
                else:
                    temp = sequence[x:x+i]

            add_stretch_rand(strands, temp)
            x+= (i-2)
            random = False
        elif sequence[x].isdigit():
            temp = sequence[x]
            i = 1
            for i in range(1, len(sequence)-x+2):
                if not sequence[x:x+i].isdigit():
                    break
                else:
                    temp = sequence[x:x+i]

            add_stretch(strands, temp)
            x+= (i-2)
        elif sequence[x] == "M":
            strands[0] += "M"
            strands[1] += "M"
            strands[2] += "M"
            strands[3] += "M"
        else:
            print("seudo-sequence character not recognized")
        x+=1
    sequences.append(strands[0] + '/' + strands[1])
    sequences.append(strands[2] + '/' + strands[3])
    return sequences


