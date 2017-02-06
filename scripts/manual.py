"""
    Manually enter DNA sequences
"""

def manual():
    """
        Description: Manually enter DNA sequence (with / separating single strands)
            Ex: >> ACT/TGA
        Returns: energy, enthalpy, entropy and temperature values for duplex
    """
    print("Enter DNA sequence separated by '/', enter nothing, 'quit' or q to finish")
    sequence = input(">> ")
    while sequence != "" and sequence != "quit" and sequence != "q":
        if sanitize_sequence(sequence):
            seq = find_melting_temperature(Sequence(sequence.split('/')[0], sequence.split('/')[1], .0004))
            print("Energy: " + str(seq.energy))
            print("Enthalpy: " + str(seq.enthalpy))
            print("Entropy: " + str(seq.entropy))
            print("Temperature: " + str(seq.temperature))
        else:
            print("Invalid Input")
        sequence = input(">> ")

def sanitize_sequence(iput):
    """
        Sanitizing Input for manual-entered sequence
    """
    if len(iput.split('/')) != 2:
        return False
    if len(iput.split('/')[0]) != len(iput.split('/')[1]):
        return False
    for x in iput:
        if x != 'C' and x != 'G' and x != 'A' and x != 'T' and x != '/':
            return False
    return True


