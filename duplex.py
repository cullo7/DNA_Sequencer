class Sequence(object):
    """ Object that holds 3' and 5' sequence for DNA duplex """

    # Set all variables so we don't crash program by trying 
    # to access a variable that hasn't be initialized
    def __init__(self, p3, p5, mol):
        self.three_prime = p3
        self.five_prime = p5
        self.set_length(len(p3))
        self.terminal = self.three_prime[self.length-1] + self.five_prime[self.length-1]
        self.initial = self.three_prime[0] + self.five_prime[0]
        self.energy = 0
        self.enthalpy = 0
        self.entropy = 0

    """ Setters """

    def set_three_prime(self, seq):
        self.three_prime = seq
        self.set_length(len(seq))

    def set_five_prime(self, seq):
        self.five_prime = seq

    def set_complementary(self, comp):
        self.complementary = comp

    # Sets molarity of dna oligonucleotide in sample
    def set_oligo_molarity(self, mol):
        self.oligo_molarity = mol

    def set_salt_molarity(self, mol):
        self.salt_molarity = mol

    def set_energy(self, en):
        self.energy = en

    def set_enthalpy(self, enth):
        self.enthalpy = enth

    def set_entropy(self, ent):
        self.entropy = ent
    
    def set_length(self, l):
        self.length = l

    def set_temperature(self, temp):
        self.temperature = temp

    """ Adders """

    def add_energy(self, en):
        self.energy = round(self.energy + en, 2)

    def add_enthalpy(self, enth):
        self.enthalpy = round(self.enthalpy + enth, 2)

    def add_entropy(self, ent):
        self.entropy = round(self.entropy + ent, 2)

    # adding tuple of enthalphy, entropy and free energy
    def add(self, tup):
        self.add_enthalpy(tup[0])
        self.add_entropy(tup[1])
        self.add_energy(tup[2])

    """ Getters """

    def get_3(self, index):
        return self.three_prime[index]

    def get_5(self, index):
        return self.five_prime[index]

    # get a sequence from 3' strand
    def get_3_s(self, start, finish):
        return self.three_prime[start:finish]

    # get a sequence from 5' strand
    def get_5_s(self, start, finish):
        return self.five_prime[start:finish]

    def get_oligo_molarity(self):
        if self.complementary:
            return 0.0001
        else:
            return 0.0004

    def get_complementary(self):
        return self.complementary

    def get_salt_molarity(self, mol):
        return self.salt_molarity

    def get_energy(self):
        return self.energy

    def get_enthalpy(self):
        return self.enthalpy

    def get_entropy(self):
        return self.entropy
    
    def get_length(self):
        return self.length

    def get_temperature(self):
        return self.temperature

    def get_terminal(self):
        return self.terminal

    def get_initial(self):
        return self.initial

    # Will probably move this to visualize
    def print_out(self):
        print("Complemtary: " + self.complementary)
        print("Symmetric: " + self.symmetry)
        print("Five Prime Sequence: " + self.five_prime)
        print("Three Prime Sequence: " + self.three_prime)
        print("Duplex Molarity: " + self.oligo_molarity)
        print("Salt Molarity: " + self.salt_molarity)
        print("Length: " + self.length)
        print("Energy: " + self.energy)
        print("Enthalpy: " + self.enthalpy)
        print("Entropy: " + self.entropy)
