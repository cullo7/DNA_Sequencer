from termcolor import colored

"""
    File for help menus and graphics
"""

def title():
    """
        Description: Prints introduction title
    """
    width = 125
    red = 'red'
    green = 'green'
    print('\n')
    print(colored('#' * 127, green))
    print('\n')
    print(colored(" /=\              ______  __   __     ___                    /=\  ", red).center(width, ' '))
    print(colored("/===\             |  _  \ | \ | |    / _ \                  /===\ ", red).center(width, ' '))
    print(colored("|===|             | | | | |  \| |   / /_\ \                 |===| ", red).center(width, ' '))
    print(colored("\===/             | |_| | | |\  |  / /___\ \                \===/ ", red).center(width, ' '))
    print(colored(" \=/              |_____/ |_| \_| /_/     \_\                \=/  ", red).center(width, ' '))
    print(colored("  X     ___   ___   _      _  ____   _   _     ___   ___      X   ", red).center(width, ' '))
    print(colored(" /=\   /  _\ /   \ | \    / | |   | | | | |   | __| | . |    /=\  ", red).center(width, ' '))
    print(colored("/===\  | |   | | | |  \  /  | | | | | | | |   | |_  | __|   /===\  ", red).center(width, ' '))
    print(colored("|===|  | |   | | | |   \/   | | __| | | | |   |  _| |  \    |===| ", red).center(width, ' '))
    print(colored("\===/  | |_  | | | | |\__/| | | |   | | | |_  | |_  | | \   \===/ ", red).center(width, ' '))
    print(colored(" \=/   \___/ \___/ |_|    |_| |_|   |_| |___| |___| |_||_\   \=/  ", red).center(width, ' '))
    print('\n')
    print(colored('#' * 127, green))

def help():
    """
        Description: Prints help menu
    """
    print("Commands:")
    print("")
    print("help: Show help menu")
    print("test: Run program on sequence and compare to expected results")
    print("show: Show sequence and expected result for that sequence")
    print("details: Explains experimental conditions and procedure")
    print("multiple: Enter menu for running and recording multiple sequences")
    print("difference: Enters menu for dining average difference of mismatch substitution")
    print("manual: Enters mode where user can manually enter sequence")
    print()

def details():
    """
        Description: Prints details about the program
    """
    print("\nParameters:\n")
    print("pH: 7")
    print("oligonucleotide: 1E-4 for self-complementary sequences and 4E-4 for all others")
    print("\nMeasurements:\n")
    print("Gibb's free energy: kcal/mol")
    print("Enthalpy: kcal/mol")
    print("Entropy: eu")
    print("Temperature: K")
    print("\nGoal:\n")
    print("To find the temperature at which half of a DNA duplex")
    print("dissociates 1/2 of its base pairs from the Gibb's energy")
    print("contribution of each base pair bond.\n")



