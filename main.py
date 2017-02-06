"""
    Compiler for a DNA duplex to calculate the dissociation and
    association melting temperature with one or two base pair
    mismatches
"""

# Library Imports
import numpy
import sys
from subprocess import call

# Local Imports

from scripts.menu import title, help, details
from scripts.difference import difference
from scripts.multiple import multiple
from scripts.tests import test
from scripts.manual import manual

"""
    Global variables
"""
bases_pairs = ["AC", "AG", "AA", "GA", "GG", "GT", "CA", "CC", "CT", "TC", "TT", "TG"]

if __name__ == '__main__':

    if len(sys.argv) > 1:
        print("Program only accepts two total command line arguments")
        print("usage: python compiler.py")
        sys.exit(0)

    call(['clear'])
    title()
    details()
    help()

    # number of samples
    with open("dna_samples/catalog.txt") as f:
        samples = int(sum(1 for _ in f) / 7)

    # Command prompt loop
    while True:
        command = input("[DNA]>> ").strip().split()
        if len(command)  < 1:
            print("Invalid command: enter help for command menu")
        elif command[0] == "help" or command[0] == "h":
            help()
        elif command[0] == "clear" or command[0] == "c":
            call(['clear'])
        elif command[0] == "multiple" or command[0] == "m":
            multiple()
        elif command[0] == "difference" or command[0] == "d":
            difference()
        elif command[0] == "manual" or command[0] == "man":
            manual()
        elif command[0] == "test" or command[0] == "t":
            if len(command) > 1:
                if command[1] == "all":
                    test_all(samples)
                elif not command[1].isdigit():
                    print("2nd argument has to be an integer")
                elif int(command[1]) > 15 or int(command[1]) < 1:
                    print("2nd argument has to be an integer between 0 and 26")
                else:
                    test(command[1])
            else:
                file_name = input("select a number from 1 and 15: ").strip()
                if file_name.isdigit() and int(file_name) <  16 and int(file_name) > 0:
                    test(file_name)
                elif file_name == "all":
                    test_all(samples)
                else:
                    print("Integer between 0 and 26 required")
                    continue
        elif command[0] == "show" or command[0] == "s":
            if len(command) > 1:
                if command[1] == "all":
                    test_all(samples)
                elif not command[1].isdigit():
                    print("2nd argument has to be an integer")
                elif int(command[1]) > 15 or int(command[1]) < 1:
                    print("2nd argument has to be an integer between 0 and 16")
                else:
                    test(command[1])
            else:
                file_name = input("select a number from 1 and 15: ").strip()
                if file_name.isdigit() and int(file_name) <  16 and int(file_name) > 0:
                    test(file_name)
                elif file_name == "all":
                    test_all(samples)
                else:
                    print("Integer between 0 and 16 required")
                    continue
        elif command[0] == "show" or command[0] == "s":
            if len(command) > 1:
                if command[1] == "all":
                    show_all(samples)
                elif not command[1].isdigit():
                    print("2nd argument has to be an integer")
                elif int(command[1]) > 15 or int(command[1]) < 1:
                    print("2nd argument has to be an integer between 0 and 16")
                else:
                    show(command[1])
            else:
                file_name = input("select a number from 1 and 15: ").strip()
                if file_name.isdigit():
                    show(file_name)
                elif file_name == "all":
                    show_all(samples)
                else:
                    print("Integer between 0 and 16 required")
                    continue
        elif command[0] == "exit" or command[0] == "e" or command[0] == "quit" or command[0] == "q":
            sys.exit()
        elif command[0] == "details" or command[0] == "d":
            details()
        else:
            print("Invalid command: enter help for command menu")
