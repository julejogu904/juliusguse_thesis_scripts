#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script recodes sequences in amino acid fasta files to different reduced amino acid alphabets.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: We chose to do all recoding to numbers, which works for IQ-TREE, even the 4-state alphabets.

#Dependencies
#NONE

import os
import sys

print('#Script: aarecode.py')
print('#Version: v20241212')
print('#Usage: python aarecode.py <input_faa> <output_faa> <alphabet>')
print('#<input_faa> must be the input protein FASTA file using the complete 20-state amino-acid alphabet. Unknown states are assumed to be X. (required)')
print('#<output_faa> must be the name of the output recoded FASTA file. (required)')
print('#<alphabet> must be the reduced amino acid alphabet. It must be D6, 9, 12, 15, or 18 (Dayhoff recoding), SR2-19 (JTT-based saturation bins), G2-19 (Grantham-based bins), or KGB2-19 (WAG-based groupings). (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 4:
    print ('Three arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

# checkpoint for input file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input file found. Proceeding.')
else:
    print('Input file not found. Exiting.')
    sys.exit(1)
#TODO: Add a check that the file is really in fasta format?

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

# Checkpoint for the reduced aa alphabet
schemes = ('D6', 'D9', 'D12', 'D15', 'D18', 'SR2', 'SR3', 'SR4', 'SR5', 'SR6', 'SR7', 'SR8', 'SR9', 'SR10', 'SR11', 'SR12', 'SR13', 'SR14', 'SR15', 'SR16', 'SR17', 'SR18', 'SR19', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'KGB2', 'KGB3', 'KGB4', 'KGB5', 'KGB6', 'KGB7', 'KGB8', 'KGB9', 'KGB10', 'KGB11', 'KGB12', 'KGB13', 'KGB14', 'KGB15', 'KGB16', 'KGB17', 'KGB18', 'KGB19', 'FYMINK', 'IVYWREL', 'IVYWRELGKP', 'SECSTR3', 'DARUMAN2L', 'DARUMAL2N','TRIAL')
if sys.argv[3] in schemes:
    print('Reduced amino acid alphabet is valid. Proceeding.')
else:
    print('Reduced amino acid alphabet is not valid. Exiting.')
    sys.exit(1)

# Define the reduced amino acid alphabets here
def recode_sequence(sequence, red_aa):
    # Define the reduced amino acid alphabets here
    translation_dicts = {
        "D6": {
            "A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
            "D": "1", "E": "1", "N": "1", "Q": "1",
            "H": "2", "K": "2", "R": "2",
            "I": "3", "L": "3", "M": "3", "V": "3",
            "F": "4", "W": "4", "Y": "4",
            "C": "5"
        },
        "D9": {
            "D": "0", "E": "0", "H": "0", "N": "0", "Q": "0",
            "I": "1", "L": "1", "M": "1", "V": "1",
            "F": "2", "Y": "2",
            "A": "3", "S": "3", "T": "3",
            "K": "4", "R": "4",
            "G": "5",
            "P": "6",
            "C": "7",
            "W": "8"
        },
        "D12": {
            "D": "0", "E": "0", "Q": "0",
            "M": "1", "L": "1", "I": "1", "V": "1",
            "F": "2", "Y": "2",
            "K": "3", "H": "3", "R": "3",
            "G": "4",
            "A": "5",
            "P": "6",
            "S": "7",
            "T": "8",
            "N": "9",
            "W": "A",
            "C": "B"
        },
        "D15": {
            "D": "0", "E": "0", "Q": "0",
            "M": "1", "L": "1",
            "I": "2", "V": "2",
            "F": "3", "Y": "3",
            "G": "4",
            "A": "5",
            "P": "6",
            "S": "7",
            "T": "8",
            "N": "9",
            "K": "A",
            "H": "B",
            "R": "C",
            "W": "D",
            "C": "E"
        }, #TODO: There's an equally scoring D15 scheme that we could use instead.
        "D18": {
            "M": "0", "L": "0",
            "F": "1", "Y": "1",
            "I": "2",
            "V": "3",
            "G": "4",
            "A": "5",
            "P": "6",
            "S": "7",
            "T": "8",
            "D": "9",
            "E": "A",
            "Q": "B",
            "N": "C",
            "H": "D",
            "K": "E",
            "R": "F",
            "W": "G",
            "C": "H"
        },
        "SR2": {
            "A": "0", "D": "0", "E": "0", "G": "0", "K": "0", "N": "0", "P": "0", "Q": "0", "R": "0", "S": "0", "T": "0",
            "C": "1", "F": "1", "H": "1", "I": "1", "L": "1", "M": "1", "V": "1", "W": "1", "Y": "1"
        },
        "SR3": {
            "A": "0", "D": "0", "E": "0", "G": "0", "N": "0", "P": "0", "S": "0", "T": "0",
            "C": "1", "H": "1", "K": "1", "Q": "1", "R": "1", "W": "1",
            "F": "2", "I": "2", "L": "2", "M": "2", "V": "2", "Y": "2"
        },
        "SR4": {
            "A": "0", "G": "0", "N": "0", "P": "0", "S": "0", "T": "0",
            "C": "1", "H": "1", "W": "1", "Y": "1",
            "D": "2", "E": "2", "K": "2", "Q": "2", "R": "2",
            "F": "3", "I": "3", "L": "3", "M": "3", "V": "3"
        },
        "SR5": {
            "A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
            "C": "1", "F": "1", "W": "1", "Y": "1",
            "D": "2", "E": "2", "N": "2",
            "H": "3", "K": "3", "Q": "3", "R": "3",
            "I": "4", "L": "4", "M": "4", "V": "4"
        },
        "SR6": {
            "A": "0", "P": "0", "S": "0", "T": "0",
            "C": "1", "W": "1",
            "D": "2", "E": "2", "G": "2", "N": "2",
            "F": "3", "H": "3", "Y": "3",
            "I": "4", "L": "4", "M": "4", "V": "4",
            "K": "5", "Q": "5", "R": "5"
        },
        "SR7": {
            "A": "0", "G": "0", "S": "0", "T": "0",
            "C": "1", "W": "1",
            "D": "2", "E": "2", "N": "2",
            "F": "3", "Y": "3",
            "H": "4", "P": "4",
            "I": "5", "L": "5", "M": "5", "V": "5",
            "K": "6", "Q": "6", "R": "6"
        },
        "SR8": {
            "A": "0", "S": "0", "T": "0",
            "C": "1", "G": "1",
            "D": "2", "E": "2", "N": "2",
            "F": "3", "Y": "3",
            "H": "4", "P": "4",
            "I": "5", "L": "5", "V": "5",
            "K": "6", "Q": "6", "R": "6",
            "M": "7", "W": "7"
        },
        "SR9": {
            "A": "0", "S": "0", "T": "0",
            "C": "1", "W": "1",
            "D": "2", "E": "2",
            "F": "3", "Y": "3",
            "G": "4", "N": "4",
            "H": "5", "Q": "5",
            "I": "6", "L": "6", "V": "6",
            "K": "7", "R": "7",
            "M": "8", "P": "8"
        },
        "SR10": {
            "A": "0", "S": "0", "T": "0",
            "C": "1", "W": "1",
            "D": "2", "E": "2",
            "F": "3", "Y": "3",
            "G": "4", "N": "4",
            "H": "5", "Q": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8", "M": "8",
            "P": "9"
        },
        "SR11": {
            "A": "0", "S": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3", "Y": "3",
            "G": "4", "N": "4",
            "H": "5", "Q": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8", "M": "8",
            "P": "9",
            "W": "A"
        },
        "SR12": {
            "A": "0", "S": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3", "Y": "3",
            "G": "4",
            "H": "5", "Q": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8", "M": "8",
            "N": "9",
            "P": "A",
            "W": "B"
        },
        "SR13": {
            "A": "0", "S": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3", "Y": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8", "M": "8",
            "N": "9",
            "P": "A",
            "Q": "B",
            "W": "C"
        },
        "SR14": {
            "A": "0", "S": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3", "L": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "M": "8",
            "N": "9",
            "P": "A",
            "Q": "B",
            "W": "C",
            "Y": "D"
        },
        "SR15": {
            "A": "0", "S": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8",
            "M": "9",
            "N": "A",
            "P": "B",
            "Q": "C",
            "W": "D",
            "Y": "E"
        },
        "SR16": {
            "A": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7", "R": "7",
            "L": "8",
            "M": "9",
            "N": "A",
            "P": "B",
            "Q": "C",
            "S": "D",
            "W": "E",
            "Y": "F"
        },
        "SR17": {
            "A": "0", "T": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7",
            "L": "8",
            "M": "9",
            "N": "A",
            "P": "B",
            "Q": "C",
            "R": "D",
            "S": "E",
            "W": "F",
            "Y": "G"
        },
        "SR18": {
            "A": "0",
            "C": "1",
            "D": "2", "E": "2",
            "F": "3",
            "G": "4",
            "H": "5",
            "I": "6", "V": "6",
            "K": "7",
            "L": "8",
            "M": "9",
            "N": "A",
            "P": "B",
            "Q": "C",
            "R": "D",
            "S": "E",
            "T": "F",
            "W": "G",
            "Y": "H"
        },
        "SR19": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4",
            "G": "5",
            "H": "6",
            "I": "7", "V": "7",
            "K": "8",
            "L": "9",
            "M": "A",
            "N": "B",
            "P": "C",
            "Q": "D",
            "R": "E",
            "S": "F",
            "T": "G",
            "W": "H",
            "Y": "I"
        },
        "G2": {
            "A": "0", "C": "0", "D": "0", "E": "0", "F": "0", "G": "0", "H": "0", "I": "0", "L": "0", "M": "0", "N": "0", "P": "0", "Q": "0", "R": "0", "S": "0", "T": "0", "V": "0", "W": "0", "Y": "0",
            "K": "1"
        },
        "G3": {
            "A": "0", "C": "0", "D": "0", "F": "0", "G": "0", "M": "0", "P": "0", "Q": "0", "R": "0", "S": "0", "T": "0", "W": "0",
            "E": "1", "H": "1", "I": "1", "L": "1", "N": "1", "V": "1", "Y": "1",
            "K": "2"
        },
        "G4": {
            "A": "0", "G": "0", "P": "0", "T": "0",
            "C": "1", "D": "1", "F": "1", "M": "1", "Q": "1", "R": "1", "S": "1", "W": "1",
            "E": "2", "H": "2", "I": "2", "L": "2", "N": "2", "V": "2", "Y": "2",
            "K": "3"
        },
        "G5": {
            "A": "0", "G": "0", "P": "0", "T": "0",
            "C": "1", "D": "1", "Q": "1",
            "E": "2", "H": "2", "I": "2", "L": "2", "N": "2", "V": "2", "Y": "2",
            "F": "3", "M": "3", "R": "3", "S": "3", "W": "3",
            "K": "4"
        },
        "G6": {
            "A": "0", "G": "0",
            "C": "1", "D": "1", "Q": "1",
            "E": "2", "H": "2", "I": "2", "L": "2", "N": "2", "V": "2", "Y": "2",
            "F": "3", "M": "3", "R": "3", "S": "3", "W": "3",
            "K": "4",
            "P": "5", "T": "5"
        },
        "G7": {
            "A": "0", "G": "0",
            "C": "1", "D": "1", "Q": "1",
            "E": "2", "H": "2", "N": "2", "Y": "2",
            "F": "3", "M": "3", "R": "3", "S": "3", "W": "3",
            "I": "4", "L": "4", "V": "4",
            "K": "5",
            "P": "6", "T": "6"
        },
        "G8": {
            "A": "0", "G": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3", "H": "3", "N": "3", "Y": "3",
            "F": "4", "M": "4", "R": "4", "S": "4", "W": "4",
            "I": "5", "L": "5", "V": "5",
            "K": "6",
            "P": "7", "T": "7"
        },
        "G9": {
            "A": "0", "G": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3", "H": "3", "N": "3", "Y": "3",
            "F": "4", "M": "4", "W": "4",
            "I": "5", "L": "5", "V": "5",
            "K": "6",
            "P": "7", "T": "7",
            "R": "8", "S": "8"
        },
        "G10": {
            "A": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3", "H": "3", "N": "3", "Y": "3",
            "F": "4", "M": "4", "W": "4",
            "G": "5",
            "I": "6", "L": "6", "V": "6",
            "K": "7",
            "P": "8", "T": "8",
            "R": "9", "S": "9"
        },
        "G11": {
            "A": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3", "H": "3", "N": "3", "Y": "3",
            "F": "4", "M": "4",
            "G": "5",
            "I": "6", "L": "6", "V": "6",
            "K": "7",
            "P": "8", "T": "8",
            "R": "9", "S": "9",
            "W": "A"
        },
        "G12": {
            "A": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3", "H": "3", "N": "3", "Y": "3",
            "F": "4", "M": "4",
            "G": "5",
            "I": "6", "L": "6",
            "K": "7",
            "P": "8", "T": "8",
            "R": "9", "S": "9",
            "V": "A",
            "W": "B"
        },
        "G13": {
            "A": "0",
            "C": "1",
            "D": "2", "Q": "2",
            "E": "3",
            "F": "4", "M": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7", "L": "7",
            "K": "8",
            "P": "9", "T": "9",
            "R": "A", "S": "A",
            "V": "B",
            "W": "C"
        },
        "G14": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4", "M": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7", "L": "7",
            "K": "8",
            "P": "9", "T": "9",
            "Q": "A",
            "R": "B", "S": "B",
            "V": "C",
            "W": "D"
        },
        "G15": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4", "M": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7", "L": "7",
            "K": "8",
            "P": "9", "T": "9",
            "Q": "A",
            "R": "B",
            "S": "C",
            "V": "D",
            "W": "E"
        },
        "G16": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7", "L": "7",
            "K": "8",
            "M": "9",
            "P": "A", "T": "A",
            "Q": "B",
            "R": "C",
            "S": "D",
            "V": "E",
            "W": "F"
        },
        "G17": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7", "L": "7",
            "K": "8",
            "M": "9",
            "P": "A",
            "Q": "B",
            "R": "C",
            "S": "D",
            "T": "E",
            "V": "F",
            "W": "G"
        },
        "G18": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4",
            "G": "5",
            "H": "6", "N": "6", "Y": "6",
            "I": "7",
            "K": "8",
            "L": "9",
            "M": "A",
            "P": "B",
            "Q": "C",
            "R": "D",
            "S": "E",
            "T": "F",
            "V": "G",
            "W": "H"
        },
        "G19": {
            "A": "0",
            "C": "1",
            "D": "2",
            "E": "3",
            "F": "4",
            "G": "5",
            "H": "6", "N": "6",
            "I": "7",
            "K": "8",
            "L": "9",
            "M": "A",
            "P": "B",
            "Q": "C",
            "R": "D",
            "S": "E",
            "T": "F",
            "V": "G",
            "W": "H",
            "Y": "I"
        },
        "KGB2": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "G": "0", "P": "0", "A": "0", "C": "0", "V": "0", "I": "0", "M": "0",
            "L": "1", "F": "1", "Y": "1", "W": "1"
        },
        "KGB3": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "G": "0", "P": "0", "A": "0", "C": "0", "V": "0", "I": "0", "M": "0",
            "L": "1", "F": "1", "Y": "1",
            "W": "2"
        },
        "KGB4": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "G": "0", "P": "0", "A": "0",
            "C": "1", "I": "1", "V": "1",
            "M": "2", "L": "2", "F": "2", "Y": "2",
            "W": "3"
        },
        "KGB5": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "G": "0", "P": "0", "A": "0",
            "C": "1", "V": "1",
            "I": "2", "M": "2", "L": "2",
            "F": "3", "Y": "3",
            "W": "4"
        },
        "KGB6": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "P": "0", "A": "0",
            "G": "1",
            "C": "2", "V": "2",
            "I": "3", "M": "3", "L": "3",
            "F": "4", "Y": "4",
            "W": "5"
        },
        "KGB7": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "N": "0", "E": "0", "D": "0", "S": "0", "T": "0", "A": "0",
            "G": "1",
            "P": "2",
            "C": "3", "V": "3",
            "I": "4", "M": "4", "L": "4",
            "F": "5", "Y": "5",
            "W": "6"
        },
        "KGB8": {
            "H": "0", "R": "0", "K": "0", "Q": "0", "S": "0", "T": "0", "A": "0",
            "N": "1", "E": "1", "D": "1",
            "G": "2",
            "P": "3",
            "C": "4", "V": "4",
            "I": "5", "M": "5", "L": "5",
            "F": "6", "Y": "6",
            "W": "7"
        },
        "KGB9": {
            "H": "0", "R": "0", "K": "0", "Q": "0",
            "N": "1", "E": "1", "D": "1",
            "A": "2", "S": "2", "T": "2", "G": "2",
            "P": "3",
            "C": "4",
            "I": "5", "V": "5",
            "M": "6", "L": "6", "F": "6",
            "Y": "7",
            "W": "8"
        },
        "KGB10": {
            "R": "0", "K": "0", "H": "0", "S": "0", "A": "0",
            "Q": "1",
            "N": "2", "E": "2", "D": "2",
            "G": "3",
            "P": "4",
            "C": "5",
            "T": "6", "I": "6", "V": "6",
            "M": "7", "L": "7", "F": "7",
            "Y": "8",
            "W": "9"
        },
        "KGB11": {
            "R": "0", "K": "0", "Q": "0",
            "N": "1", "G": "1",
            "E": "2", "D": "2",
            "A": "3", "S": "3", "T": "3",
            "P": "4",
            "C": "5",
            "I": "6", "V": "6",
            "H": "7", "M": "7", "L": "7",
            "F": "8",
            "Y": "9",
            "W": "A"
        },
        "KGB12": {
            "R": "0", "K": "0", "Q": "0",
            "E": "1", "D": "1",
            "N": "2", "A": "2", "S": "2", "T": "2",
            "G": "3",
            "P": "4",
            "C": "5",
            "I": "6", "V": "6",
            "H": "7",
            "M": "8", "L": "8",
            "F": "9",
            "Y": "A",
            "W": "B"
        },
        "KGB13": {
            "R": "0", "K": "0",
            "Q": "1", "E": "1",
            "D": "2",
            "N": "3", "G": "3",
            "H": "4", "A": "4",
            "S": "5", "T": "5",
            "P": "6",
            "C": "7",
            "I": "8", "V": "8",
            "M": "9", "L": "9",
            "F": "A",
            "Y": "B",
            "W": "C"
        },
        "KGB14": {
            "R": "0",
            "K": "1",
            "Q": "2", "E": "2",
            "D": "3",
            "N": "4", "G": "4",
            "H": "5", "A": "5",
            "S": "6", "T": "6",
            "P": "7",
            "C": "8",
            "I": "9", "V": "9",
            "M": "A", "L": "A",
            "F": "B",
            "Y": "C",
            "W": "D"
        },
        "KGB15": {
            "R": "0",
            "K": "1",
            "Q": "2", "E": "2",
            "D": "3",
            "N": "4", "G": "4",
            "H": "5", "A": "5",
            "S": "6", "T": "6",
            "P": "7",
            "C": "8",
            "I": "9", "V": "9",
            "M": "A",
            "L": "B",
            "F": "C",
            "Y": "D",
            "W": "E"
        },
        "KGB16": {
            "R": "0",
            "K": "1",
            "Q": "2",
            "E": "3",
            "D": "4",
            "N": "5", "G": "5",
            "H": "6", "A": "6",
            "S": "7", "T": "7",
            "P": "8",
            "C": "9",
            "I": "A", "V": "A",
            "M": "B",
            "L": "C",
            "F": "D",
            "Y": "E",
            "W": "F"
        },
        "KGB17": {
            "R": "0",
            "K": "1",
            "Q": "2",
            "E": "3",
            "D": "4",
            "N": "5", "G": "5",
            "H": "6", "A": "6",
            "S": "7",
            "T": "8",
            "P": "9",
            "C": "A",
            "I": "B", "V": "B",
            "M": "C",
            "L": "D",
            "F": "E",
            "Y": "F",
            "W": "G"
        },
        "KGB18": {
            "R": "0",
            "K": "1",
            "Q": "2",
            "E": "3",
            "D": "4",
            "N": "5", "G": "5",
            "H": "6", "A": "6",
            "S": "7",
            "T": "8",
            "P": "9",
            "C": "A",
            "I": "B",
            "V": "C",
            "M": "D",
            "L": "E",
            "F": "F",
            "Y": "G",
            "W": "H"
        },
        "KGB19": {
            "R": "0",
            "K": "1",
            "Q": "2",
            "E": "3",
            "D": "4",
            "N": "5", "G": "5",
            "H": "6",
            "A": "7",
            "S": "8",
            "T": "9",
            "P": "A",
            "C": "B",
            "I": "C",
            "V": "D",
            "M": "E",
            "L": "F",
            "F": "G",
            "Y": "H",
            "W": "I"
        },
        "FYMINK": {
            "F": "0", "Y": "0", "M": "0", "I": "0", "N": "0", "K": "0",
            "A": "1",
            "R": "2",
            "D": "3",
            "C": "4",
            "Q": "5",
            "E": "6",
            "G": "7",
            "H": "8",
            "L": "9",
            "P": "A",
            "S": "B",
            "T": "C",
            "W": "D",
            "V": "E"
        },
        "TRIAL": {
            "A": "1", "R": "2", "N": "1", "D": "0", "C": "3", "Q": "0", "E": "0", "G": "3", "H": "2", "I": "2", "L": "2", "K": "1", "M": "2", "F": "0", "P": "0", "S": "2", "T": "0", "W": "2", "Y": "1", "V": "2"
        },
        "IVYWREL": {
            "I": "0", "V": "0", "Y": "0", "W": "0", "R": "0", "E": "0", "L": "0",
            "A": "1",
            "N": "2",
            "D": "3",
            "C": "4",
            "Q": "5",
            "G": "6",
            "H": "7",
            "K": "8",
            "M": "9",
            "F": "A",
            "P": "B",
            "S": "C",
            "T": "D"
        },
        "IVYWRELGKP": {
            "I": "0", "V": "0", "Y": "0", "W": "0", "R": "0", "E": "0", "L": "0", "G": "0", "K": "0", "P": "0",
            "A": "1",
            "N": "2",
            "D": "3",
            "C": "4",
            "Q": "5",
            "H": "6",
            "M": "7",
            "F": "8",
            "S": "9",
            "T": "A"
        },
        "SECSTR3": {
            "H": "0",
            "E": "1",
            "C": "2"
        },
        "DARUMAN2L": {
            "0": "N",
            "1": "Y"
        },
        "DARUMAL2N": {
            "N": "0",
            "Y": "1"
        }
    }
    translation_dict = translation_dicts.get(red_aa)
    new_sequence = ""
    for aa in sequence:
        if aa == 'X':
            new_sequence += '?'
        else:
            new_sequence += translation_dict.get(aa, aa)
    return new_sequence

print('Recoding sequences.')
with open(sys.argv[1], "r") as goat, open(sys.argv[2], "w") as sheep:
    for line in goat:
        if line.startswith(">"):
            sheep.write(line)  # Keep header lines
        else:
            recoded_sequence = recode_sequence(line.strip(), sys.argv[3])
            sheep.write(recoded_sequence + "\n")  # Perform recoding on sequence lines

print('All done!')
