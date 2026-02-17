#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 MEDEA Lab
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script parses the multi-line combined FASTA and horizontal PSIPRED output of S4PRED to produce a single-line fasta with the secondary structure only, converting any residues that weren't predicted with sufficient confidence to an ambiguous state X.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio":"https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

print('#Script: parses4pred.py')
print('#Version: vYYYYMMDD')
print('#Usage: python parses4pred.py <datasets> <input_ext> <output_ext> <threshold>')
print('#<datasets> must be the directory containing the S4PRED combined FASTA-like and horizontal PSIPRED output files with <input_ext>. (trailing slash optional) (required)')
print('#<input_ext> must be the extension of the S4PRED files in <datasets> that will be parsed. (leading dot optional) (required)')
print('#<output_ext> must be the extension of the FASTA files where the secondary structure-converted sequences will be written in single-line format. (leading dot optional) (required)')
print('#<threshold> must be an integer between 0 and 9 (inclusive) of the minimum confidence score for a predicted secondary structure to be included. If no state probability reaches <threshold>, it is treated as ambiguous ("X"). 0 will always retain the secondary structure predicted. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 5:
    print ('Four arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

datasetsdir = sys.argv[1]
input_ext = sys.argv[2]
output_ext = sys.argv[3]
#threshold = float(sys.argv[4])

# Check if the extensions start with a dot, otherwise add them.
if input_ext.startswith('.') == False: #This is not absolutely necessary, since os.path.splitext will detect the extension anyway. It's more of a precaution against double dots.
    input_ext = str('.' + input_ext)
if output_ext.startswith('.') == False:
    output_ext = str('.' + output_ext)

#Checkpoint for datasets directory existence and trailing slash. Convert to abspath to make sure there are no issues when called through doggo_wag.
if os.path.exists(datasetsdir) == True:
    print ('Datasets directory found. Proceeding.')
    datasetsdir = os.path.abspath(datasetsdir)
    datasetsdir = os.path.join(datasetsdir, '')
else:
    print ('Datasets directory not found. Exiting.')
    sys.exit(1)

#Check if files with a given extension exist in the datasets directory and create a list of them.
filenames = []
for fname in os.listdir(datasetsdir):
    if fname.endswith(input_ext):
        fname = os.path.join(datasetsdir, fname)
        filenames.append(fname)
if len(filenames) > 0:
    print('File(s) with the input extension found in the datasets directory. Proceeding.')
else:
    print('No files with the input extension found in the datasets directory. Exiting.')
    sys.exit(1)

#Check if threshold can be an integer between 0 and 9 (included).
flag1 = True
try:
    threshold = int(sys.argv[4])
except:
    flag1 = False
if flag1 and (0 <= threshold <= 9):
    print('Threshold given is valid. Proceeding.')
else:
    print('Threshold given is invalid. Exiting.')
    sys.exit(1)
#Create a string of the threshold without the decimal point that can be added to the outfile names.
threshold_str = str(threshold).replace(".", "")

#Remove any previous output files with the same name.
print('Removing files with names identical to the output.')
removal = ('rm -r *' + output_ext + ' 2> /dev/null')
os.system(removal)

# Loop to parse each file with input_ext in datasets.
print('Parsing S4PRED output files.')
mapping = {"H": "0", "E": "1", "C": "2", "X": "?"}
for infile in filenames:
    #Separate the fasta file stem to use for the temporary single-line alignment files.
    filestem = str(os.path.basename(infile).split(os.extsep, 1)[0])
    records = []
    with open(infile, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    # Each record starts with ">" and contains:
    #   line 0: >header
    #   line 1: sequence
    #   line 2: predicted SS
    #   following lines: "Conf:" blocks until next ">"
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue

        header_line = lines[i]           # includes ">"
        header_text = header_line[1:]    # remove ">", keep everything else
        aa_seq = lines[i + 1]
        ss_pred = lines[i + 2]

        # Collect confidence digits across multiple "Conf:" lines
        j = i + 3
        conf_digits = []
        while j < len(lines) and not lines[j].startswith(">"):
            if lines[j].startswith("Conf:"):
                conf_digits.append(lines[j][6:].strip())
            j += 1
        conf_string = "".join(conf_digits)

        if len(conf_string) != len(aa_seq):
            raise ValueError(
                f"Conf length {len(conf_string)} != seq length {len(aa_seq)} in {header_text}")

        # Apply threshold
        new_pred = []
        for aa, ss, conf in zip(aa_seq, ss_pred, conf_string):
            original_symbol = mapping[ss]

            # Apply threshold
            if int(conf) >= threshold:
                new_symbol = original_symbol
            else:
                new_symbol = mapping["X"]

            new_pred.append(new_symbol)

        # Create SeqRecord
        seqrec = SeqRecord(
            Seq("".join(new_pred)),
            id=header_text, # Original FASTA is retained whole as id
            description="" # Nothing more to add to the header
        )
        records.append(seqrec)

        i = j  # move to next record
    # Write output fasta
    with open(str(filestem + output_ext + "_" + threshold_str), "w") as outfile:
        SeqIO.write(records, outfile, "fasta")

print('All done!')
