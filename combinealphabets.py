#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 MEDEA Lab
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script combines the (reduced) alphabets from two sets of unaligned or aligned sequences.

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

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

print('#Script: combinealphabets.py')
print('#Version: vYYYYMMDD')
print('#Usage: python combinealphabets.py <datasets1> <input_ext1> <states1> <datasets2> <input_ext2> <states2> <output_ext>')
print('#<datasets1> must be the directory containing the FASTA files with <input_ext1> and alphabet <states1>. (trailing slash optional) (required)')
print('#<input_ext1> must be the extension of the FASTA files in <datasets1> whose alphabet <states1> will be combined. (leading dot optional) (required)')
print('#<states1> must be an integer corresponding to the number of states in the first alphabet for FASTA files in <datasets1> with <input_ext1>. The combined number of states cannot exceed 36. Gaps must be in the same positions for all dataset pairs. Ambiguous states (must be "?") in one or both datasets are retained in the output. (required)')
print('#<datasets2> must be the directory containing the FASTA files with <input_ext2> and alphabet <states2>. (trailing slash optional) (required)')
print('#<input_ext2> must be the extension of the FASTA files in <datasets2> whose alphabet <states2> will be combined. (leading dot optional) (required)')
print('#<states2> must be an integer corresponding to the number of states in the second alphabet for FASTA files in <datasets2> with <input_ext2>. The combined number of states cannot exceed 36. Gaps must be in the same positions for all dataset pairs. Ambiguous states (must be "?") in one or both datasets are retained in the output. (required)')
print('#<output_ext> must be the extension of the FASTA files where the secondary structure-converted sequences will be written in single-line format. (leading dot optional) (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 8:
    print ('Seven arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

datasetsdir1 = sys.argv[1]
input_ext1 = sys.argv[2]
#states1 = int(sys.argv[3])
datasetsdir2 = sys.argv[4]
input_ext2 = sys.argv[5]
#states2 = int(sys.argv[6])
output_ext = sys.argv[7]

# Check if the extensions start with a dot, otherwise add them.
if input_ext1.startswith('.') == False: #This is not absolutely necessary, since os.path.splitext will detect the extension anyway. It's more of a precaution against double dots.
    input_ext1 = str('.' + input_ext1)
if input_ext2.startswith('.') == False:
    input_ext2 = str('.' + input_ext2)
if output_ext.startswith('.') == False:
    output_ext = str('.' + output_ext)

#Checkpoint for datasets directories existence and trailing slash. Convert to abspath to make sure there are no issues when called through doggo_wag.
if os.path.exists(datasetsdir1) == True:
    print ('First datasets directory found. Proceeding.')
    datasetsdir1 = os.path.abspath(datasetsdir1)
    datasetsdir1 = os.path.join(datasetsdir1, '')
else:
    print ('First datasets directory not found. Exiting.')
    sys.exit(1)
if os.path.exists(datasetsdir2) == True:
    print ('Second datasets directory found. Proceeding.')
    datasetsdir2 = os.path.abspath(datasetsdir2)
    datasetsdir2 = os.path.join(datasetsdir2, '')
else:
    print ('Second datasets directory not found. Exiting.')
    sys.exit(1)

#Check if files with a given extension exist in the datasets directory and create a list of them.
filenames1 = []
for fname1 in os.listdir(datasetsdir1):
    if fname1.endswith(input_ext1):
        fname1 = os.path.join(datasetsdir1, fname1)
        filenames1.append(fname1)
if len(filenames1) > 0:
    print('File(s) with the first input extension found in the first datasets directory. Proceeding.')
else:
    print('No files with the first input extension found in the first datasets directory. Exiting.')
    sys.exit(1)

filenames2 = []
for fname2 in os.listdir(datasetsdir2):
    if fname2.endswith(input_ext2):
        fname2 = os.path.join(datasetsdir2, fname2)
        filenames2.append(fname2)
if len(filenames2) > 0:
    print('File(s) with the second input extension found in the second datasets directory. Proceeding.')
else:
    print('No files with the second input extension found in the second datasets directory. Exiting.')
    sys.exit(1)

#Iterate over each file in the list for the first alphabet.
for iterfname1 in filenames1:
    corrfname1 = str(datasetsdir2 + os.path.basename(iterfname1).split(os.extsep, 1)[0] + input_ext2) #Find the corresponding dataset for the second alphabet, with which they share a file stem.
    if not os.path.exists(corrfname1):
        print ('Corresponding file for ' + iterfname1 + ' not found. Exiting.' )
        sys.exit(1)
# Then vice versa. By doing it reciprocally we ensure a one to one correspondence between files.
for iterfname2 in filenames2:
    corrfname2 = str(datasetsdir1 + os.path.basename(iterfname2).split(os.extsep, 1)[0] + input_ext1) #Find the corresponding dataset for the second alphabet, with which they share a file stem.
    if not os.path.exists(corrfname2):
        print ('Corresponding file for ' + iterfname2 + ' not found. Exiting.' )
        sys.exit(1)

#Check if states1 and states2 are integers between 2 and 18 (included). The reason is that we need to have reduced alphabets (so no single-state) and the combined number of states cannot exceed 36 (IQ-TREE morphological states are 0-9, then A-Z).
#We define the number of states instead of automatically detecting in case they don't all appear in all datasets e.g., a helical protein won't contain any E characters for sheet.
flag1 = True
try:
    states1 = int(sys.argv[3])
except:
    flag1 = False
if flag1 and (2 <= states1 <= 18):
    print('First number of states given is valid. Proceeding.')
else:
    print('First number of states given is invalid. Exiting.')
    sys.exit(1)
flag2 = True
try:
    states2 = int(sys.argv[6])
except:
    flag2 = False
if flag2 and (2 <= states2 <= 18):
    print('Second number of states given is valid. Proceeding.')
else:
    print('Second number of states given is invalid. Exiting.')
    sys.exit(1)

#Check if the combined number of states does not exceed 36.
combined_states = states1 * states2
if combined_states <= 36:
    print('Combined number of states is valid. Proceeding.')
else:
    print('Combined number of states is not valid. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print('Removing files with names identical to the output.')
removal = ('rm -r *' + output_ext + ' 2> /dev/null')
os.system(removal)

# Sort file lists to ensure consistent pairing
filenames1.sort()
filenames2.sort()

# Prepare IQ-TREE alphabet symbols (0-9 then A-Z)
alphabet_symbols = [str(i) for i in range(10)] + [chr(c) for c in range(ord('A'), ord('Z')+1)]

# Define valid alphabets for states1 and states2
alphabet1 = alphabet_symbols[:states1]
alphabet2 = alphabet_symbols[:states2]

# Create mapping for combined alphabet (cartesian product of states1 x states2)
combined_mapping = {}
combined_alphabet = []
counter = 0
for c1 in alphabet1:
    for c2 in alphabet2:
        combined_mapping[(c1, c2)] = alphabet_symbols[counter]
        combined_alphabet.append(alphabet_symbols[counter])
        counter += 1

# Sanity check
#if len(combined_mapping) != combined_states:
#    print("Combined alphabet mapping size does not match expected number of states. Exiting.")
#    sys.exit(1)

# Print combined alphabet mapping
mapping_strings = []
for c1 in alphabet1:
    for c2 in alphabet2:
        symbol = combined_mapping[(c1, c2)]
        mapping_strings.append(f"{c1}{c2}={symbol}")
print("Combined alphabet states are:", ", ".join(mapping_strings))

print('Combining alphabets.')
# Iterate through file pairs and combine sequences
for lf1, lf2 in zip(filenames1, filenames2):
    filestem = str(os.path.basename(lf1).split(os.extsep, 1)[0])
    records1 = list(SeqIO.parse(lf1, "fasta"))
    records2 = list(SeqIO.parse(lf2, "fasta"))

    if len(records1) != len(records2):
        print("Datasets", lf1, "and", lf2, "do not contain the same number of sequences. Exiting.")
        sys.exit(1)

    combined_records = []
    for r1, r2 in zip(records1, records2):
        if r1.id != r2.id:
            print("Sequence IDs", r1.id, "and", r2.id, "do not match. Exiting.")
            sys.exit(1)
        if len(r1.seq) != len(r2.seq):
            print("Sequence lengths for IDs", r1.id, "and", r2.id, "do not match. Exiting.")
            sys.exit(1)

        # Check sequence characters are within allowed alphabets
        allowed1 = set(alphabet1) | {"-", "?"}
        allowed2 = set(alphabet2) | {"-", "?"}

        for ch1 in set(str(r1.seq)):
            if ch1 not in allowed1:
                print(f"Character {ch1} in {r1.id} exceeds alphabet1 ({states1} states). Exiting.")
                sys.exit(1)
        for ch2 in set(str(r2.seq)):
            if ch2 not in allowed2:
                print(f"Character {ch2} in {r2.id} exceeds alphabet2 ({states2} states). Exiting.")
                sys.exit(1)

        # Combine sequences. Unmatched gaps will give an error. If either of the two sequences has an ambiguous character ("?"), the combined_seq also gets an "?" for that position.
        combined_seq = []
        for pos, (lc1, lc2) in enumerate(zip(str(r1.seq), str(r2.seq)), start=1):
            if lc1 == "-" and lc2 != "-":
                print(f"Gap mismatch at position {pos} in {r1.id}. Exiting.")
                sys.exit(1)
            if lc2 == "-" and lc1 != "-":
                print(f"Gap mismatch at position {pos} in {r1.id}. Exiting.")
                sys.exit(1)

            if lc1 == "-" and lc2 == "-":
                combined_seq.append("-")
            elif lc1 == "?" or lc2 == "?":
                combined_seq.append("?")
            else:
                combined_seq.append(combined_mapping[(lc1, lc2)])

        combined_seq_str = "".join(combined_seq)

        # Build new SeqRecord
        new_record = SeqRecord(Seq(combined_seq_str), id=r1.id, description="")
        combined_records.append(new_record)

    # Write combined file
    with open(str(filestem + output_ext), "w") as outfile:
        SeqIO.write(combined_records, outfile, "fasta")

print('All done!')
