#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 MEDEA Lab
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script trims an alignment (gappy and/or parsimony informative or constant sites) while being alphabet-agnostic. Thus, it can be used on morphological character alignments or non-standard amino acids or mixed alignments.``

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that works on all files in the working directory with a given extension.

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

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

print('#Script: trimmorphological.py')
print('#Version: vYYYYMMDD')
print('#Usage: python trimmorphological.py <datasets> <input_ext> <output_ext> <data_type> <mode> <gap_threshold>')
print('#<datasets> must be the directory containing the alignments with <input_ext> to be trimmed in FASTA format. (trailing slash optional) (required)')
print('#<input_ext> must be the extension of the alignments to be trimmed. (leading dot optional) (required)')
print('#<output_ext> must be the extension of the aligned FASTA files. (leading dot optional) (required)')
print('#<data_type> must be "dna" for DNA, "aa" for amino acids, or "morpho" for morphological data. (required)')
print('#<mode> must be "gaps" (remove positions with gaps over a fraction; ambiguous characters are treated as gaps), constant (remove constant sites), "nonparsimony" (remove non-parsimony informative variable sites), "parsimony" (remove parsimony informative sites). All modes can be combined divided by dashes, as long as at least one of constant, nonparsimony, or parsimony remains. Gaps are removed first if specified with other types. (required)')
print('#<gap_threshold> must be a float of the minimum fraction of gaps in a position to be removed by a <mode> including "gaps".')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 6:
    print ('Five arguments found. Proceeding.')
elif len(sys.argv) == 7:
    print ('Six arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

datasetsdir = sys.argv[1]
input_ext = sys.argv[2]
output_ext = sys.argv[3]
data_type = sys.argv[4]
mode = sys.argv[5]
#gap_threshold = float(sys.argv[6])

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

# Define valid alphabets and ambiguous/missing data
if data_type == "dna":
    print('Data type DNA found. Proceeding.')
    valid = set("ATCG")
    unknown = set("N")
elif data_type == "aa":
    print('Data type amino acids found. Proceeding.')
    valid = set("ACDEFGHIKLMNPQRSTVWYUO") # 20 aa + U (selenocysteine) + O (pyrrolysine)
    unknown = set("X")
elif data_type == "morpho":
    print('Data type morphological found. Proceeding.')
    valid = set("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    unknown = set("?")
else:
    print('Data type is not valid. Exiting.')
    sys.exit(1)

#Checkpoint that modes (dash-delimited) are valid, and at least one type of site (constant, nonparsimony, parsimony) will remain.
valid_modes = {"gaps", "constant", "nonparsimony", "parsimony"}
modeset = set(mode.split("-"))
min_modes = {"constant", "nonparsimony", "parsimony"}
if not modeset.issubset(valid_modes):
    print("Mode given is invalid. Exiting.")
    sys.exit(1)
if min_modes.issubset(modeset):
    print("At least one of constant, nonparsimony, or parsimony must remain. Exiting.")
    sys.exit(1)
print("Mode given is valid. Proceeding.")

# If "gaps" was specified, then we need a gap threshold argument
if "gaps" in modeset:
    if len(sys.argv) < 7:
        print("Gap threshold for gaps mode not given. Exiting.")
        sys.exit(1)
    flag = True
    try:
        gap_threshold = float(sys.argv[6])
    except:
        flag = False
    if flag and (0.0 <= gap_threshold <= 1.0):
        print('Gap threshold given is valid. Proceeding.')
    else:
        print('Gap threshold given is invalid. Exiting.')
        sys.exit(1)

#Remove any previous output files with the same name.
print('Removing files with names identical to the output.')
removal = ('rm -r *' + output_ext + ' 2> /dev/null')
os.system(removal)

# For each alignment go over it twice, once for gappy sites and a second one for constant, nonparsimony, and parsimony.
for untrimmed in filenames:

    alignment = AlignIO.read(untrimmed, "fasta")
    outfile = str(os.path.basename(untrimmed).split(os.extsep, 1)[0] + output_ext)

    # Trim gaps (if the option has been picked)
    if "gaps" in modeset:
        kept_columns = []
        for col_idx in range(alignment.get_alignment_length()):
            column = alignment[:, col_idx]
            gap_count = sum(1 for c in column if c == '-' or c in unknown)
            gap_fraction = gap_count / len(column)
            if gap_fraction >= gap_threshold:
                continue

            kept_columns.append(col_idx)

        new_records = []
        for rec in alignment:
            new_seq = ''.join(rec.seq[i] for i in kept_columns)
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))
        alignment = MultipleSeqAlignment(new_records)

    # Trim constant, nonparsimony, parsimony (if any of the options has been picked)
    if "constant" in modeset or "nonparsimony" in modeset or "parsimony" in modeset:
        kept_columns = []
        for col_idx in range(alignment.get_alignment_length()):
            column = alignment[:, col_idx]

            # remove gaps and ambiguous from consideration. If only gaps and ambiguous characters are found in a site it will be removed anyway, as a backup if "gaps" is not picked.
            filtered = [c for c in column if c in valid]
            if not filtered:
                continue  # all unknowns/gaps

            # constant sites
            if "constant" in modeset:
                if len(set(filtered)) == 1: #If there's only one character in the set i.e., the site is constant, remove it.
                    continue

            # nonparsimony sites
            if "nonparsimony" in modeset:
                counts = {}
                for c in filtered:
                    counts[c] = counts.get(c, 0) + 1
                if len(counts) > 1:  # variable
                    informative_states = [s for s in counts if counts[s] >= 2]
                    if len(informative_states) < 2:
                        continue

            # parsimony sites
            if "parsimony" in modeset:
                counts = {}
                for c in filtered:
                    counts[c] = counts.get(c, 0) + 1
                informative_states = [s for s in counts if counts[s] >= 2]
                if len(informative_states) >= 2:
                    continue

            kept_columns.append(col_idx)

        # Rebuild alignment with new SeqRecords
        new_records = []
        for rec in alignment:
            new_seq = ''.join(rec.seq[i] for i in kept_columns)
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))
        alignment = MultipleSeqAlignment(new_records)

    # write output
    AlignIO.write(alignment, outfile, "fasta")

print('All done!')
