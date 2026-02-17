#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 Julius Guse & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This is the WhereDoGGo? wrapper for predicting, trimming and recoding secondary structure for phylogenetic reconstruction.

#NOTE 1: All code was written and tested on ARM macOS. Please report any issues.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)
#2) S4PRED (https://anaconda.org/bioconda/s4pred)
    #2.1) PyTorch (https://anaconda.org/pytorch/pytorch)
#3) BMGE (https://anaconda.org/bioconda/bmge)

import argparse #Read and parse command line arguments to control paths and parameters for pipeline steps
import shutil #Copy, move, or delete files and directories (e.g., create temporary files or organize alignment outputs)(shell utilities)
import sys #Access to system functions such as script termination (sys.exit()) or argument lists
import os #Manage file paths and directory structures (e.g., input/output folders for analyses)
import tarfile #Unpack or create compressed archive files (tar.gz) to store large alignments
import subprocess #Start external programs or tools from within the script

#Absolut path of the script doggo_wag.py
script_path = os.path.dirname(os.path.abspath(__file__))

##Check if required non-standard libraries are installed.
import importlib.util

nonstandardlibraries = {"Bio" : "https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

##Check if required external programs are installed.
externalprograms = {"bmge":
                        {"check": ["which", "bmge"], "link": "https://anaconda.org/bioconda/bmge"},
                    "pytorch":
                        {"check": ["python", "-c", "import torch; print(torch.__version__)"], "link": "https://anaconda.org/pytorch/pytorch"}
                    }
for extprg, data in externalprograms.items():
    try:
        output = subprocess.check_output(data["check"], universal_newlines=True).strip()
    except subprocess.CalledProcessError:
        print('External program', extprg, 'not installed. Download it from: ', data["link"],'. Exiting.')
        sys.exit(1)

##Checkpoint if s4pred directory and requirements exists
s4pred_dir = os.path.join(script_path, 's4pred')
if not os.path.isdir(s4pred_dir):
    print('Directory s4pred not found in script directory. Exiting.')
    sys.exit(1)
s4pred_req = {"run_model.py": "file",
              "network.py":   "file",
              "utilities.py": "file",
              "weights":      "dir"}
for name, kind in s4pred_req.items():
    if kind == "file" or kind == "dir":
        pass
    else:
        print('The s4pred directory exists but is missing one or more required files.')
        print('Required: run_model.py, network.py, utilities.py, weights')
        print('Missing: ', name)
        sys.exit(1)

##Check if internal scripts are in the PATH.
try:
    transferindels = (subprocess.check_output("which transferindels.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script transferindels.py not found in PATH. Exiting.')
    sys.exit(1)



print(r"""
                           ▄              ▄
                          ▌▒█           ▄▀▒▌
                          ▌▒▒█        ▄▀▒▒▒▐
                         ▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐
                       ▄▄▀▒░▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐
                     ▄▀▒▒▒░░░▒▒▒░░░▒▒▒▀██▀▒▌
                    ▐▒▒▒▄▄▒▒▒▒░░░▒▒▒▒▒▒▒▀▄▒▒▌
                    ▌░░▌█▀▒▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐
                   ▐░░░▒▒▒▒▒▒▒▒▌██▀▒▒░░░▒▒▒▀▄▌
                   ▌░▒▄██▄▒▒▒▒▒▒▒▒▒░░░░░░▒▒▒▒▌
                  ▌▒▀▐▄█▄█▌▄░▀▒▒░░░░░░░░░░▒▒▒▐
                  ▐▒▒▐▀▐▀▒░▄▄▒▄▒▒▒▒▒▒░▒░▒░▒▒▒▒▌
                  ▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒▒▒░▒░▒░▒▒▐
                   ▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒░▒░▒░▒░▒▒▒▌
                   ▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▒▄▒▒▐
                    ▀▄▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▄▒▒▒▒▌
                      ▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀
                        ▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀
                           ▒▒▒▒▒▒▒▒▒▒▀▀
__        ___                   ____         ____  ____      ___
\ \      / / |__   ___ _ __ ___|  _ \  ___  / ___|/ ___| ___|__ \
 \ \ /\ / /| '_ \ / _ \ '__/ _ \ | | |/ _ \| |  _| |  _ / _ \ / /
  \ V  V / | | | |  __/ | |  __/ |_| | (_) | |_| | |_| | (_) |_|
   \_/\_/  |_| |_|\___|_|  \___|____/ \___/ \____|\____|\___/(_)
                    """)

parser = argparse.ArgumentParser(description="Henlo, am new doggo, i wag fram happy to help you!")
parser.add_argument("-p", "--prediction", required=True, choices=["in_pred","out_pred"], default="in_pred", help="Secondary structure prediction within the script = 'in_pred' (default) or extern = 'out_pred'. If 'out_pred': A directory containing the predictions must be provided as INPUT. (required)")
parser.add_argument("-sfx", "--pred_suffix", required=False, default=".s4predout")
parser.add_argument("-i", "--input", required=True, help="INPUT must be a directory containing both *_einsi.tar.gz and removemultiples.tar.gz files from doggo_sniff. (required)")
parser.add_argument("-bp", "--bmge_params", required=False, default="-t AA", help="Additional parameters to pass to BMGE trimming. Default: '-t AA' (optional)")
parser.add_argument("--timeout", type=int, required=False, default=3600)
parser.add_argument("--keep_temp", required=False, action="store_true")
parser.add_argument("--log_level", required=False, default="INFO", choices=["DEBUG","INFO","WARN","ERROR"])
parser.add_argument("-o", "--output", required=False, default="wag_output", help="OUTPUT directory where all processed results, concatenations, and archives will be stored. Default: directory 'wag_output' will be created (optional)")
args = parser.parse_args()

print('Henlo, am new doggo. I wag happy for structure nao. I speak info messages in hooman lingo.' + '\n')

#if internal prediction is wanted
internal_pred = (args.prediction == "in_pred")

##Checkpoint for input dataset existence.
if not os.path.exists(args.input) or not os.path.isdir(args.input):
    print('Input directory not found or is not a directory. Exiting.')
    sys.exit(1)
else:
    args.input = os.path.abspath(args.input)

##Checkpoint if input directory contains 'einsi.tar.gz' (always) and 'removemultiples.tar.gz' (for internal prediction)
t_archives = os.listdir(args.input)
if not any(t.endswith("einsi.tar.gz") for t in t_archives):
    print('No einsi.tar.gz found in the input directory. Exiting.')
    sys.exit(1)
if internal_pred:
    if not any(t.endswith("removemultiples.tar.gz") for t in t_archives):
        print('No removemultiples.tar.gz files found in the input directory. Exiting.')
        sys.exit(1)

##Checkpoint for output directory
if args.output == parser.get_default("output"):
    print('/////','No output directory provided, default directory "wag_output" will be created in the input path..')
    args.output = os.path.abspath(os.path.join(args.input, str(args.output)))
    try:
        os.makedirs(args.output, exist_ok=True)
    except OSError as e:
        print('Failed to create output directory:', str(e), 'Exiting.')
        sys.exit(1)
else:
    print('/////','Output will be stored in', args.output)

#Remove any previous output directories with the same name as new output directories
print ('/////','Removing directories with names identical to the new output directories..')
old_dirs = ["predout", "parsed_pred", "indeltransfer", "preconcatenation", "trimming", "concatenation", "temp", "logs"]
for d in old_dirs:
    d_path = os.path.join(args.output, d)
    if os.path.exists(d_path):
        try:
            shutil.rmtree(d_path)
        except OSError as e:
            print('Could not remove', d,':', e,'. Exiting.')

#Create new directories inside the output directory for every upcoming output
predout_dir = os.path.join(args.output, "predout")
os.makedirs(predout_dir, exist_ok=True)
indeltrans_dir = os.path.join(args.output, "indeltransfer")
os.makedirs(indeltrans_dir, exist_ok=True)
trimmed_dir = os.path.join(args.output, "trimming")
os.makedirs(trimmed_dir, exist_ok=True)
precon_dir = os.path.join(args.output, "preconcatenation")
os.makedirs(precon_dir, exist_ok=True)
concat_dir = os.path.join(args.output, "concatenation")
os.makedirs(concat_dir, exist_ok=True)
log_dir = os.path.join(args.output, "logs")
os.makedirs(log_dir, exist_ok=True)

#Create temporary directory inside the output directory for files from decompressed .tar.gz's
temp_dir = os.path.abspath(os.path.join(args.output, "temp"))
os.makedirs(temp_dir, exist_ok=True)
print('/////', 'Created output directories and temporary directory..')

##Checkpoint for storage space before extracting tar.gz files
if shutil.disk_usage(temp_dir).free < 1024**3: # 1GB minimum
    print('/////', 'Warning: Less than 1GB free space available.')
    sys.exit(1)

#Creating common-log-file
common_log = os.path.join(log_dir, "common_log_file.log")

#Decompress and save tar.gz files
if internal_pred:
    print('/////', 'Extracting alignments and unaligned files..')
    removemultiples_tar = [f for f in os.listdir(args.input) if f.endswith("removemultiples.tar.gz")]
    einsi_tar = [f for f in os.listdir(args.input) if f.endswith("einsi.tar.gz")]
    for ext in removemultiples_tar + einsi_tar:
        archive = os.path.join(args.input, ext)
        try:
            with tarfile.open(archive, "r:gz") as tar:
                for member in tar.getmembers():
                    extracted_path = os.path.abspath(os.path.join(temp_dir, member.name))
                    safety_check = os.path.commonpath([temp_dir, extracted_path]) != temp_dir
                    if safety_check:
                        print('Unsafe member in archive:', member.name, 'Exiting.')
                        sys.exit(1)
                tar.extractall(path=temp_dir)
        except (tarfile.TarError, OSError) as e:
            print('Failed to extract', archive + ':', str(e), 'Exiting.')
            sys.exit(1)
if not internal_pred:
    einsi_tar = [f for f in os.listdir(args.input) if f.endswith("einsi.tar.gz")]
    for ext in einsi_tar:
        archive = os.path.join(args.input, ext)
        try:
            with tarfile.open(archive, "r:gz") as tar:
                for member in tar.getmembers():
                    extracted_path = os.path.abspath(os.path.join(temp_dir, member.name))
                    safety_check = os.path.commonpath([temp_dir, extracted_path]) != temp_dir
                    if safety_check:
                        print('Unsafe member in archive:', member.name, 'Exiting.')
                        sys.exit(1)
                tar.extractall(path=temp_dir)
        except (tarfile.TarError, OSError) as e:
            print('Failed to extract', archive + ':', str(e), 'Exiting.')
            sys.exit(1)

#S4PRED preparation and prediction
#Create directory for demultiplied and einsi files (only one directory) and collect filenames for s4pred prediction
if internal_pred:
    pot_ua_dirs  = [os.path.join(temp_dir, d) for d in os.listdir(temp_dir) if d.endswith("removemultiples") and os.path.isdir(os.path.join(temp_dir,d))]
    pot_al_dirs = [os.path.join(temp_dir, d) for d in os.listdir(temp_dir) if d.endswith("einsi") and os.path.isdir(os.path.join(temp_dir,d))]
    
    if len(pot_ua_dirs) != 1:
        print('Expected exactly one removemultiples folder in', temp_dir + '. Exiting.')
        sys.exit(1)
    if len(pot_al_dirs) != 1:
        print('Expected exactly one einsi folder in', temp_dir + '. Exiting.')
        sys.exit(1)
    
    uafiles_dir = pot_ua_dirs[0] #takes only the first element of the list, in case there are more than 1 directory
    alfiles_dir = pot_al_dirs[0] #takes only the first element of the list, in case there are more than 1 directory
    
    ua_files = sorted([f for f in os.listdir(uafiles_dir) if f.endswith(".faademultiplied")])
    al_files = sorted([f for f in os.listdir(alfiles_dir) if f.endswith(".einsi")])
    
    if not ua_files:
        print('No unaligned (demultiplied) files found. Exiting.')
        sys.exit(1)
    if not al_files:
        print('No alignment (einsi) files found. Exiting.')
        sys.exit(1)
    
    if ua_files and al_files:
        print('/////', 'Collected', len(ua_files),  'unaligned files.')
        print('/////', 'Collected', len(al_files),  'aligned files.')
    
    #Path to the run_model.py from the s4pred-env
    run_model_path = os.path.join(s4pred_dir, "run_model.py")
    
    #Prediciting secondary structure
    print('/////', 'Running s4pred predictions..')
    
    for uafile in ua_files[:3]:
        uafile_path = os.path.join(uafiles_dir, uafile)
        predout_path = os.path.join(predout_dir, uafile.replace('.faademultiplied', '.s4predout'))
        pred_cmd = [sys.executable, run_model_path, '-t', 'fas', '-c', uafile_path] #s4pred command
    
        with open(predout_path, 'w') as pred_fh, open(common_log, 'a') as log_fh:
            try:
                result = subprocess.run(pred_cmd, stdout=pred_fh, stderr=log_fh, check=True, timeout=args.timeout)
            except subprocess.CalledProcessError as e:
                print('s4pred prediction failed: ', e)
                sys.exit(1)
    print('/////', 'Predicted and stored to secondary structure to', predout_dir)
    
#Parsing the '.s4predout' files to extract only the accession and secondary structure
print('/////', 'Parsing s4pred predictions..')
if not internal_pred:
    for file in os.listdir(args.input):
        src = os.path.join(args.input, file)
        dst = os.path.join(predout_dir, file)
        if os.path.isfile(src):
            shutil.copyfile(src, dst)

parsed_dir = os.path.join(args.output, "parsed_pred")
os.makedirs(parsed_dir, exist_ok=True)
pred_files = sorted([f for f in os.listdir(predout_dir) if f.endswith(args.pred_suffix)])
for predfile in pred_files:
    if not predfile:
        print('Found nothing to parse. Exiting.')
        sys.exit(1)
    pars_infile = os.path.join(predout_dir, predfile)
    pars_outfile = os.path.join(parsed_dir, predfile.replace(args.pred_suffix, '.faapred'))

    with open(pars_infile, "r") as fh, open(pars_outfile,"w") as out:
        for header in fh:
            if header.startswith(">"):
                seq = fh.readline() #read and ignore
                ss = fh.readline() #read and write
                conf1 = fh.readline() #read and ignore
                conf2 = fh.readline() #read and ignore
                conf3 = fh.readline() #read and ignore
            out.write(header)
            out.write(ss)

#Indeltransfer from einsi files to faas4pred files to create aligned secondary structure
print('/////', 'Transferring indels from aligned einsi files to predicted structures..')
#transferindels.py needs <unaligned> <unaligned_ext> <alignments> <alignments_ext> <output_ext> as input
transferindels_py = os.path.join(script_path, "transferindels.py")
pot_al_dirs = [os.path.join(temp_dir, d) for d in os.listdir(temp_dir) if d.endswith("einsi") and os.path.isdir(os.path.join(temp_dir,d))]
alfiles_dir = pot_al_dirs[0]
al_files = sorted([f for f in os.listdir(alfiles_dir) if f.endswith(".einsi")])
if not al_files:
    print('No alignment (einsi) files found.')
    sys.exit(1)
faas4_files = [f for f in os.listdir(parsed_dir) if f.endswith('.faapred')]
if not faas4_files:
    print('No ".faas4pred files" found. Exiting.')
    sys.exit(1)
for faas4file in faas4_files:
    transfer_outpath = os.path.join(indeltrans_dir, faas4file.replace(".faapred", '.alnfaapred'))
    transfer_cmd = [sys.executable, transferindels_py, parsed_dir, '.faapred', alfiles_dir, '.einsi', '.alnfaapred']

    with open(transfer_outpath, "w") as transfer_fh, open(common_log, "a") as log_fh:
        try:
            result = subprocess.run(transfer_cmd, cwd=indeltrans_dir, stdout=transfer_fh, stderr=log_fh, check=True, timeout=args.timeout)
        except subprocess.CalledProcessError:
            print('Indel transfer failed for', faas4file,'. See', common_log,'..')
            sys.exit(1)

#BMGE gap-only trimming of the aligned secondary structure sequences
print('/////', 'Using BMGE to trim alignments for gaps only..')
aligned_secstr_files = [f for f in sorted(os.listdir(indeltrans_dir)) if f.endswith('.alnfaapred')]
if not aligned_secstr_files:
    print('No ".alnfaapred files" found. Exiting.')
    sys.exit(1)
trimm_log = os.path.join(trimmed_dir, "trimming.log")
for aligned_secstr_file in aligned_secstr_files:
    aligned_secstr_path = os.path.join(indeltrans_dir, aligned_secstr_file)
    trimmed_outpath = os.path.join(trimmed_dir, aligned_secstr_file.replace('.alnfaapred', '.alnfaapred.trimmed'))
    trimm_cmd = ['bmge', '-i', aligned_secstr_path, '-t', 'AA', '-h', '1', '-w', '1', '-of', trimmed_outpath]

    with open(common_log, "a") as log_fh, open(trimm_log, "a") as trim_fh:
        try:
            result = subprocess.run(trimm_cmd, cwd=trimmed_dir, stdout=trim_fh, stderr=log_fh, check=True, timeout=args.timeout)
        except subprocess.CalledProcessError:
            print('Error arose while running BMGE... Exiting.')
            sys.exit(1)

#This is for running the preconcatenation script (responsible for generating a file with the dataset names to be concatenated in the next step)
print ('/////', 'Running preconcatenation via preconcatenation.sh!')
preconcat_sh = os.path.join(script_path, "preconcatenation.sh")
trimmed_files = [f for f in sorted(os.listdir(trimmed_dir)) if f.endswith('.alnfaapred.trimmed')]
if not trimmed_files:
    print('No ".trimmed files" found. Exiting.')
    sys.exit(1)
for trimmed_file in trimmed_files:
    preconcat_outpath = os.path.join(precon_dir, trimmed_file.replace('.alnfaapred.trimmed', '.preconcatenated'))
    precon_cmd = ["bash", preconcat_sh, trimmed_dir, '.trimmed']

    with open(preconcat_outpath, "w") as precon_fh, open(common_log, "a") as log_fh:
        try:
            result = subprocess.run(precon_cmd, cwd=precon_dir, stdout=precon_fh, stderr=log_fh, check=True, timeout=args.timeout)
        except subprocess.CalledProcessError:
            print('Error during preconcatenation.sh script. Exiting.')
            sys.exit(1)

##This is for running the concatenation script and creating the final .concatenation fasta file with concatenated marker genes
print ('/////', 'Concatenating secondary structure of markers via concatenation.py.')
concat_py = os.path.abspath(os.path.join(script_path, "concatenation.py"))
dataset_pass = os.path.abspath(os.path.join(precon_dir, 'dataset.pass'))
concat_file = os.path.abspath(os.path.join(concat_dir, 'concatenation.fasta'))
concat_log = os.path.abspath(os.path.join(concat_dir, 'concat.log'))
with open(concat_log, "w") as concat_fh, open(common_log, "a") as log_fh:
    try:
        result = subprocess.run([sys.executable, '-u', concat_py, trimmed_dir, dataset_pass, concat_file], cwd=concat_dir, stdout=concat_fh, stderr=log_fh, check=True, timeout=args.timeout)
    except subprocess.CalledProcessError:
        print('Error during concatenation. Exiting.')
        sys.exit(1)
print('/////', 'Parsing, Indeltransfer, Trimming and Concatenation successful!')

#Temporary directory cleanup
if not args.keep_temp:
    try:
        shutil.rmtree(temp_dir)
        print('/////', 'Temporary files and directories deleted.')
    except OSError as e:
        print('Warning: Failed to cleanup temp directory:', str(e))
else:
    print('Temporary directories and files kept!')
print('/////', 'All done!')