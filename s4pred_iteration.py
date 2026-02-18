#!/usr/bin/env python3

import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", required=True)
parser.add_argument("-e", "--input_ext", required=True)
parser.add_argument("-m", "--run_model_path", required=True)
args = parser.parse_args()

input_dir  = os.path.abspath(args.input_dir)
run_model  = os.path.abspath(args.run_model_path)

for fname in os.listdir(input_dir):
    if not fname.endswith(args.input_ext):
        continue

    in_file  = os.path.join(input_dir, fname)

    cmd = [
        run_model,
        "-t", "fas",
        "-t2", "horiz",
        "-c", in_file,
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

print("Done.")