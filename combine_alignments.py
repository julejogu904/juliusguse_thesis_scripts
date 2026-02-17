#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 Julius Guse & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#Combines files that were processed with transferindels.py to check and compare.

#NOTE 1: All code was written and tested on ARM macOS. Please report any issues.

import sys, pathlib
from collections import defaultdict

if len(sys.argv) < 4:
    sys.exit(f"Usage: {sys.argv[0]} <outdir> <dir1> <dir2> <dir3> {sys.argv[5]}")

outdir = pathlib.Path(sys.argv[1]); outdir.mkdir(parents=True, exist_ok=True)
dirs = [pathlib.Path(d) for d in sys.argv[2:5]]
outext = sys.argv[5]

def group_by_prefix(d: pathlib.Path):
    g = defaultdict(list)
    for p in d.iterdir():
        if p.is_file():
            pref = p.name.split('.', 1)[0]
            g[pref].append(p)
    return g

groups = [group_by_prefix(d) for d in dirs]

# 1) gleiche Anzahl Dateien?
counts = [sum(len(v) for v in g.values()) for g in groups]
for prefix, count in zip(dirs, counts):
    if len(set(counts)) != 1:
        print(f"Warnung: unterschiedliche Dateizahlen {counts}", file=sys.stderr)

# 2) gemeinsame Präfixe
common = set(groups[0])
for g in groups[1:]:
    common &= set(g)
if not common:
    sys.exit("Keine gemeinsamen Präfixe gefunden.")

# Mehrfachtreffer je Präfix abfangen
bad = [k for g in groups for k,v in g.items() if len(v) != 1]
if bad:
    sys.exit(f"Mehr als eine Datei pro Präfix gefunden: {sorted(set(bad))}")

def read_fasta(path: pathlib.Path):
    hdr, seq = None, []
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s: continue
            if s.startswith(">"):
                if hdr: yield hdr, "".join(seq)
                hdr, seq = s, []
            else:
                seq.append(s)
    if hdr: yield hdr, "".join(seq)

def fasta_map(path: pathlib.Path):
    return dict(read_fasta(path))

for pref in sorted(common):
    files = [g[pref][0] for g in groups]
    maps = [fasta_map(f) for f in files]
    # 3) gemeinsame Header
    headers = set(maps[0])
    for m in maps[1:]:
        headers &= set(m)
    if not headers:
        print(f"Warnung: keine gemeinsamen Header für Präfix {pref}", file=sys.stderr)
        continue

    out = outdir / f"{pref}.{outext}"
    with out.open("w") as fh:
        fh.write(f"# ===== {pref} =====\n")
        for h in sorted(headers):
            fh.write(f"{h}\n")
            for m in maps:
                fh.write(m[h] + "\n")
            fh.write("\n")