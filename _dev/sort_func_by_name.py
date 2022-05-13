#!/usr/bin/env python
import sys

lines = []
with open(sys.argv[1], "r") as fh:
    grouping = False
    last = []
    for line in fh:
        line = line.rstrip()
        if grouping:
            last.append(line)
            if not line:
                grouping = False
                lines.append(last)
        else:
            if line.startswith("def") or line.startswith("@"):
                last = [line]
                grouping = True

for line in sorted(lines, key=lambda x: x[0]):
    print(*line, sep="\n")
    print("")
