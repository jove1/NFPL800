#!/usr/bin/env python3

from pathlib import Path
import sys

for fname in sys.argv[1:]:
    p_in = Path(fname)
    if p_in.stem.endswith("_problem"):
        continue

    p_out = p_in.with_stem(p_in.stem + "_problem")
    
    with (
        p_in.open() as f_in, 
        p_out.open("w") as f_out 
    ):

        NORMAL, PROBLEM, SOLUTION = range(3)
        state = NORMAL
        for n, line in enumerate(f_in, 1):
            if line.startswith("#PROBLEM"):
                state = PROBLEM
                f_out.write(line)

            elif line.startswith("#END"):
                state = NORMAL
                f_out.write(line)

            elif line.startswith("#SOLUTION"):
                state = SOLUTION

            elif state == NORMAL:
                f_out.write(line)

            elif state == SOLUTION:
                pass

            else:
                assert state == PROBLEM
                assert line.startswith("#")
                f_out.write(line[1:])
