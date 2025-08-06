
import os
import subprocess
from configs import FASTTREE_EXE, FASTTREE_INPUT, FASTTREE_OUTPUT

def run_fasttree():
    if not os.path.exists(FASTTREE_INPUT):
        raise FileNotFoundError(f"Input file not found: {FASTTREE_INPUT}")
    cmd = [FASTTREE_EXE, "-nt", FASTTREE_INPUT]
    with open(FASTTREE_OUTPUT, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

if __name__ == "__main__":
    run_fasttree()
