from sys import argv, exit
import glob
import numpy as np

if __name__ == "__main__":
    if len(argv) < 2:
        ext = '*.fit'
    else:
        ext = argv[1]
    files = np.sort(glob.glob(ext))
    with open("names.txt", 'wt') as f:
        for i in files:
            print(f"{i}\t;\t", file=f)
        print(f"File \'names.txt\' with {len(files)} entries has been created")
    exit(0)
