from sys import argv, exit
import glob
import numpy as np
from astropy.io import fits

if __name__ == "__main__":
    if len(argv) < 2:
        ext = '*.fit'
    else:
        ext = argv[1]
    files = np.sort(glob.glob(ext))
    with open("names.txt", 'wt') as f:
        for i in files:
            with fits.open(i) as hdu:
                hdr = hdu[0].header
                if 'OBJNAME' in hdr:
                    objname = hdr['OBJNAME']
                elif 'OBJECT' in hdr:
                    objname = hdr['OBJECT']
                else:
                    objname = ''
            print(f"{i}\t;\t{objname}", file=f)
        print(f"File \'names.txt\' with {len(files)} entries has been created")
    exit(0)
