# PyYAP
Yet Another Pipeline in Python for echelle data reduction

Further development of the PySEX pipeline, originally written by Vadim Krushinsky (Ural Federal University, Ekaterinburg, Russia) for reducing simple echelle spectra.
Compared to the original, this version contains many improvements, and many procedures were written from scratch. For correct use, the pipeline requires Python v3.10+
and the following modules:
- internal (sys, os, warnings, math, shutil, datetime, types, datetime, types, time, copy, logging, multiprocessing).
- external (astropy, numpy, matplotlib, scipy, scikit-image, sklearn, astroquery, argparse).

The current version of the pipeline is designed for the reduction of echelle spectra coming from the Medium Resolution Echelle Spectrograph (MRES) of the 2.4-m Thai
National Telescope at Doi Inthanon (Thailand). The package was independently tested on machines running Debian Linux Sid, macOS (Homebrew), and Windows 10-11.

To install the package, unpack the archive in the preferred place and set its path in the first line of the file reduce4me.py (variable Pkg_path).

Examples of usage:
- lazy regime:
python3  /path/to/the/package/reduce4me.py /current/work/directory/with/FITSfiles 
- advanced regime (extract data using APEX algorithm with the subtraction of scattered and without removing of temporary files):
python3  /path/to/the/package/reduce4me.py /current/work/directory/with/FITSfiles --device=mres --method=APEX --sl=True --strip=False

To get help run this:
python3  /path/to/the/package/reduce4me.py /current/work/directory/with/FITSfiles  --help

Warning:
1. By default, the package is configured to work with data from the original version of MRES. The upgraded version of MRES (with an Andor iKon-M CCD) has code umres.
2. Before the first run, it is mandatory to create in the directory with FITS-files an ASCII file names.txt of the following format:

name-of-bias-file.fits   ;  Bias

name-of-flat-file.fits   ;  Flat

name-of-obj-file.fits   ;  HD 92554*

name-of-arc-file.fits   ; ThAr
*Objects' names must be recognizable by Simbad

The minimum set of types must include bias, fits, object, thar

Play with the program or contact me if you need help: eugene AT narit.or.th

Eugene Semenko

Last update 27 May 2025, NARIT
