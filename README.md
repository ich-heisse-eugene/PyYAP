# PyYAP
_**Py**thon-based **Y**et **A**nother **P**ipeline_ for echelle data reduction.

The idea for PyYAP was inspired by the pipeline PySEX, developed by Vadim Krushinsky (Ural Federal University, Ekaterinburg, Russia) to reduce simple echelle spectra.
PyYAP retains some legacy from the original code, yet, compared to it, this pipeline offers significant improvements, with most of the procedures written from scratch.

For correct use, the pipeline requires Python v3.10+ and the following modules:
- `sys, os, warnings, math, shutil, datetime, types, datetime, types, time, copy, logging, multiprocessing`
- `astropy, numpy, matplotlib, scipy, scikit-image, sklearn, astroquery, argparse`.

The package was independently tested on machines running Debian Linux Sid, macOS (Homebrew), and Windows 10-11.

The current version of the pipeline is primarily designed for the reduction of echelle spectra obtained with the Medium-Resolution Echelle Spectrograph (MRES) of the 2.4-m Thai
National Telescope at Doi Inthanon (Thailand). However, it was successfully adapted to other devices, including both classic slit (MAESTRO) and fibre-fed (e.g., eShel, Whoppshel) echelle spectrographs.

To install this pipeline:
1) Unpack the archive in the preferred place, and
2) Set the path to the directory with the pipeline to the variable `Pkg_path` in the first line of the file `reduce4me.py`.

Two regimes of usage:
- lazy regime (default configuration corresponds to MRES before December 2024, i.e. `--device=mres`):
`python3  /path/to/PyYAP/reduce4me.py /current/work/directory/with/FITSfiles`
- advanced regime (extract data using APEX algorithm with the subtraction of scattered and without removing temporary files):
`python3  /path/to/PyYAP/reduce4me.py /current/work/directory/with/FITSfiles --device=mres --method=APEX --sl=True --strip=False`

To get help, run this:
`python3  /path/to/PyYAP/reduce4me.py /current/work/directory/with/FITSfiles  --help`

WARNING:
1. By default, the package is configured to handle the spectra from the original version of MRES used before 16 December 2024. The upgraded version of MRES (with an Andor iKon-M CCD) uses the `umres` code.
2. Before the first run, it is mandatory to create in the directory with FITS-files an ASCII file `names.txt` of the following format:

name-of-bias-file.fits   ;  Bias

name-of-flat-file.fits   ;  Flat

name-of-obj-file.fits   ;  HD 92554*

name-of-arc-file.fits   ; ThAr
*Objects' names **must** be recognisable by Simbad. By default, PyYAP corrects the wavelength scale of scientific objects to the barycentre.

The fastest way to create this file is to run the programme `create_names_txt.py` from the pipeline's root directory: `python3 /path/to/PyYAP/create_names_txt.py "*.fits"`. After that, verify the content of the second column of `names.txt`. The minimum set of types **must** include bias, fits, object, and thar. Optionally, you can use sky in case of observations of the daytime solar spectrum.

Play with the program or contact me if you need help: eugene AT narit.or.th

Eugene Semenko

Last update 17 December 2025, NARIT
