# PyYAP
Yet Another Pipeline in Python for echelle data reduction

A further development of PySEX pipeline, originally written by Vadim Krushinsky (Ural Federal University, Ekaterinburg, Russia) for reduction of simple echelle spectra. Compare to the original, this version contains tonnes of improvement, many procedures were written from scratch. For correct using, the pipeline requires Python v3.8+ and the following modules:
- internal (sys, os, warnings, math, shutil, datetime, types, datetime, types, time, copy, logging).
- external (astropy, numpy, matplotlib, PyAstronomy, scipy, scikit-image, sklearn, astroquery, argparse).

Current version of pipeline is designed to reduce echelle spectra coming from the Medium Resolution Echelle Spectrograph (MRES) of the 2.4-m Thai National Telescope. The package was tested on a machine running Debian Linux Sid.
