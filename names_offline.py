#!/usr/bin/env python3
"""
Programme for offline updating the raw FITS files obtained with MRES@TNT, Doi Inthanon, Chiang Mai, Thailand
Author: Eugene Semenko
Last update: 20 Dec 2025

README1st
---------
To run this code call it as:
$ python3 names_offline.py names_offline.txt --device=mres [or umres]

The ASCII file names_offline.txt must have the following format

code  ;  file name ;  object name  ; RA   ;  DEC   ; epoch

where code is a character symbol and corresponds to
c -- calibration i.e. ThAr spectrum [mandatory]
b -- bias frame [mandatory]
f -- flat field spectrum [mandatory]
d -- dark [optionally]
s -- science object (sky or any other scientific object) [mandatory]

IMPORTANT. The file names_offline.txt must use ';' as delimiter, lines starting from '#' are considered comments and will be ignored
Coordinates and epoch are mandatory only for s-type of files. Leave fields empty for the rest of types
All fields can contain spaces, but they will be stripped by this code. The maximum length of the filed variables including optional spaces:
code -- 10 Unicode symbols. The meaningful is the first non-space symbol.
file name -- maximum 50 Unicode symbols
object name -- maximum 20 Unicode symbols, but only 16 of them will be used by the code. I.e. 16 character for the name plus maximum 4 spaces.
RA -- Right ascension (ICRS). Maximum 20 symbols. Recommended format is hexadecimal (HH:MM:SS.s) of decimal degress (DD.DDDD)
DEC -- Declination (ICRS). maximum 20 symbols. Recommended format is hexadecimal (+/-DD:MM:SS.s) of decimal degress (DD.DDDD)
EPOCH -- Epoch of the coordinates. maximum 20 symbols. Recommended format is decimal year with fractions, e.g., 2000.0 or 2025.45

For example:
s  ;  polaris_spec.fits  ; Polaris  ; 02:31:49.09  ; +89:15:50.79 ; 2000.0
c  ;  thar_spec.fits     ; ThAr     ;    ;   ;
b  ;  bias1.fits         ; Bias     ;    ;   ;
...
"""
import argparse
import numpy as np
from sys import argv, exit
import time
from astropy.io import fits
import astropy.time as Time
from astropy import coordinates as coord, units as u


if __name__ == "__main__":
    ver = "2025.12"
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file with required information", type=str, default="")
    parser.add_argument("--device", help="Device: mres|umres. Use 'mres' for observations with the Newton CCD 2048x512 pix, or 'umres' otherwise (iKon-M 1Kx1K pix)", type=str, default="umres")
    # parser.add_argument("--forceupdate", help="Force updating data", action="store_true")
    args = parser.parse_args()

    if args.device == 'mres' or args.device == 'umres':
        obsname = 'TNO'         # Thai National Observatory, Doi Inthanon
        obslat = 18.573828      # Latitude of the observatory
        obslon = 98.4817485     # Longitude of the observatory, E
        obsalt = 2549.          # Altitude of the observatory
        gain = 1.16              # Electronic gain in e-/ADU. Andor Newton
        rdnoise = 4.5           # CCD readout noise
        if args.device == 'umres':   # Andor iKon-M
            gain = 1.13         # Electronic gain in e-/ADU
            rdnoise = 3.2       # CCD readout noise

    imtype, files, objnames, ra0, dec0, epoch0 = np.loadtxt(args.input.strip(), unpack=True, \
                             usecols=(0,1,2,3,4,5), delimiter=';', comments='#', \
                             dtype={'names': ('imtype', 'filename', 'objname', 'ra', 'dec', 'epoch'), \
                             'formats': ('U10', 'U50', 'U20', 'U20', 'U20', 'U20')})
    for ii in range(len(files)):
        files[ii] = files[ii].strip()
        objnames[ii] = objnames[ii].strip()
        print(f"File: {files[ii]}\tObjname: {objnames[ii]}")
        with fits.open(files[ii].strip(), mode='update') as hdu:
            hdr = hdu[0].header
            if hdr['NAXIS'] == 2:
                data = hdu[0].data.copy()
            elif hdr['NAXIS'] == 3:
                data = hdu[0].data[0].copy()
            if data.dtype.name != 'uint32':
                print(f"Scale data from {hdu[0].data.dtype.name} to 'float32'")
                hdu[0].data = np.float32(data)
            if 'PIPELINE' not in hdr or hdr['PIPELINE'].lower().find("pyyap") == -1:
                if 'DATE-OBS' in hdr:
                    if hdr['DATE-OBS'].find('T') != -1:
                        tm_start = Time.Time(hdr['DATE-OBS'])
                    elif 'UT' in hdr:
                        tm_start = Time.Time(hdr['DATE-OBS']+'T'+hdr['UT'])
                elif 'FRAME' in hdr:
                    tm_start = Time.Time(hdr['FRAME'])
                elif 'DATE' in hdr:
                    hdr.set('DATE-OBS', hdr['DATE'], 'Copy of DATE')
                    tm_start = Time.Time(hdr['DATE'])
                if 'EXPOSURE' in hdr and 'EXPTIME' not in hdr:
                    texp = hdr['EXPOSURE'] * u.s
                    hdr.set('EXPTIME', hdr['EXPOSURE'], 'Exposure (s)')
                elif 'EXPOSURE' in hdr and 'EXPTIME' in hdr:
                    texp = hdr['EXPTIME'] * u.s
                elif 'EXPTIME' in hdr:
                    texp = hdr['EXPTIME'] * u.s
                else:
                    texp = 0 * u.s
                tm_mid = tm_start + texp/2.
                tm_end = tm_start + texp
                ut = tm_mid.ymdhms[3] + tm_mid.ymdhms[4]/60. + tm_mid.ymdhms[5]/3600.
                ut_end = tm_end.ymdhms[3] + tm_end.ymdhms[4]/60. + tm_end.ymdhms[5]/3600.
                if 'DATE' not in hdr and 'DATE-OBS' in hdr:
                    hdr.set('DATE', hdr['DATE-OBS'].split(".")[0], 'Copy of DATE-OBS')
                hdr.set('DISPAXIS', 1, 'Keyword for IRAF')
                if 'GAIN' not in hdr:
                    hdr.set('GAIN', gain, "e-/ADU")
                if 'RDNOISE' not in hdr:
                    hdr.set('RDNOISE', rdnoise, "e-")
                hdr.set('OBSGEO-B', obslat, 'Latitude of the observatory')
                hdr.set('OBSGEO-L', obslon, 'Longitude of the observatory')
                hdr.set('OBSGEO-H', obsalt, 'Altitude of the observatory')
                hdr.set('OBSERVAT', obsname, 'Observatory')

                if (imtype[ii].strip()[0].lower() == "f"):
                    hdr.set('IMAGETYP', 'FLAT', '')
                elif (imtype[ii].strip()[0].lower() == "b"):
                    hdr.set('IMAGETYP', 'BIAS', '')
                elif (imtype[ii].strip()[0].lower() == "c"):
                    hdr.set('IMAGETYP', 'THAR', '')
                elif (objnames[ii].lower() == "sky" or imtype[ii].strip()[0].lower() == "s"):
                    hdr.set('IMAGETYP', 'OBJ', '')
                elif (imtype[ii].strip()[0].lower() == "d"):
                        hdr.set('IMAGETYP', 'DARK', '')

                if (imtype[ii].strip()[0].lower() == "s" and objnames[ii].lower() != "sky"):
                    if 'RA' in hdr and 'DEC' in hdr:
                        if hdr['RA'] != '' and hdr['DEC'] != '':
                            ra = hdr['RA']; dec = hdr['DEC']
                    elif 'MOUNTRA' in hdr and 'MOUNTDE' in hdr:
                        if hdr['MOUNTRA'] != '' and hdr['MOUNTDE'] != '':
                            ra = hdr['MOUNTRA']; dec = hdr['MOUNTDE']
                            hdr.set('RA', ra, 'RA in hours')
                            hdr.set('DEC', dec, 'DEC in degrees')
                    else:
                            hdr.set('RA', ra0[ii], 'RA')
                            hdr.set('DEC', dec0[ii], 'DEC')
                    if 'EPOCH' not in hdr:
                        hdr.set('EPOCH', epoch0[ii], 'EPOCH of coordinates')
                    observat = coord.EarthLocation.from_geodetic(obslon, obslat, obsalt * u. m)
                    dateobs = np.char.replace(tm_mid.fits, 'T', ' ')
                    dateobs = Time.Time(dateobs, scale='utc', location=observat)
                    hdr.set('JD', f"{dateobs.jd:.4f}", 'Julian Date')
                    hdr.set('IMAGETYP', 'OBJ', '')
                if 'OBJNAME' not in hdr and 'OBJECT' not in hdr:
                    hdr['OBJNAME'] = objnames[ii]
                elif 'OBJNAME' not in hdr:
                    hdr['OBJNAME'] = hdr['OBJECT']
                if hdr['OBJNAME'].strip() == '':
                    hdr.set('OBJNAME', objnames[ii], '')
                if 'DATE-OBS' in hdr:
                    hdr['DATE-OBS'] = tm_mid.fits
                else:
                    hdr.set('DATE-OBS', tm_mid.fits, 'Mid-exposure time')
                if 'UT' in hdr:
                    hdr['UT'] = ut
                else:
                    hdr.set('UT', ut, '')
                if 'UTEND' in hdr:
                    hdr['UTEND'] = ut_end
                else:
                    hdr.set('UTEND', ut_end, '')
                hdr.append(('PIPELINE', f"PyYAPv{ver}", 'Pipeline version'), end=True)
                hdu[0].header = hdr
                hdu.flush()
                print("File %s has been updated" %(files[ii].strip()))
    print("...done")
    exit(0)
