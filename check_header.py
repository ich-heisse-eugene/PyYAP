import numpy as np
from sys import argv, exit
from astroquery.simbad import Simbad
import time
from astropy.io import fits
import astropy.time as Time
import astropy.units as u
from astropy import coordinates as coord, units as u
from PyAstronomy.pyasl import helcorr, baryCorr

def fill_headers(file_names, device):
    if device == 'mres':
        obsname = 'TNO'     # Thai National Observatory, Doi Inthanon
        obslat = 18.573828  # Latitude of the observatory
        obslon = 98.4817485 # Longitude of the observatory, E
        obsalt = 2549.      # Altitude of the observatory
        gain = 0.55            # Electronic gain in e-/ADU
        rdnoise = 2.0       # CCD readout noise
    elif device == 'eshel_ccs':
        obsname = 'CCS'     # NARIT provincial observatory, Chachoengsao
        obslat = 13.593682  # Latitude of the observatory
        obslon = 101.256209 # Longitude of the observatory, E
        obsalt = 11.        # Altitude of the observatory
        gain = 0.95            # Electronic gain in e-/ADU
        rdnoise = 7.0       # CCD readout noise
    elif device == 'eshel_krt':
        obsname = 'KRT'     # NARIT provincial observatory, Korat
        obslat = 14.873460  # Latitude of the observatory
        obslon = 102.02883  # Longitude of the observatory, E
        obsalt = 245.       # Altitude of the observatory
        gain = 0.55            # Electronic gain in e-/ADU
        rdnoise = 7.0       # CCD readout noise
    elif device == 'eshel_tno':
        obsname = 'TNO'     # Thai National Observatory, Doi Inthanon
        obslat = 18.573828  # Latitude of the observatory
        obslon = 98.4817485 # Longitude of the observatory, E
        obsalt = 2549.        # Altitude of the observatory
        gain = 0.95         # Electronic gain in e-/ADU
        rdnoise = 7.0       # CCD readout noise
    elif device == 'eshel2':
        obsname = 'TNO'     # Thai National Observatory, Doi Inthanon
        obslat = 18.573828  # Latitude of the observatory
        obslon = 98.4817485 # Longitude of the observatory, E
        obsalt = 2549.        # Altitude of the observatory
        gain = 0.55          # Electronic gain in e-/ADU
        rdnoise = 7.0       # CCD readout noise
    elif device == 'maestro':
        obsname = 'Terskol' # Terskol Observatory, Mt. Elbrus, Russia
        obslat = 43.272777  # Latitude of the observatory
        obslon = 42.500000  # Longitude of the observatory, E
        obsalt = 3100.      # Altitude of the observatory
        gain = 1.0          # Electronic gain in e-/ADU
        rdnoise = 4.0       # CCD readout noise

    files, objnames = np.loadtxt(file_names, unpack=True, usecols=(0,1), dtype=str, delimiter=';')
    simbad_session = Simbad()
    simbad_session.add_votable_fields('ra(H;ICRS;J2000)', 'dec(D;ICRS;J2000)')
    simbad_session.remove_votable_fields('coordinates')
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
            if 'DATE-OBS' in hdr:
                tm_start = Time.Time(hdr['DATE-OBS'])
            elif 'FRAME' in hdr:
                tm_start = Time.Time(hdr['FRAME'])
            if 'EXPOSURE' in hdr:
                texp = hdr['EXPOSURE'] * u.s
                hdr.set('EXPTIME', hdr['EXPOSURE'], 'Exposure (s)')
            else:
                texp = hdr['EXPTIME'] * u.s
            tm_mid = tm_start + texp/2.
            tm_end = tm_start + texp
            ut = tm_mid.ymdhms[3] + tm_mid.ymdhms[4]/60. + tm_mid.ymdhms[5]/3600.
            ut_end = tm_end.ymdhms[3] + tm_end.ymdhms[4]/60. + tm_end.ymdhms[5]/3600.
            if 'DATE' not in hdr:
                hdr.set('DATE', hdr['DATE-OBS'].split(".")[0], 'Copy of DATE-OBS')
            hdr.set('DISPAXIS', 1, 'Keyword for IRAF')
            hdr.set('GAIN', gain, '')
            hdr.set('RDNOISE', rdnoise, '')
            hdr.set('OBSGEO-B', obslat, 'Latitude of the observatory')
            hdr.set('OBSGEO-L', obslon, 'Longitude of the observatory')
            hdr.set('OBSGEO-H', obsalt, 'Altitude of the observatory')
            hdr.set('OBSERVAT', obsname, 'Thai National Observatory')

            if (objnames[ii].lower() == "flat"):
                hdr.set('IMAGETYP', 'FLAT', '')
            elif (objnames[ii].lower() == "bias"):
                hdr.set('IMAGETYP', 'BIAS', '')
            elif (objnames[ii].lower() == "thar"):
                hdr.set('IMAGETYP', 'THAR', '')
            elif (objnames[ii].lower() == "sky"):
                hdr.set('IMAGETYP', 'OBJ', '')
            elif (objnames[ii].lower() == "dark"):
                    hdr.set('IMAGETYP', 'DARK', '')

            if (objnames[ii].lower() != "flat") and (objnames[ii].lower() != "thar") \
               and (objnames[ii].lower() != "bias") and (objnames[ii].lower() != "sky") and \
               (objnames[ii].lower() != "dark"):
               query_result = simbad_session.query_object(objnames[ii])
               ra = query_result['RA_H_ICRS_J2000'][0]
               dec = query_result['DEC_D_ICRS_J2000'][0]
               if 'RA' in hdr:
                   hdr['RA'] = ra
               else:
                   hdr.set('RA', ra, 'RA in hours')
               if 'DEC' in hdr:
                   hdr['DEC'] = dec
               else:
                   hdr.set('DEC', dec, 'DEC in degrees')
               if 'EPOCH' in hdr:
                   hdr['EPOCH'] = 2000.
               else:
                   hdr.set('EPOCH', 2000., 'EPOCH of coordinates')
               star = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
               observat = coord.EarthLocation.from_geodetic(obslon, obslat, obsalt * u. m)
               dateobs = np.char.replace(tm_mid.fits, 'T', ' ')
               dateobs = Time.Time(dateobs, scale='utc', location=observat)
               ltt_bary = dateobs.light_travel_time(star)
               bjd = dateobs.jd + ltt_bary.value
               bcr, hjd = helcorr(observat.lon.degree, observat.lat.degree, obsalt, star.ra.degree, star.dec.degree, dateobs.jd)
               hdr.set('HJD', hjd, 'Heliocentric JD')
               hdr.set('BJD', bjd, 'Barycentric JD')
               hdr.set('BARYCORR', bcr, 'Barycentric correction')
               hdr.set('IMAGETYP', 'OBJ', '')
            if 'OBJNAME' in hdr:
                hdr['OBJNAME'] = objnames[ii]
            else:
                hdr.set('OBJNAME', objnames[ii], '')
            if 'DATE-OBS' in hdr:
                hdr['DATE-OBS'] = tm_mid.fits
            else:
                hdr.set('DATE-OBS', tm_mid.fits, '')
            if 'UT' in hdr:
                hdr['UT'] = ut
            else:
                hdr.set('UT', ut, '')
            if 'UTEND' in hdr:
                hdr['UTEND'] = ut_end
            else:
                hdr.set('UTEND', ut_end, '')
            hdu[0].header = hdr
            hdu.flush()
            print("File %s has been updated" %(files[ii].strip()))
    print("...done")
    return None
