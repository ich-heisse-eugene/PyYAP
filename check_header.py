import numpy as np
from sys import argv, exit
from astroquery.simbad import Simbad
import time
from astropy.io import fits
import astropy.time as Time
from astropy import coordinates as coord, units as u

def fill_headers(file_names, device):
    if device == 'mres' or device == 'umres':
        obsname = 'TNO'         # Thai National Observatory, Doi Inthanon
        obslat = 18.573828      # Latitude of the observatory
        obslon = 98.4817485     # Longitude of the observatory, E
        obsalt = 2549.          # Altitude of the observatory
        gain = 1.0              # Electronic gain in e-/ADU. Andor Newton
        rdnoise = 3.0           # CCD readout noise
        if device == 'umres':   # Andor iKon-M
            gain = 1.13         # Electronic gain in e-/ADU
            rdnoise = 3.2       # CCD readout noise
    elif device == 'eshel_tno':
        obsname = 'TNO'         # Thai National Observatory, Doi Inthanon
        obslat = 18.573828      # Latitude of the observatory
        obslon = 98.4817485     # Longitude of the observatory, E
        obsalt = 2549.          # Altitude of the observatory
        gain = 0.95             # Electronic gain in e-/ADU
        rdnoise = 7.0           # CCD readout noise
    elif device == 'eshel_zdnc':
        obsname = 'Zdanice'     # ASA 0.8m telescope, Ždánice
        obslat = 49.06574       # Latitude of the observatory
        obslon = 17.03784       # Longitude of the observatory, E
        obsalt = 250            # Altitude of the observatory
        gain = 1                # Electronic gain in e-/ADU
        rdnoise = 1             # CCD readout noise
    ######## LEGACY instruments ############
    # elif device == 'eshel_ccs': # Legacy device. Not available anymore
    #     obsname = 'CCS'         # NARIT provincial observatory, Chachoengsao
    #     obslat = 13.593682      # Latitude of the observatory
    #     obslon = 101.256209     # Longitude of the observatory, E
    #     obsalt = 11.            # Altitude of the observatory
    #     gain = 0.95             # Electronic gain in e-/ADU
    #     rdnoise = 7.0           # CCD readout noise
    # elif device == 'eshel_krt': # Legacy device. Not available anymore
    #     obsname = 'KRT'         # NARIT provincial observatory, Korat
    #     obslat = 14.873460      # Latitude of the observatory
    #     obslon = 102.02883      # Longitude of the observatory, E
    #     obsalt = 245.           # Altitude of the observatory
    #     gain = 0.95             # Electronic gain in e-/ADU
    #     rdnoise = 7.0           # CCD readout noise
    # elif device == 'maestro':
    #     obsname = 'Terskol'     # Terskol Observatory, Mt. Elbrus, Russia
    #     obslat = 43.272777      # Latitude of the observatory
    #     obslon = 42.500000      # Longitude of the observatory, E
    #     obsalt = 3100.          # Altitude of the observatory
    #     gain = 1.0              # Electronic gain in e-/ADU
    #     rdnoise = 4.0           # CCD readout noise

    files, objnames = np.loadtxt(file_names, unpack=True, usecols=(0,1), dtype=str, delimiter=';')
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
                if hdr['DATE-OBS'].find('T') != -1:
                    tm_start = Time.Time(hdr['DATE-OBS'])
                elif 'UT' in hdr:
                    tm_start = Time.Time(hdr['DATE-OBS']+'T'+hdr['UT'])
            elif 'FRAME' in hdr:
                tm_start = Time.Time(hdr['FRAME'])
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
                hdr.set('GAIN', gain, 1.0)
            if 'RDNOISE' not in hdr:
                hdr.set('RDNOISE', rdnoise, 1.0)
            hdr.set('OBSGEO-B', obslat, 'Latitude of the observatory')
            hdr.set('OBSGEO-L', obslon, 'Longitude of the observatory')
            hdr.set('OBSGEO-H', obsalt, 'Altitude of the observatory')
            hdr.set('OBSERVAT', obsname, 'Observatory')

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

            if (objnames[ii].lower() != "flat") and (objnames[ii].lower() != "thar") and (objnames[ii].lower() != "bias") and (objnames[ii].lower() != "sky") and (objnames[ii].lower() != "dark"):
                if 'RA' in hdr and 'DEC' in hdr:
                    if hdr['RA'] != '' and hdr['DEC'] != '':
                        ra = hdr['RA']; dec = hdr['DEC']
                elif 'MOUNTRA' in hdr and 'MOUNTDE' in hdr:
                    if hdr['MOUNTRA'] != '' and hdr['MOUNTDE'] != '':
                        ra = hdr['MOUNTRA']; dec = hdr['MOUNTDE']
                else:
                    print(f"Requesting information for {objnames[ii]} from SIMBAD")
                    simbad_session = Simbad()
                    try:
                        query_result = simbad_session.query_object(objnames[ii])
                    except Exception as e:
                        print(f"Error: {e}. Leave the wavelength scale uncorrected")
                        ra = ''; dec = ''
                    else:
                        coo = coord.SkyCoord(ra=query_result['ra'], dec=query_result['dec'])
                        ra = coo.ra.to(u.hourangle).to_string(sep=":", precision=2)[0]
                        dec = coo.dec.to(u.degree).to_string(sep=":", precision=2)[0]
                        hdr.set('RA', ra, 'RA in hours')
                        hdr.set('DEC', dec, 'DEC in degrees')
                if 'EPOCH' not in hdr:
                    hdr.set('EPOCH', 2000., 'EPOCH of coordinates')
                observat = coord.EarthLocation.from_geodetic(obslon, obslat, obsalt * u. m)
                dateobs = np.char.replace(tm_mid.fits, 'T', ' ')
                dateobs = Time.Time(dateobs, scale='utc', location=observat)
                # if ra != '' and dec != '':
                #     star = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
                #     ltt_bary = dateobs.light_travel_time(star)
                #     bjd = dateobs.jd + ltt_bary.value
                #     bcr = star.radial_velocity_correction(obstime=dateobs)
                #     hdr.set('BJD', bjd, 'Barycentric JD')
                #     hdr.set('BARYCORR', bcr.to(u.km/u.s).value, 'Barycentric correction')
                # else:
                hdr.set('JD', dateobs.jd, 'Julian Date')
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
            hdu[0].header = hdr
            hdu.flush()
            print("File %s has been updated" %(files[ii].strip()))
    print("...done")
    return None
