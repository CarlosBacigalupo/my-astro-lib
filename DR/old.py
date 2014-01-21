#!/usr/bin/env python

#Python routine to calibrated RHEA data for bias, darks, bad pixel masks. This calibration includes the flat field files and combines them into a master nightly flat if at all available. 

import numpy as np
import pylab as pl
import scipy, pyfits, commands,os
from pyraf import iraf


def ask(text,def_val):
    temp=raw_input(text+' = '+str(def_val)+' =')
    if len(temp)==0:
        return def_val
    else:
        return type(def_val)(temp)

iraf.images(_doprint=0)     # load noao
if os.path.exists('flat.fits'): os.system('rm flat.fits')
#Assumes all images in the current directory are from the same camera, and that there are biases and darks here.
lis=commands.getoutput('ls *.fits')
lis=lis.split('\n')
path=commands.getoutput('pwd')+'/'

#get the calibration files
try: 
    bias=pyfits.getdata('../calibration_frames/mbias.fits')
    dark=pyfits.getdata('../calibration_frames/mdark.fits')
    dark_exposure=pyfits.getheader('../calibration_frames/mdark.fits')['EXPTIME']
    bad=pyfits.getdata('../calibration_frames/badpix.fits')
except Exception:
    print 'Unable to import calibration files. Assumes they are in a calibration_frames directory in the root directory of the data set.'
    sys.exit()

print 'Successfully imported calibration files'

if os.path.exists('calibrated/'): os.system('rm -R calibrated/')
os.system('mkdir calibrated')

fl=ask('Does this night have any flats to be combined?','y')

if fl=='y':
    if os.path.exists('calibrated_flats/'): os.system('rm -R calibrated_flats/')
    os.system('mkdir calibrated_flats')
    print 'Started calibrating flats'
    flats=[]
    #create master flat for the night from a combination of all flat frames for the night.
    for i in lis:
        try: f=pyfits.open(i)
        except Exception: print 'Unable to open image '+i
        h=f[0].header
        exp=h['EXPTIME']
        imtype=h['IMGTYPE']
        if imtype=='Flat':
            d=f[0].data
            #calibration
            d-=bias
            d-=dark*(exp/dark_exposure)
            d*=bad
            h.update('BIASED','True','Set to true if image is de-biased')
            h.update('DARKED','True','Set to true if image is dark subtracted')
            h.update('BADPIX','True','True if image has bad pixels as NaNs')
            pyfits.writeto('calibrated_flats/'+i,d,header=h)
            flats.append('calibrated_flats/'+i)
    print 'Finished calibrating flats'
    print 'combining flats into nightly master flat'
    np.savetxt('flatlist',flats,fmt='%30s',delimiter='\n')
    iraf.imcombine(input='@flatlist',output='flat.fits',combine='median')
    f=pyfits.open('flat.fits',mode='update')
    flat_temp=f[0].data
    flat_temp[flat_temp==0]=np.nan
    f.flush()
    print max(flat_temp.flatten())
    iraf.imarith(operand1='flat.fits',op='/',operand2=max(flat_temp.flatten()),result='flat.fits')
    f=pyfits.open('flat.fits',mode='update')
    flat_temp=f[0].data
    flat_temp[flat_temp==0]=np.nan
    f.flush()
    flat=pyfits.getdata('flat.fits')
    os.system('rm flatlist')

print 'Started calibrating images'
for i in lis:
    try: f=pyfits.open(i)
    except Exception: print 'Unable to open image '+i
    h=f[0].header
    exp=h['EXPTIME']
    imtype=h['IMGTYPE']
    if imtype=='Light' or imtype=='HgAr':
        d=f[0].data
        #calibration calibration
        d-=bias
        d-=dark*(exp/dark_exposure)
        d*=bad
        h.update('BIASED','True','Set to true if image is de-biased')
        h.update('DARKED','True','Set to true if image is dark subtracted')
        h.update('BADPIX','True','True if image has bad pixels as NaNs')
#        if imtype=='Light':
#            d/=flat
#            h.update('FLATED','True','True if image has been flat fielded')
        pyfits.writeto('calibrated/'+i,d,header=h)
    
print 'Finished calibrating images'
