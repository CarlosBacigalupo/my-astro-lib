from pyraf import iraf
import os
import numpy as np
from constants import *
import pylab as pl


def find_nearest(array,value):
    '''Function used to find the neareast element of an array to a given value. Usage: find_nearest(array,value) '''
    idx = (np.abs(array-value)).argmin()
    return idx 
    
    
def stitch(specLambda, specFlux):
    '''Function used to stitch spectral orders together. usage: stitch(specLambda, specFlux) where specLambda is a list of arrays of wavelengths and specFlux is the corresponding list of arrays of flux values.
    THIS FUNCTION REQUIRES SUBSTANTIAL UPDATE AND MODIFICATIONS. USE WITH CAUTION.'''
    fullLambda = []
    fullFlux = []
    
    #start from the highest order (lowest wavelengths)
    for i in range(len(specLambda)-1,0,-1):
        fullLambda = np.append(fullLambda, specLambda[i])
        fullFlux = np.append(fullFlux, specFlux[i])
        
#         fullIdx = find_nearest(fullLambda,min(specLambda[i-1]))
#         newIdx = find_nearest(specLambda[i-1],max(fullLambda))
        
        #find the mid point of the overlap region and the corresponding closest array elements to it
        midPoint = (max(fullLambda) - min(specLambda[i-1]))/2 + min(specLambda[i-1])
        fullIdx = find_nearest(fullLambda,midPoint)
        newIdx = find_nearest(specLambda[i-1],midPoint)
        
        #Now cut the lowest order up to the mid point and the highest wavelength order from the mid point onwards
        fullLambda = fullLambda[:fullIdx]
        fullFlux = fullFlux[:fullIdx]
        specLambda = specLambda[:i-1] + (specLambda[i-1][newIdx:],) + specLambda[i:]
        specFlux = specFlux[:i-1] + (specFlux[i-1][newIdx:],) + specFlux[i:]
       

# 
#         print len(fullLambda[fullIdx:])
#         print len(specLambda[i-1][:newIdx])
        
    #This is necessary since there is no need to find the mid point on the last order.
    midPoint = (max(fullLambda) - min(specLambda[0]))/2 + min(specLambda[0])
    fullIdx = find_nearest(fullLambda,midPoint)
    newIdx = find_nearest(specLambda[0],midPoint)
    
    fullLambda = fullLambda[:fullIdx]
    fullFlux = fullFlux[:fullIdx]
    specLambda = (specLambda[0][newIdx:],) + specLambda[1:]
    specFlux = (specFlux[0][newIdx:],) + specFlux[1:]
    
    fullLambda = np.append(fullLambda, specLambda[0])
    fullFlux = np.append(fullFlux, specFlux[0])
    
    return fullLambda, fullFlux
    
def divide_flat_full(specOrderStar, specLambdaStar, specFluxStar, specOrderFlat, specLambdaFlat, specFluxFlat, booPlot = False):
    
    if min(specLambdaStar[0])<min(specLambdaFlat[0]):
        #remove from low end star x 2
        map = specLambdaStar[0]>=min(specLambdaFlat[0])
        specLambdaStar = (specLambdaStar[0][map],) + specLambdaStar[1:]
        specFluxStar = (specFluxStar[0][map],) + specFluxStar[1:]
        
    elif min(specLambdaStar[0])>min(specLambdaFlat[0]):
        #remove from low end flat x 2
        map = specLambdaFlat[0]>=min(specLambdaStar[0])
        specLambdaFlat = (specLambdaFlat[0][map],) + specLambdaFlat[1:]
        specFluxFlat = (specFluxFlat[0][map],) + specFluxFlat[1:]
    
    if max(specLambdaStar[0])>max(specLambdaFlat[0]):
        #remove from high end star x 2
        map = specLambdaStar[0]<=max(specLambdaFlat[0])
        specLambdaStar = (specLambdaStar[0][map],) + specLambdaStar[1:]
        specFluxStar = (specFluxStar[0][map],) + specFluxStar[1:]

    elif max(specLambdaStar[0])<max(specLambdaFlat[0]):
        #remove from high end flat x 2
        map = specLambdaFlat[0]<=max(specLambdaStar[0])
        specLambdaFlat = (specLambdaFlat[0][map],) + specLambdaFlat[1:]
        specFluxFlat = (specFluxFlat[0][map],) + specFluxFlat[1:]

    if len(specLambdaStar[0])>len(specLambdaFlat[0]):
        #reduce length star
        specLambdaStar = (specLambdaStar[0][:len(specLambdaFlat[0])],) + specLambdaStar[1:]
        specFluxStar = (specFluxStar[0][:len(specFluxFlat[0])],) + specFluxStar[1:]
        
    elif len(specLambdaStar[0])<len(specLambdaFlat[0]):
        #reduce lenght flat
        specLambdaFlat = (specLambdaFlat[0][:len(specLambdaStar[0])],) + specLambdaFlat[1:]
        specFluxFlat = (specFluxFlat[0][:len(specFluxStar[0])],) + specFluxFlat[1:]
    
    
    flatLambda = (specLambdaFlat[0],)
    flatFlux = (specFluxStar[0]/(specFluxFlat[0]/max(specFluxFlat[0])),)
    
    if booPlot:
        pl.plot(specLambdaFlat[0],specFluxFlat[0]/max(specFluxFlat[0]))
        pl.plot(specLambdaStar[0],specFluxStar[0]/max(specFluxStar[0]))
        pl.show()
        
    for i in range(1, len(specOrderStar)):
        
        if min(specLambdaStar[i])<min(specLambdaFlat[i]):
            #remove from low end star x 2
            map = specLambdaStar[i]>=min(specLambdaFlat[i])
            specLambdaStar = specLambdaStar[:i] + (specLambdaStar[i][map],) + specLambdaStar[i+1:]
            specFluxStar = specFluxStar[:i] + (specFluxStar[i][map],) + specFluxStar[i+1:]
            
        elif min(specLambdaStar[i])>min(specLambdaFlat[i]):
            #remove from low end flat x 2
            map = specLambdaFlat[i]>=min(specLambdaStar[i])
            specLambdaFlat = specLambdaFlat[:i] + (specLambdaFlat[i][map],) + specLambdaFlat[i+1:]
            specFluxFlat = specFluxFlat[:i] + (specFluxFlat[i][map],) + specFluxFlat[i+1:]
        
        if max(specLambdaStar[i])>max(specLambdaFlat[i]):
            #remove from high end star x 2
            map = specLambdaStar[i]<=max(specLambdaFlat[i])
            specLambdaStar = specLambdaStar[:i] + (specLambdaStar[i][map],) + specLambdaStar[i+1:]
            specFluxStar = specFluxStar[:i] + (specFluxStar[i][map],) + specFluxStar[i+1:]
    
        elif max(specLambdaStar[i])<max(specLambdaFlat[i]):
            #remove from high end flat x 2
            map = specLambdaFlat[i]<=max(specLambdaStar[i])
            specLambdaFlat = specLambdaFlat[:i] + (specLambdaFlat[i][map],) + specLambdaFlat[i+1:]
            specFluxFlat = specFluxFlat[:i] + (specFluxFlat[i][map],) + specFluxFlat[i+1:]
    
        if len(specLambdaStar[i])>len(specLambdaFlat[i]):
            #reduce length star
            specLambdaStar = specLambdaStar[:i] + (specLambdaStar[i][:len(specLambdaFlat[i])],) + specLambdaStar[i+1:]
            specFluxStar = specFluxStar[:i] + (specFluxStar[i][:len(specFluxFlat[i])],) + specFluxStar[i+1:]           
        elif len(specLambdaStar[i])<len(specLambdaFlat[i]):
            #reduce lenght flat
            specLambdaFlat = specLambdaFlat[:i] + (specLambdaFlat[i][:len(specLambdaStar[i])],) + specLambdaFlat[i+1:]
            specFluxFlat = specFluxFlat[:i] + (specFluxFlat[i][:len(specFluxStar[i])],) + specFluxFlat[i+1:]
        
        flatFlux = flatFlux + (specFluxStar[i]/(specFluxFlat[i]/max(specFluxFlat[i])),)
        flatLambda = flatLambda + (specLambdaFlat[i],)
        
    return flatLambda, flatFlux
    
def check_if_file_exists(filename):
    '''Function used to delete a particular file if it exists'''
    if os.path.isfile(filename): os.remove(filename)
    return filename

def medianCombine(inFiles, outFilename):
    '''Function that performs a median combination of frames and outputs a result. usage: medianCombine(inFiles, outFilename) where inFiles is a list of filenames and outFilename is the desired name of the result.'''
    #Start by loading all the frames into a data cube
    try: 
        for k in range(0,len(inFiles)):
            im = pyfits.getdata(inFiles[k])
            if (k == 0):
                cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
            cube[:,:,k] = im
    except Exception: print 'Failed to load files into memory. Are you perhaps trying to combine too many files?'
    #median all the values along the image direction
    med_im=np.median(cube,axis=2)
    #Write resulting file to disk
    try: pyfits.writeto(outFilename,med_im)
    except Exception: print 'Failed to write median combine result to disk'
    
def averageCombine(inFiles, outFilename):
    '''Function that performs an average combination of frames and outputs a result. usage: medianCombine(inFiles, outFilename) where inFiles is a list of filenames and outFilename is the desired name of the result.'''
    #Start by loading all the frames into a data cube
    try: 
        for k in range(0,len(inFiles)):
            im = pyfits.getdata(inFiles[k])
            if (k == 0):
                cube=np.empty([im.shape[0],im.shape[1],len(inFiles)])
            cube[:,:,k] = im
    except Exception: print 'Failed to load files into memory. Are you perhaps trying to combine too many files?'
    #Average all the values in the image direction
    med_im=np.average(cube,axis=2)
    try: pyfits.writeto(outFilename,med_im)
    except Exception: print 'Failed to write average combine result to disk'
    
def batchMedian():
    '''THIS FUNCTION HAS TOO MANY THINGS HARD CODED. Probably written for some sort of one time job...'''
    start = range(1,721,5)
    stop = range(5,721,5)
    os.chdir('/media/win-desktop/6hrs/Darks')
    
#    mink=10
#    maxk=15
#    files = 720/5
    for j in range(720/5-1):
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFiles, outFilename)

def batchDark():
    '''THIS FUNCTION HAS TOO MANY THINGS HARD CODED. Probably written for some sort of one time job...'''
    files = range(1,721)

    os.chdir('/media/win-desktop/6hrs/Darks')

#    mink=10
#    maxk=15
#    files = 720/5
    for j in len(files):
        
        inFile = pyfits.getdata(files[j],0)
        inFiles = ['Hg_20s_%03d.FIT' % k for k in range(start[j],start[j+1])]
        outFilename = 'med' + str(start[j]) + '_' + str(start[j+1]) + '.fits' 
        medianCombine(inFile, darkFile, outFilename)

def calibrateImage(inFileName, biasFileName, darkFileName, flatFileName, outFileName, pixelmaskFileName=''):
    '''Function used to fully calibrate an image. usage: calibrateImage(inFileName, biasFileName, darkFileName, flatFileName, outFileName, pixelmaskFileName='') where all inputs are self explanatory'''
    #Open the science image into memory
    hdulist = pyfits.open(inFileName)
    inExpTime=hdulist[0].header['EXPOSURE']
    
    inFile = hdulist[0].data
   
    #Bias subtract
    biasFile = pyfits.getdata(biasFileName) 
    outFile=inFile-biasFile
    
    #Darks
    hdulist = pyfits.open(darkFileName)
    darkExpTime = hdulist[0].header['EXPOSURE']
    
    darkFile=hdulist[0].data
    darkFile=darkFile/darkExpTime*inExpTime
    
    outFile-=darkFile
    
    #Flats
    flatFile=pyfits.getdata(flatFileName)
    flatFile/=np.max(flatFile.flatten())
    #outFile=outFile/flatFile
    
    #Bad Pixels
    if pixelmaskFileName!='':
        badpix=pyfits.getdata(pixelmaskFileName)
        outFile*=badpix
    
    #write result to file
    pyfits.writeto(outFileName, outFile)

def subtractDark(inFileName, darkFileName, outFilename):
    '''Function that is a wrapper for a dark subtraction. Can be really used to subtract any frame from any other frame. Usage: subtractDark(inFileName, darkFileName, outFilename)'''
    inFile = pyfits.getdata(inFileName) 
    darkFile = pyfits.getdata(darkFileName) 
    
    outFile = inFile - darkFile
    
    pyfits.writeto(outFilename, outFile)
    
def normalise_image(inFile):
    '''Function that normalises an image using the maximum value as a normalising constant. Single input is the file name'''
    outFile = inFile/np.max(inFile)
    return outFile
    
                              
