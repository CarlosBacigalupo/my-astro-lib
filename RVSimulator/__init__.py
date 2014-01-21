
import numpy as np
from PyAstronomy import pyTiming
import csv
import pylab as plt
from matplotlib.ticker import ScalarFormatter
import pyfits as pf
from math import sqrt, cos, pi
from scipy import signal, interpolate, optimize

def clean_flux(flux, xDef = 1, lambdas = []):
	'''clean a flux array to cross corralate to determine RV shift
		eliminates NaNs
		moving median to reduce peaks
		optional: increase resolution by xDef times
		
	'''	
	
	#Copy to output in case of no resampling
	fluxHD = flux
	newLambdas = lambdas
	
	#if enought data -> resample
	if ((xDef>1) and (len(lambdas)>0)):
		fFluxHD = interpolate.interp1d(lambdas,flux) 
		newLambdas = np.arange(min(lambdas), max(lambdas),(max(lambdas)-min(lambdas))/len(lambdas)/xDef)
		fluxHD = fFluxHD(newLambdas)
	
	#clean NaNs and median outliers	
	fluxHD[np.isnan(fluxHD)] = 0
	fluxNeat = fluxHD	
	fluxMed = signal.medfilt(fluxHD,5)
	w = np.where(abs((fluxHD-fluxMed)/np.maximum(fluxMed,50)) > 0.4)
	for ix in w[0]:
		fluxNeat[ix] = fluxMed[ix]
	
	#remove trailing zeros, devide by fitted curve (flatten) and apply tukey window
	fluxNeat = np.trim_zeros(fluxNeat,'f') 
	newLambdas = newLambdas[-len(fluxNeat):]
	fluxNeat = np.trim_zeros(fluxNeat,'b') 
	newLambdas = newLambdas[:len(fluxNeat)]
	
	fFluxNeat = optimize.curve_fit(cubic, newLambdas, fluxNeat, p0 = [1,1,1,1])
	fittedCurve = cubic(newLambdas, fFluxNeat[0][0], fFluxNeat[0][1], fFluxNeat[0][2], fFluxNeat[0][3])
# 	plt.plot(fittedCurve)
# 	plt.plot(fluxNeat)
# 	plt.show()
# 	plt.plot(fluxNeat/fittedCurve-1)
# 	plt.show()
	
	fluxFlat = fluxNeat/fittedCurve-1
	
	fluxWindow = fluxFlat * tukey(0.1, len(fluxFlat))
	
	return newLambdas, fluxWindow

def tukey(alpha, N):

	tukey = np.zeros(N)
	for i in range(int(alpha*(N-1)/2)):
		tukey[i] = 0.5*(1+cos(pi*(2*i/alpha/(N-1)-1)))
	for i in range(int(alpha*(N-1)/2),int((N-1)*(1-alpha/2))):
		tukey[i] = 1
	for i in range(int((N-1)*(1-alpha/2)),int((N-1))):
		tukey[i] = 0.5*(1+cos(pi*(2*i/alpha/(N-1)-2/alpha+1)))
	
	return tukey
def quad(x,a,b,c):
	
	return a*x**2+b*x+c

def cubic(x,a,b,c,d):
	
	return a*x**3+b*x**2+c*x+d


def resolving_power():
	'''
    Produces different plots of RV precision calculations.
    
    Parameters
    ----------
    None
         
    Returns
    -------
    Nothing
        
    Notes
    -----
    Outputs RV precision as a function of wavelength.
    '''
	
	R = 45000
	
	Lambda1 = np.arange(4718,4903)
	Lambda2 = np.arange(5649,5873)
	Lambda3 = np.arange(6481,6739)
	Lambda4 = np.arange(7590,7890)
	Lambda = np.hstack((Lambda1, Lambda2, Lambda3, Lambda4))
	
	deltaLamda = 1
	# deltaLamda = Lambda /R
	
	rv = c/Lambda * deltaLamda
	
	#########Plotting mess ahead
	# ax1 = fig.add_subplot(111, axisbg = 'black')
	#        im = np.sqrt(im) #Remove this line for Hg
	#    
	#        labels = np.array([0])
	#         
	#        for line in open('solar2.txt'):
	#            Lambda = float(str(line).split()[0]) #Wavelength
	#            labels = np.vstack((labels,np.array([Lambda])))
	#
	#        labels=labels[1:]       
	#        im=mpimg.imread('solar.png')
	#        color=colorTable
	#        print random.randrange(-30,-10) random()
	#        plt.subplots_adjust(bottom = 0.1)
	#        plt.plot( x,-z, "o", markersize=7, color=colorTable, markeredgewidth=1,markeredgecolor='g', markerfacecolor='None' )
	# plt.imshow(im,extent=[-imWidth/2 , imWidth/2 , -imHeight/2 , imHeight/2])
	# plt.set_cmap(cm.Greys_r)
	# ax1.scatter(CCDX, -CCDY ,s=8, color=colorTable , marker='o', alpha =.5)
	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.scatter(Lambda, rv)
	plt.axis()
	y_formatter = ScalarFormatter(useOffset=False)
	ax1.yaxis.set_major_formatter(y_formatter)
	plt.title('RV resolution at R=' + str(R))
	plt.ylabel('delta RV [m/s]')
	plt.xlabel('wavelength (nm)')       
	plt.show()

def W(Lambda, A0):
	'''
    Calculates the weight function form Bouchy 2001.
    
    Parameters
    ----------
    Lambda : np.array
        n x 1 np.array with wavelength bins
        
    AO : np.array
        n x 1 np.array with counts
    
            
    Returns
    -------
    W : np.array
        n x 1 np.array weights as a function of pixel.
        
    Notes
    -----
    Lambda and A0 should be equal length.
    Uses:
    W(i) = Lambda(i)**2 (dA0(i)/dLambda(i))**2 / A0(i)
    Assumes noise free detector. (No sigma_D**2 term in the denominator).
    dA0(i)/dLambda(i) simplified as discrete DeltaY/DeltaX.
    '''

	dA0dL = np.zeros(len(A0)-1)
	
	for i in range(len(A0)-2): #compute partial derivative
		dA0dL[i] = (A0[i+1] - A0[i])/(Lambda[i+1] - Lambda[i])

	#compute W (removing last term from Lambda and A0 as dA0dL has n-1 terms.
	W = Lambda[:-1]**2 * dA0dL**2 / A0[:-1]
	
	return W


def Q(W, A0):
	'''
    Calculates the Q factor of a spectrum from W(weight) and A0(flux) form Bouchy 2001.
    
    Parameters
    ----------
    W : np.array
        n x 1 np.array weight 
        
    AO : np.array
        n x 1 np.array with flux counts
    
            
    Returns
    -------
    Q : float
        Quality factor. 
        
    Notes
    -----

    '''
	A0[np.isnan(A0)] = 0
	Q = sqrt(np.sum(W)/np.sum(A0))
	
	return Q

def QdRV(Lambda, A0):
	
	W1 = W(Lambda, A0)
	W1[np.isnan(W1)] = 0
	Q = Q(W1, A0)
	dRV = c/sqrt(np.sum(W1))
	
	return Q, dRV
def extract_single_from_2dfdr(fileName, Y):
	'''
    Extracts a single fibre from a reduced HERMES fits file.
    
    Parameters
    ----------
    fileName :  str
        Name of the fits file to be extracted
    
    Y : int
        Y-coord of the fibre to be extracted
     
    Returns
    -------
    Lambda : np.array
        n x 1 np.array with wavelength bins
        
    AO : np.array
        n x 1 np.array with counts
     
    Notes
    -----
    The reduced file is simply a compilation of the extracted fibres. 
    Tramlines, extraction methods and other steps are computed in 2dfdr.  
    '''
	a = pf.open(fileName)

	A0 = a[0].data[Y,:]             	
	Lambda = extract_HERMES_wavelength(fileName)
	
	return Lambda, A0

def extract_all_from_2dfdr(fileName):
	'''
    Extracts all fibres from a reduced HERMES fits file.
    
    Parameters
    ----------
    fileName :  str
        Name of the fits file to be extracted
     
    Returns
    -------
    Lambda : np.array
        n x 1 np.array with wavelength bins
        
    AO : np.array
        n x 1 np.array with counts
     
    Notes
    -----
    The reduced file is simply a compilation of the extracted fibres. 
    Tramlines, extraction methods and other steps are computed in 2dfdr.  
    '''
	a = pf.open(fileName)

	A = a[0].data             	
	
	return A


def extract_HERMES_wavelength(fileName):

	a = pf.open(fileName)

	CRVAL1 = a[0].header['CRVAL1'] # / Co-ordinate value of axis 1                    
	CDELT1 = a[0].header['CDELT1'] #  / Co-ordinate increment along axis 1             
	CRPIX1 = a[0].header['CRPIX1'] #  / Reference pixel along axis 1                   
	
	#Creates an array of offset wavelength from the referece px/wavelength
	Lambda = CRVAL1 - (CRPIX1 - (np.arange(int(CRPIX1)*2)) -1)* CDELT1

	return Lambda



	

