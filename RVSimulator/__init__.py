
import sys
import isomod
import amodesmod as amod
import noisemod as noise
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

def read_fibre_table(fileName):
#Reads a 2df fits file and returns the data from the fibre HDU (fibre table)

	a = pf.open(fileName)
	for i in range(1,len(a)):
		fibreType = a[i].header['EXTNAME']
		if fibreType=='FIBRES' :b = np.array(a[2].data)
	return b


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


def load_2dfdr_fibre_table(fileName):
	
	a = pf.open(fileName)	
	fibreTable = a[2].data
	
	return fibreTable

	
def movingaverage(data, window_size):
	# Smoothing function. It takes in the timeseries that is to be smoothed and the size of the window of points to be averaged over.

	window = numpy.ones(int(window_size))/float(window_size)
	return numpy.convolve(data, window, 'same')

def main():
	"""
	This program reads in an isochrone produced from http://stev.oapd.inaf.it/cgi-bin/cmd, along with the mass of the star to be considered.
		From this, it draws the other necessary parameters for the star (luminousity and temperature) from the isochrone. It then calculates the expected mode frequencies and amplitudes using equations
        from Huber et.al (2011), Kjeldsen & Bedding (1995) and Lecture Notes on Stellar Oscillations by J. Christensen-Dalsgaard [5th Ed. 2003/03].
	The program then sums the sines of these frequencies along with noise of specified amplitude to produced the ideal signal. From this we can then use the program to test various scenarios with different sampling rates and observation length and see the results.
	The programs uses the Lomb-Scargle method to analyse the sampled signal.Many modifications can be made in order to organise your output and how the program runs.
		The purpose of this program is to be able to run simlutions on measuring mode frequenceis of different stars to find suitable targets for the HERMES spectrographs to search for exoplanets around.
	"""

	# Uses the isochrone module to read model. It then returns and sets the values of Luminousity and Temperature.
	# The isochrone may not have the precise mass you are looking for, and will instead return the values corresponding to the closest mass found.
	isomod.set_data(infile, mass)     #sets data of isochrone to be read and the mass to be found.
	lumin_t,temp_t=isomod.run()       #runs isomod.py to find the corresponding luminousity and temperature for a giant star.
	lumin = 10**float(lumin_t)        #Converts the log given in the model to be given in units of solar luminousty.
	temp = 10**float(temp_t)          #Converts the log given in the model to give the effective temperature.
	
	# Uses the amodes2 module to calculates the frequencies and amplitudes of different modes.
	amod.setconst(mass, lumin, temp)  #set data to calculate all frequencies and amplitudes for specified n and l ranges.
	data = amod.modefreq()            #(n_min=-1, n_max=1, l_max=3, l_min=0). This function returns 2 list: The first list contains the calculated frequencies (data[0]) and the second list contains their corresponding amplitudes (data[1]).
	nu_max = amod.return_nu_max()
	
	timelist = (days*24*3600)/timeSampling #calculates the size of time array neccesary for the specified sampling rate and observation times.
	
	print "Mass(Mo):", mass, ", Luminousity(Lo):",lumin, ", Temperature(K):", temp, "nu_max:", nu_max
	print "Peak frequency(microHz)(ls)	Peak frequency(microsHz)(smooth)	Peak amplitude(ls)	Peak amplitude(smooth)	Peak STD	Good Star?"
	
	with open(outfile + ".csv", 'wb') as csvfile:
		output = csv.writer(csvfile, delimiter=',')
		header1 = ["Mass(Mo):", mass, "Luminousity(Lo):", lumin, "Temperature(K):", temp, "nu_max:", nu_max]
		header2 = ["Peak frequency(microHz)(ls)", "Peak frequency(microsHz)(smooth)", "Peak amplitude(ls)", "Peak amplitude(smooth)", "Good Star?"]
		output.writerow(header1)
		output.writerow(header2)
		for i in range(iterations):
			time = timeSampling*numpy.arange(timelist)                              #defines time domain
			noise.setdata(1e-6*data[0].flatten(), data[1].flatten(), time, obsErr)  #noise.setdata(frequencies, amplitudes, time, observational error):
			timeseries = noise.make_noise()
			#timeseries = noise.add_planet(series, 1, 100)
			
			# Algorithm that takes out chunks of data from timeseries
			sample = 0
			block = (obsPeriod*3600)/timeSampling # number of measurements taken while observing
			spaceInt = 0
			newdata = []
	
			while (sample < len(timeseries)):
				#block = numpy.random.random_integers(20, 60)
				if spaceInt == len(spaceLength):
					spaceInt=0
				for i in numpy.arange(block):
					if sample + i < len(timeseries):
						newdata.append(timeseries[sample + i])
				sample = sample + block
				#space = numpy. random. random_integers(60, 180)
				space = ((spaceLength[spaceInt]) *3600)/timeSampling #number of measurements that occur outside of the observing period
				spaceInt = spaceInt + 1
				for j in numpy.arange(space):
					if sample + j < len(timeseries):
						newdata.append(0)
				sample = sample + space
			
			#This section plots the new time series for the non-continous observing, along with a descrete fourier transform that the new time series.
			plt.figure()
			fig1 = plt.subplot(111)
			plt.plot(time,newdata)
			fig1.set_title("Signal")
			fig1.set_xlabel("Time (s)")
			fig1.set_ylabel("Amplitude (m/s)")
			'''
			fig2 = plt.subplot(212)
			ftvel = numpy.abs(numpy.fft.fft(newdata))
			plt.plot(ftvel)
			fig2.set_title("FFT of signal")
			fig2.set_xlabel("Frequency ($\mu Hz$)")
			fig2.set_ylabel("Amplitude")
			plt.draw()
			'''
			#The Lomb Scargle module we found:
			t = time
			x = newdata
			# PyAstronomy stuffs
			ts=pyTiming.pyPeriod.TimeSeries(t,x,error=numpy.ones(numpy.size(x))*0.00000001)
			ft=pyTiming.pyPeriod.Fourier(ts)
			#print ft.FAP
			#Lets hard-wire a frequency separation of 0.5 cycles per data series length. For evenly sampled data,
			#this would mean that every second point is independent
			fsep = 0.5/(numpy.max(t)-numpy.min(t))
			ls=pyTiming.pyPeriod.Gls(ts,freq=numpy.arange(fsep,ts.returnNyquist(),fsep ))
	
			#ls.plot()
			goodStar = False                       #A variable that it used it determine whether the star we are looking at is suitable
			smoothls = movingaverage(ls.power, 30) #smooths lomb-scargle
			medls = numpy.median(smoothls)         #Calculates the median of the smoothed lomb-scargle
			medlsVal = [medls]
			if (numpy.max(smoothls) > 2*medls):    #determines whether the star used is a good star, by checking whether the peak power is abover a certain threshold. In this case the threshold is 2* the median amplitude of the Lomb-Scargle periodogram
				goodStar = True
			'''
			stdList = [] #list of standard deviations for smoothls
			for i in range(len(ls.freq)):	#This for loop calculated the standard deviation over range of 50 points in smoothls
				binSet = []
				if i < 25:
					for j in range(50):
						binSet.append(smoothls[j])
				elif i > len(ls.freq)-25:
					this = len(ls.freq)-50
					for j in range(50):
						binSet.append(smoothls[this + j])
				else:
					for j in range(i-25, i+25):
						binSet.append(smoothls[j])
				std = numpy.std(binSet)
				stdList.append(std)
	
			plt.figure()
			fig1 = plt.subplot(111)
			plt.plot(ls.freq, stdList)
			fig1.set_title("Standard deviation of each frequency")
			fig1.set_xlabel("Time (s)")
			fig1.set_ylabel("STD")
			'''
			plt.figure()
			fig2 = plt.subplot(111)
			plt.plot(ls.freq, smoothls)
			plt.plot(ls.freq, numpy.ones(len(smoothls))*medls, ':')
			plt.plot(ls.freq, numpy.ones(len(smoothls))*medls*2, ':')	
			fig2.set_title("Smoothed Lomb-Scargle Periodogram")
			fig2.set_xlabel("Frequency (Hz)")
			fig2.set_ylabel("Power (probably (m/s)^2)")
			plt.show(block=False)
	
	
			freqInd = 0
			for j in numpy.arange(len(ls.power)):
				if ls.power[j] == max(ls.power):
					freqInd = j
			smoothInd = 0
			for k in numpy.arange(len(smoothls)):
				if smoothls[k] == max(smoothls):
					smoothInd = k
			print (ls.freq[freqInd])*1000000, (ls.freq[smoothInd])*1000000, max(ls.power), max(smoothls), goodStar
			printresults = [(ls.freq[freqInd])*1000000, (ls.freq[smoothInd])*1000000, max(ls.power), max(smoothls), goodStar]
			output.writerow(printresults)
		plt.show()
		output.writerow(medlsVal)
		print medls
		
# Sets all variables
infile = "test_iso2.csv"			#sys.argv[1], File containing isochrone
mass = 1 #float(sys.argv[1])
lumin = 0
temp = 0
outfile = 'out.csv' #sys.argv[2]
timeSampling = 900 				#Sampling rate given in seconds eg. sampling every 10 minutes = 600 seconds.
obsErr = 10					#The amount of noise you wish to add to your timeseries, in meters per second
days = 4 					#How long in days you will be observing straight
obsPeriod = 0.5 				#Hours of observing at night
spaceLength = [1.0,1.0,1.0,1.0,1.0,16]		#List of lengths of spaces between observation periods in hours.
medlsVal = []  					#Median amplitude value of the smoothed lomb scargle. This is determined in the program
iterations = 1					#Number of times you would like ot run the script

