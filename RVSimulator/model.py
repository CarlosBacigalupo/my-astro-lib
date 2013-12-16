##!/usr/bin/env python

# import sys
import isomod
import amodesmod as amod
import noisemod as noise
import numpy
from matplotlib import pyplot
from PyAstronomy import pyTiming
import csv


'''
	This program reads in an isochrone produced from http://stev.oapd.inaf.it/cgi-bin/cmd, along with the mass of the star to be considered.
		From this, it draws the other necessary parameters for the star (luminousity and temperature) from the isochrone and calculates the expected mode frequencies and amplitudes using equations
        from Huber et.al (2011), Kjeldsen & Bedding (1995) and Lecture Notes on Stellar Oscillations by J. Christensen-Dalsgaard [5th Ed. 2003/03]. The number of frequencies calculated can be specified by the user.
	The program then sums the sines of these frequencies along with some noise to produced the ideal signal. From this we can then use the program to test various scenarios with different sampling rates and observation length and see the results.
	The programs used the Lomb-Scargle method to analyse the sampled signal.Many modifications can be made in order to organise your output and how the program runs.
		The purpose of this program is to be able to run simlutions on measuring mode frequenceis of different stars to find suitable targets for the HERMES spectrographs to search for exoplanets around.
'''

def movingaverage(data, window_size):
	"""
	Defines the smoothing function used later on. It take in the time series and the size od the range of points to be averaged over.
	"""
	window= numpy.ones(int(window_size))/float(window_size)
	return numpy.convolve(data, window, 'same')


# Sets all variables
infile = "test_iso2.csv"	#sys.argv[1], File containing isochrone
#mass = float(sys.argv[1])
mass = 1.395
lumin = 0
temp = 0
#outfile = sys.argv[2]
outfile = 'out'
timeSampling = 2400 				#Sampling rate given in seconds eg. sampling every 10 minutes = 600 seconds.
days = 4 						#How long in days you will be observing straight
obsPeriod = 7 					#Hours of observing at night
spaceLength = [1,1,1,17] 		#List of lengths of spaces between observation periods in hours.
medlsVal = [] 					#Median amplitude value of the smoothed lomb scargle. This is determined in the program
iterations = 100				#Number of times you would like to run the script

																			
# Using the isochrone module to read model. It then returns and sets the values of Luminousity and Temperature.
# The isochrone may not have the precise mass you are looking for, and will instead return the values corresponding to the closest mass found.
isomod.set_data(infile, mass)    #sets data of isochrone to be read and the mass to be found.
lumin_t,temp_t=isomod.run()      #runs isomod.py to find the corresponding luminousity and temperature for a giant star.
lumin = 10**float(lumin_t)       #Converts the log given in the model to be given in units of solar luminousty.
temp = 10**float(temp_t)         #Converts the log given in the model to give the effective temperature.

# Uses the amodes2 module to calculates the frequencies and amplitudes of different modes.
amod.setconst(mass, lumin, temp)  #set data to calculate all frequencies and amplitudes for specified n and l ranges.
data = amod.modefreq()            #was (5,9,3) amod.modefreq(n_min, n_max, l_max, l_min=0). This function returns 2 list: The first list contains the calculated frequencies (data[0]) and the second list contains their corresponding amplitudes (data[1]).
nu_max = amod.return_nu_max()


timelist = (days*24*3600)/timeSampling

header1 = ['Mass[Solar Masses]:',str(mass) ,'Luminosity[Solar Luminosities]:',str(lumin),'Temperature[K]:',str(temp),'nu_max[microHz]:',str(nu_max)]
header2 = ['Peak frequency[microHz](ls)','Peak frequency[microsHz](smooth)','Peak amplitude[m/s](ls)','Peak amplitude[m/s](smooth)','Peak STD','Good Star?']
print ', '.join(header1)
print ', '.join(header2) 

with open(outfile + ".csv", 'wb') as csvfile:
	output = csv.writer(csvfile, delimiter=',')
	output.writerow(header1)
	output.writerow(header2)
	
	
	
	for i in range(iterations):
		time = timeSampling*numpy.arange(timelist)                          #defines time domain
		noise.setdata(1e-6*data[0].flatten(), data[1].flatten(), time, 10)  #noise.setdata(frequencies, amplitudes, time, observational error):
		timeseries = noise.make_noise()
		#timeseries = noise.add_planet(series, 1, 100)
		
		# Algorithm that takes out chunks of data from timeseries, not points
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
		pyplot.figure()
		fig1 = pyplot.subplot(111)
		pyplot.plot(time,newdata)
		fig1.set_title("Signal")
		fig1.set_xlabel("Time (s)")
		fig1.set_ylabel("Amplitude (m/s)")
		'''
		fig2 = pyplot.subplot(212)
		ftvel = numpy.abs(numpy.fft.fft(newdata))
		pyplot.plot(ftvel)
		fig2.set_title("FFT of signal")
		fig2.set_xlabel("Frequency ($\mu Hz$)")
		fig2.set_ylabel("Amplitude")
		pyplot.draw()
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
		smoothls = movingaverage(ls.power, 50) #smooths lomb-scargle
		medls = numpy.median(smoothls)         #Calculates the median of the smoothed lomb-scargle
		medlsVal = [medls]
		if (numpy.max(smoothls) > 2*medls):    #determines whether the star used is a good star, by checking whether the peak power is abover a certain threshold
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
x
		pyplot.figure()
		fig1 = pyplot.subplot(111)
		pyplot.plot(ls.freq, stdList)
		fig1.set_title("Standard deviation of each frequency")
		fig1.set_xlabel("Time (s)")
		fig1.set_ylabel("STD")
		
		'''
		doPlot=False
		if doPlot:
			pyplot.figure()
			fig2 = pyplot.subplot(111)
			pyplot.plot(ls.freq, smoothls)
			pyplot.plot(ls.freq, numpy.ones(len(smoothls))*medls, ':')
			pyplot.plot(ls.freq, numpy.ones(len(smoothls))*medls*2, ':')	
			fig2.set_title("Smoothed Lomb-Scargle Periodogram")
			fig2.set_xlabel("Frequency (Hz)")
			fig2.set_ylabel("Power (probably (m/s)^2)")
			pyplot.show(block=False)
		
		freqInd = 0
		for j in numpy.arange(len(ls.power)):
			if ls.power[j] == max(ls.power):
				freqInd = j
		smoothInd = 0
		for k in numpy.arange(len(smoothls)):
			if smoothls[k] == max(smoothls):
				smoothInd = k
		print (ls.freq[freqInd])*1000000, (ls.freq[smoothInd])*1000000, max(ls.power), max(smoothls), goodStar
		printresults = [(ls.freq[freqInd])*1000000, (ls.freq[smoothInd])*1000000, max(ls.power), max(smoothls), 0 ,goodStar]
		output.writerow(printresults)
	if doPlot: pyplot.show()
	output.writerow(medlsVal)
	print medls



