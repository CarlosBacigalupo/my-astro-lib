import sys
import numpy as np
import math
from PyAstronomy import pyTiming
import pylab as plt
# import isomod
# import amodesmod as amod
import noisemod as noise



class star():
    
    mass = 0  #[M_sun]
    delta_nu = 0 #[microHz]
    A_max = 0 #[m/s]
    not_a_giant_flag = True 
    isoMass = 0  #[M_sun]
    isoLogTemp = 0  #temperature found by the isochrone 
    isoLogLumin = 0
    isoTemp = 0
    isoLumin = 0
    luminLimit = 1
    isochrone_mass = []
    isochrone_logLumin = []
    isochrone_logTemp = []
    isochrone_fileName = 'iso_Z0_02_t4e9.dat'
    solar_delta_nu = 0 #[microHz]
    solar_temp = 0 #[Kelvin]
    solar_A_max = 0  #[m/s]
    solar_nu_max = 0 #[microHz]
    output = True
    
    def __init__(self):
        self.load_solar_data()
    
    def load_solar_data(self, solar_delta_nu = 135.1, solar_temp = 5777., solar_A_max = 0.2, solar_nu_max = 3090.):
        self.solar_delta_nu = solar_delta_nu
        self.solar_temp = solar_temp
        self.solar_A_max = solar_A_max
        self.solar_nu_max = solar_nu_max

            
    def read_isochrone(self, isochrone_fileName = ''):
        '''
        Reads in the isochrone model in order ot be able to later determine the luminousity and temperature of the input mass in the model.py program. 
        The input for this function is the isochrone in .csv format.
        '''
        if isochrone_fileName != '':
            self.isochrone_fileName = isochrone_fileName
        
        if (self.isochrone_fileName == ''):
            print 'No isochrome file name specified.'
            print 'Usage: <class_name>.isochrone_fileName = filename '
            return
        
        try:
            a = np.loadtxt(self.isochrone_fileName ,skiprows=8)
        except:
            print 'Could not open file ' + str(self.isochrone_fileName)
            print sys.exc_info()[0]
            return
        
        self.isochrone_mass = a[:,1]
        self.isochrone_logTemp=a[:,2]
        self.isochrone_logLumin=a[:,4]
        self.isochrone_BMag = a[:,6]
        self.isochrone_VMag = a[:,7]
        
        if self.output==True:print 'Isochrone ' + self.isochrone_fileName + ' loaded.'
        
        
    def populate_stellar_parameters(self):
        '''
        This function finds the nearest mass match in the isochrone to the input mass. 
        It also checks whether the luminousty is above the declared luminousity threshold.
        This threshold is to make sure the returned values correspond to a giant star and not a main sequence star.
        The input for this function is a mass to be found and the isochrones list of mass, luminousity and temperature.
        '''
        if len(self.isochrone_mass)==0:
            print 'Isochrone not loaded. Use <class_name>.read_isochrone().'
            return
        
        if self.mass==0:
            print 'Reference stellar mass not set. Use <class_name>.mass = <stellar_mass>.'
            return
        
        distance_array = np.abs( self.isochrone_mass - self.mass )
        closest_point = np.min(distance_array)
        closest_point_index = np.where( distance_array == closest_point )[0][0]    
           
        self.isoMass = self.isochrone_mass[closest_point_index]
        self.isoLogTemp = self.isochrone_logTemp[closest_point_index]   
        self.isoLogLumin = self.isochrone_logLumin[closest_point_index]
        self.isoTemp = 10**self.isoLogTemp
        self.isoLumin = 10**self.isoLogLumin
        if self.output==True:print 'Mass, Luminosity and Temperature read from isochrone.'
        
        self.calculate_A_max()
        self.calculate_delta_nu()
        self.calculate_nu_max()
        self.calculate_freq_modes()
        self.populate_giant_flag()
        
        
        
    def populate_giant_flag(self):
        
        if self.isoLumin < self.luminLimit:
            self.not_a_giant_flag = True
            if self.output==True:print 'Not a giant flag set'


    def calculate_delta_nu(self):
        '''
        Inputs: Star Mass (Solar masses), luminousity (Solar luminousty) and temperature (Kelvin)
        Optionally: frequency spacing of the sun or it is set to the value of  135.1 micro-hertz ([[http://%28Frohlich%20et%20al%20%281997%29%29|Frohlich et al (1997)]])
        Outputs: Frequency difference between modes in micro-hertz
        '''

        if ((self.isoLumin==0) or (self.isoMass==0) or (self.isoTemp==0)):
            print 'Stellar parameters not set. Use <class_name>.populate_stellar_parameters().'
            return
        
        delta_nu = (((self.isoMass ** 0.5) * ((self.isoTemp / self.solar_temp) ** 3.)) / (self.isoLumin ** 0.75)) * self.solar_delta_nu        
        self.delta_nu = delta_nu
        if self.output==True:print 'Calculated delta_nu.'
        
        
    def calculate_A_max(self):
        '''
        Inputs: Star Mass (Solar masses), luminousity (Solar luminousity)
        Outputs: Amplitude in  meters/second
        '''
        
        if ((self.isoLumin==0) or (self.isoMass==0) or (self.isoTemp==0)):
            print 'Stellar parameters not set. Use <class_name>.populate_stellar_parameters().'
            return
        
        A_max = ((self.isoLumin ** 0.84) / (self.isoMass ** 1.32)) * self.solar_A_max
        self.A_max = A_max
        if self.output==True:print 'Calculated A_max in m/s'
        
        
    def calculate_nu_max(self):
        '''
        Inputs: Star Mass (Solar masses), luminousity (Solar luminousty) amd temperature (Kelvin)
        Optionally: Maximum frequency of the sun in micro-hertz, or it is set to 3090 micro-hertz ([[http://download.springer.com/static/pdf/639/art%253A10.1023%252FA%253A1004969622753.pdf?auth66=1355541218_4d08390331e483ae1e0f7db56e202198&ext=.pdf|Frohlich et al (1997)]])
        Outputs: The frequency at which the highest amplitude occurs in micro-hertz (nu_max)
        '''
        
        if ((self.isoLumin==0) or (self.isoMass==0) or (self.isoTemp==0)):
            print 'Stellar parameters not set. Use <class_name>.populate_stellar_parameters().'
            return
        
        nu_max = ((self.isoMass * ((self.isoTemp / self.solar_temp) ** 3.5)) / self.isoLumin) * (self.solar_nu_max)
        self.nu_max = nu_max
        if self.output==True:print 'Calculated nu_max.'
   
    def calculate_nu_amp(self, nu):
        """
        Inputs: The maximum amplitude of the star in question in meters/second, the maximum frequency in micro-hertz of the star and the frequency at which the amplitude is being calculated.
        Outputs: The velocity amplitude of the the star at the specified frequency
        """
        nu_amp = (self.A_max * math.exp((-16. * math.log(2.) * ((nu - self.nu_max) ** 2.)) / (self.nu_max ** 2.)))   
        return nu_amp
       
    def calculate_freq_modes(self, n_low=-1, n_high=-1, l_min=0, l_max=3):
        '''
         It iterates through all n and l values to return the frequency in mirco-hertz and amplitude in meters/second
    
        Parameters
        ----------
        n_low = -1: int
            Minimum value of the radial nodes.
            
        n_high = -1: int
            Maximum value of the radial nodes.
    
        l_min = 0: int
            Minimum number of surface nodes.
        
        l_max = 3: int
            Maximum number of surface nodes.
      
        Returns
        -------
        modefreq : np.array
            n x 2 np.array with frequency, amplitude
          
        Notes
        -----
        '''  
    
        
        if ((self.nu_max==0) or (self.delta_nu==0)):
            print 'Asterosesimological parameters not calculated.'
            print 'Use <class_name>.calculate_delta_nu()'
            print 'and <class_name>.calculate_nu_max().'
            return
       
        nulist = []
        amplist = []
        
        #If n_low or n_high is -1, then we compute this based on nu_max and D_nu
        if (n_low == -1):
               n_low = int(0.5*self.nu_max/self.delta_nu)
        if (n_high == -1):
               n_high = int(np.ceil(1.5*self.nu_max/self.delta_nu))
        if (l_max == 4):
               raise UserWarning
           
        for n in range(n_low,n_high):
            for l in range (l_min, l_max+1):
                nu = self.delta_nu*(n+(l/2.))
                A = self.calculate_nu_amp(nu)
                nulist.append(nu)
                amplist.append(A)

        self.modeFreq = nulist
        self.modeAmp = amplist
        if self.output==True:print 'Frequency modes calculated. See <class_name>.modeFreq and  <class_name>.modeAmp'

   
   
   
class ObsSimulation():
    
    outfile = 'out.csv'
    obs_samplingRate = 900  #How long will it take to produce a data point (in secs)
    obs_days = 4                     #Number of days that the observing run will last
    obs_error = 10                    #The amount of noise you wish to add to your timeseries, in meters per second
    obs_dailyWindow =  [1,-0.5,1,-0.5,1,-0.5,1,-0.5,1,-0.5,1,-16] #Sequence of intervals over an observing day (~ 24hrs). Negative means no observing. 
    iterations = 2                    #Number of times you would like ot run the script
    smoothness = 30
    
    
    def __init__(self):
        print 'ObsSimulation initialised'
        
 
    def create_noissyFlux(self):
        """
        Creates a time series from the input frequencies and amplitudes along with specified observational error
        """
        nu = np.multiply(self.star.modeFreq, 1e-6)
        A_vel = self.star.modeAmp

        self.noissyFlux = self.obs_error * np.random.normal(size=len(self.timeSequence)) #Adds noise due to equipment
        
        for i in range(len(nu)):
            As = (np.random.normal()/np.sqrt(2))*A_vel[i] #Adds intrinsic measurement noise.
            Ac = (np.random.normal()/np.sqrt(2))*A_vel[i]
            newMode = As * np.sin(nu[i] * 2 * np.pi * self.timeSequence) + Ac * np.cos(nu[i] * 2 * np.pi * self.timeSequence)
            self.noissyFlux = self.noissyFlux + newMode  
#             plt.plot(newMode)
#         plt.plot(self.cleanFlux) 
#         plt.show()
            

    def create_cleanFlux(self):
        """
        Creates a time series from the input frequencies and amplitudes
        """
        nu = np.multiply(self.star.modeFreq, 1e-6)
        A_vel = self.star.modeAmp
        
#         plt.ion()
        self.cleanFlux = np.zeros(len(self.timeSequence))
        for i in range(len(nu)):
            newMode = A_vel[i] * np.sin(nu[i] * 2 * np.pi * self.timeSequence) + A_vel[i] * np.cos(nu[i] * 2 * np.pi * self.timeSequence)
            self.cleanFlux = self.cleanFlux +  newMode
#             plt.plot(newMode)
#         plt.plot(self.cleanFlux) 
#         plt.show()
     
     
    def create_flux(self):
        
        self.timeSequence = np.arange(0,self.obs_days*24*3600, self.obs_samplingRate)                              #defines time domain
        self.create_cleanFlux()
        self.create_noissyFlux()        
        
#         plt.plot(self.cleanFlux)
#         plt.plot(self.noissyFlux)
#         plt.show()
    
    
    def create_lomb_scargle(self):
        
        t = self.timeSequence
        x = self.noissyFlux
        
        ts = pyTiming.pyPeriod.TimeSeries(t,x,error=np.ones(np.size(x))*1e-8)
        ft = pyTiming.pyPeriod.Fourier(ts)
        
        #Lets hard-wire a frequency separation of 0.5 cycles per data series length. For evenly sampled data,
        #this would mean that every second point is independent
        fsep = 0.5/(np.max(t)-np.min(t))
        self.ls = pyTiming.pyPeriod.Gls(ts,freq=np.arange(0,ts.returnNyquist()+fsep,fsep))
        self.noissyFlux_LS = self.ls.power


    def remove_sections(self):
        
        self.obs_dailyWindow = [1,-0.5,1,-0.5,1,-0.5,1,-0.5,1,-0.5,1,-16] #Sequence of intervals over 24hrs. Negative means no observing. 
        
        dayMap = []
        daily_window_sampled = np.multiply(self.obs_dailyWindow, 3600/self.obs_samplingRate)
        for i in daily_window_sampled:
            if i>0:
                dayMap = np.hstack((dayMap, np.repeat(True,abs(i))))
            elif i<0:
                dayMap = np.hstack((dayMap, np.repeat(False,abs(i))))
        
        fullMap = np.tile(dayMap, int(len(self.timeSequence)/len(dayMap)))
        fullMap = np.hstack((fullMap, dayMap[:len(self.timeSequence)-len(fullMap)]))
        
        self.cleanFlux = np.multiply(self.cleanFlux, fullMap.astype('int'))
        self.noissyFlux =  np.multiply(self.noissyFlux, fullMap.astype('int'))
        
    def plot(self):
        
        #This section plots the new time series for the non-continous observing, along with a descrete fourier transform that the new time series.
        plt.subplot(321)
        plt.title("Noissy Flux")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (m/s)")  
        plt.plot(self.timeSequence,self.noissyFlux)
        
        plt.subplot(322)
        plt.title("FFT of signal")
        plt.xlabel("Frequency ($\mu Hz$)")
        plt.ylabel("Amplitude")
        ftvel = np.abs(np.fft.fft(self.noissyFlux))
        plt.plot(ftvel)
        
        plt.subplot(323)
        plt.title("LS Noisy Flux")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (m/s)")  
        ftvel = np.fft.ifft(self.noissyFlux_LS)
        plt.plot(self.timeSequence, ftvel)
                
        plt.subplot(324)
        plt.title("FFT of signal")
        plt.xlabel("Frequency ($\mu Hz$)")
        plt.ylabel("Amplitude")
        plt.plot(self.noissyFlux_LS)

        plt.subplot(325)
        plt.title("LS Smoothed Noissy Flux")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude (m/s)")  
        ftvel = np.abs(np.fft.ifft(self.smoothNoissyFlux_LS))
        plt.plot(self.timeSequence, ftvel)
        
        plt.subplot(326)
        plt.title("FFT of signal")
        plt.xlabel("Frequency ($\mu Hz$)")
        plt.ylabel("Amplitude")
        plt.plot(self.smoothNoissyFlux_LS)
        
        plt.show()

    
              
    def check_star(self):
        
        self.goodStar = False  #A variable that it used it determine whether the star we are looking at is suitable
        if (np.max(self.smoothNoissyFlux_LS) > 2 * np.median(self.smoothNoissyFlux_LS)):    #determines whether the star used is a good star, by checking whether the peak power is abover a certain threshold. In this case the threshold is 2* the median amplitude of the Lomb-Scargle periodogram
            self.goodStar = True
       
    def moving_average(self, input_array = [], smoothness = []):
        '''
        Smoothing function. It takes in the time series and the size od the range of points to be averaged over.
        '''
        
        input_array = self.noissyFlux_LS
        smoothness = self.smoothness
        window = np.ones(int(smoothness))/float(smoothness)
        return np.convolve(input_array, window, 'same')    
                                       
    def calculate_std_dev(self):         
        '''
        stdList = [] #list of standard deviations for smoothNoissyFlux_LS
        for i in range(len(ls.freq)):    #This for loop calculated the standard deviation over range of 50 points in smoothNoissyFlux_LS
            binSet = []
            if i < 25:
                for j in range(50):
                    binSet.append(smoothNoissyFlux_LS[j])
            elif i > len(ls.freq)-25:
                this = len(ls.freq)-50
                for j in range(50):
                    binSet.append(smoothNoissyFlux_LS[this + j])
            else:
                for j in range(i-25, i+25):
                    binSet.append(smoothNoissyFlux_LS[j])
            std = np.std(binSet)
            stdList.append(std)
 
        plt.figure()
        fig1 = plt.subplot(111)
        plt.plot(ls.freq, stdList)
        plt.title("Standard deviation of each frequency")
        fig1.set_xlabel("Time (s)")
        fig1.set_ylabel("STD")
        
        plt.figure()
        fig2 = plt.subplot(111)
        plt.plot(self.ls.freq, self.smoothNoissyFlux_LS)
        plt.plot(self.ls.freq, np.ones(len(self.smoothNoissyFlux_LS)) * np.median(self.smoothNoissyFlux_LS), ':')
        plt.plot(self.ls.freq, np.ones(len(self.smoothNoissyFlux_LS)) * np.median(self.smoothNoissyFlux_LS)*2, ':')    
        plt.title("Smoothed Lomb-Scargle Periodogram")
        fig2.set_xlabel("Frequency (Hz)")
        fig2.set_ylabel("Power (probably (m/s)^2)")
        plt.show()
        '''
        
    def write_output(self):
 
        print "Mass(M_sun):", self.star.isoMass, ", Luminousity(L_sun):",self.star.isoLumin, ", Temperature(K):", self.star.isoTemp, "nu_max(microHz):", self.star.nu_max

        print "Peak frequency(microHz)(smooth_ls)    Peak frequency(microsHz)(smooth)    Peak amplitude (smooth_ls)    Peak amplitude (smooth)    Median Amlpitude (smooth)    Good Star?"
        
        for i in self.results:
            print i
     
    def save_results(self):
                
        Col1 = (self.ls.freq[np.where(self.noissyFlux_LS == max(self.noissyFlux_LS))][0]) * 1e6
        Col2 = (self.ls.freq[np.where(self.smoothNoissyFlux_LS == max(self.smoothNoissyFlux_LS))][0]) * 1e6
        Col3 = max(self.noissyFlux_LS)
        Col4 = max(self.smoothNoissyFlux_LS)
        Col5 = np.median(self.smoothNoissyFlux_LS)
        Col6 = self.goodStar
                
        if len(self.results)>0:
            self.results = np.append(self.results, [[ Col1, Col2, Col3, Col4, Col5, Col6]],0)
        else:
            self.results = np.array([[ Col1, Col2, Col3, Col4, Col5, Col6]])
                
    def clean_data(self):
        
        self.results = []
        self.star = []
        
    def start_simulation(self, mass = 1):
        
        self.clean_data()
        
        self.star = star()
        self.star.read_isochrone() 
        self.star.mass = mass
        self.star.populate_stellar_parameters()

        for i in range(self.iterations):            
            #This block creates and shapes the simulated signal
            self.create_flux()
            self.remove_sections()
            self.create_lomb_scargle()     
            self.smoothNoissyFlux_LS = self.moving_average() #smooths input function (noisy lomb-scargle by default).

            #complete data and output
            self.check_star()
            self.calculate_std_dev()
            self.save_results()    
   
   
# a = ObsSimulation()
# a.start_simulation()
# a.plot()
# a.write_output()



