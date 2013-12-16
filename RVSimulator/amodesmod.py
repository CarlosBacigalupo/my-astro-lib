#!/usr/bin/env python

import sys
import csv
import math
import numpy

"""
Reads in 4 values in the following order: Mass, Luminousity, Temperature and the name of the output file.
The Mass and Luminousity are given in solar units, and temperature is given in Kelvin, and the output file is in the format .csv
"""

D_nu = 0
A_max = 0
nu_max = 0
n_low = 0
n_high = 0

def Delta_nu(mass, lumin, temp, delta_sun=135.1):
    """
    Inputs: Star Mass (Solar masses), luminousity (Solar luminousty) and temperature (Kelvin)
        Optionally: frequency spacing of the sun or it is set to the value of  135.1 micro-hertz ([[http://%28Frohlich%20et%20al%20%281997%29%29|Frohlich et al (1997)]])
    Outputs: Frequency difference between modes in micro-hertz
    """
    return (((mass ** 0.5) * ((temp / 5777.) ** 3.)) / (lumin ** 0.75)) * delta_sun

def vel_Amp(mass, lumin, lumin_exp=0.84, mass_exp=1.32, Amp_velsun=0.2):
    """
    Inputs: Star Mass (Solar masses), luminousity (Solar luminousity)
    Optionally: Luminousty exponent and  mass exponent
    Outputs: Amplitude in  meters/second
    """
    return ((lumin ** lumin_exp) / (mass ** mass_exp)) * Amp_velsun

def find_nu_max(mass, lumin, temp, nu_max_sun=3090.):
    """
    Inputs: Star Mass (Solar masses), luminousity (Solar luminousty) amd temperature (Kelvin)
    Optionally: Maximum frequency of the sun in micro-hertz, or it is set to 3090 micro-hertz ([[http://download.springer.com/static/pdf/639/art%253A10.1023%252FA%253A1004969622753.pdf?auth66=1355541218_4d08390331e483ae1e0f7db56e202198&ext=.pdf|Frohlich et al (1997)]])
    Outputs: The frequency at which the highest amplitude occurs in micro-hertz (nu_max)
    """
    return ((mass * ((temp / 5777.) ** 3.5)) / lumin) * (nu_max_sun)

def setconst(mass, lumin, temp):
    """
    Sets the constants mass, luminousity adn temperature.
    """
    global D_nu
    global A_max
    global nu_max
    D_nu = Delta_nu(mass, lumin, temp)
    A_max = vel_Amp(mass, lumin)
    nu_max = find_nu_max(mass, lumin, temp)
    print "D_nu (microhertz):" + str(D_nu) + " A_max(m/s): " + str(A_max) + " nu_max(microhertz): " + str(nu_max)


def nu_amp(A_max, nu_max, nu):
    """
    Inputs: The maximum amplitude of the star in question in meters/second, the maximum frequency in micro-hertz of the star and the frequency at which the amplitude is being calculated.
    Outputs: The velocity amplitude of the the star at the specified frequency
    """
    return (A_max * math.exp((-16. * math.log(2.) * ((nu - nu_max) ** 2.)) / (nu_max ** 2.)))


def modefreq(n_low=-1, n_high=-1, l_min=0, l_max=3):
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

    nulist = []
    amplist = []
    
    #If n_low or n_high is -1, then we compute this based on nu_max and D_nu
    if (n_low == -1):
	       n_low = int(0.5*nu_max/D_nu)
    if (n_high == -1):
	       n_high = int(numpy.ceil(1.5*nu_max/D_nu))
    if (l_max == 4):
	       raise UserWarning
       
    for n in range(n_low,n_high):
        for l in range (l_min, l_max+1):
            nu = D_nu*(n+(l/2.))
            A = nu_amp(A_max, nu_max, nu)
            nulist.append(nu)
            amplist.append(A)
            
    modefreq = numpy.array([nulist, amplist])
    
    return modefreq


def return_nu_max():
    """
    Inputs: nothing.
    Outputs: peak frequency.
    """
    return nu_max

