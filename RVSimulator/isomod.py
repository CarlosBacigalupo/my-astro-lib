#!/usr/bin/env python
import sys
import numpy


#defines variables
the_infile = ""
the_value = 0
luminthreshold = 1  #checked if luminousity of the star is above this value. This is to make sure the module returns a giant star, not a main sequence one.


#global mass, lumin, temp, matchlumin, matchtemp

matchlumin = numpy.empty(0)
matchtemp = numpy.empty(0)


def readmod(infile):            #Reads in model data
	"""
	Readmod reads in the isochrone model in order ot be able to later determine the luminousity and temperature of the input mass in the model.py program. 
	The input for this function is the isochrone in .csv format.
	"""
	a = numpy.loadtxt(infile,delimiter=',',skiprows=13)
	mass = a[:,2]
	lumin=a[:,4]
	temp=a[:,5]
	return mass, lumin, temp

def find(value, mass, lumin, temp):
	"""
	This function finds the nearest mass match in the isochrone to the input mass. It also checks whether the luminousty is above the declared luminousity threshold.
		This threshold is to make sure the returned values correspond to a giant star and not a main sequence star.
	The input for this function is a mass to be found and the isochrones list of mass, luminousity and temperature.
	"""
	matchlumin = numpy.empty(0)
	matchtemp = numpy.empty(0)
	nearest = find_nearest(mass, value)
	matchlumin = numpy.append(matchlumin, lumin[nearest])
	matchtemp = numpy.append(matchtemp, temp[nearest])
	if lumin[nearest] < luminthreshold:
		print "Not a giant"
	return matchlumin, matchtemp

def set_data(the_file, the_mass):
	"""
	This function sets the data in this module in order to make running all other functions in the module simple. 
	"""
	global the_infile
	global the_value
	the_infile = the_file
	the_value = the_mass


def find_nearest(array,value):
	"""
	Finds the index of the nearest mass to that input value.
	"""
	idx = (numpy.abs(array-value)).argmin()
	return idx


def run():
	"""
	Runs module to find and return the luminousity and temperature of input mass.
	"""
	mass, lumin, temp=readmod(the_infile)
	matchlumin, matchtemp = find(the_value, mass, lumin, temp)
 	return matchlumin,matchtemp

