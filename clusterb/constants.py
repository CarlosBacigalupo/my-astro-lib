SOLAR_MASS = 1.9891e30 #Kg.
NU_MAX_SUN = 3090. # +/- 30 microHz   (Huber at al. 2011)
TEMP_EFF_SUN = 5777. #K
SOLAR_AMAG = 4.83 #Mag
SOLAR_LUMIN = 3.839e26 #W J/s kg*m^2/s^3

header_fixed = ['ID','Parallax','ParallaxSig','CI','CISig','vMag','vMagSig']

header_add = ['d','Lumin','Temp','Mass','nuMax','PMax','PMaxHrs','PMaxDays']

header = header_fixed + header_add

#the lines bellow shoud be redundant with header (not yet)
#from dataset
IDCol = 0
parallaxCol = 1
paralaxSigmaCol = 2
CICol = 3
CIColSigma = 4
vMagCol = 5
vMagSigmaCol = 6

#added columns
dCol = 7
luminCol = 8
tempCol = 9
massCol = 10
nuMaxCol = 11
PMaxCol = 12
PMaxHrsCol = 13
PMaxDaysCol = 14
