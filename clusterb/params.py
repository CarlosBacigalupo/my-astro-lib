from tools import  *


DATA_DIR = 'targets/47_tuc/'



tidal_rad = sex2dec(0,42,19)
cluster_centre = (sex2dec(0,24,05.2)* 15, - sex2dec(72,04,51)) #Harris 97 (RA deg, Dec deg)
cluster_centre_sex = ('00 24 05.20', '-72 04 51.00') #Harris 97 (RA deg, Dec deg)

cluster_pm = (1.730, -2.733) # pmRA * cos(Dec), pmDec http://www.astro.yale.edu/dana/gl_2012_J2000.cat1 

pp_cutoff = 20 # based on observation (mas/yr)
pp_cluster = -18 # harris 97 (mas/yr)
mag_max = 13
mag_min = 12
