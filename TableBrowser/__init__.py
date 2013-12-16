import numpy as np
import pylab as p
import pandas as pd
import math

# from params import *
# from tools import *


def load_47tuc():
    global DATA_DIR, tidal_rad, pp_cutoff, cluster_centre_sex, mag_max, mag_min, cluster_pm
    
    '''Fields
      0-RAJ2000     1-DEJ2000     2-pmRA      3-pmDE     4-Jmag    5-Hmag    6-Kmag   7-b1mag  8-b2mag   9-imag 10-Vmag 11-pm_tot 12-ID 13-R(dist from centre)
    deg         deg      mas/yr    mas/yr     mag     mag     mag    mag    mag    mag 
    '''
    main = np.loadtxt(DATA_DIR + 'all.tsv')
    
    print str(main.shape[0]) + ' stars imported.' + str(main.shape)
    print ''
    
    
    ###Clean
    print 'Cleaning 99.999 data placeholders, this is probably wrong'
    main[main[:,4]==99.999,4]=0
    main[main[:,5]==99.999,5]=0
    main[main[:,6]==99.999,6]=0
    main[main[:,7]==99.99,7]=0
    main[main[:,8]==99.99,8]=0
    main[main[:,9]==99.99,9]=0
    
    
    
    ###Add Columns
    print 'Adding V_mag column'
    main_temp = np.zeros((main.shape[0],main.shape[1]+1))
    main_temp[:,:-1] = main
    main_temp[:,-1] = 2.0 * ( (main_temp[:,4] - main_temp[:,6]) + 0.14) + 0.382 * np.exp( ((main_temp[:,4] - main_temp[:,6]) -  0.2) /0.50) + main_temp[:,6]
    main = main_temp    
    print ' Created V_mag ' + str(main_temp.shape)
    print ''
    
    print 'Adding pm_tot column'
    main_temp = np.zeros((main.shape[0],main.shape[1]+1))
    main_temp[:,:-1] = main
    main_temp[:,-1] = np.sqrt((cluster_pm[0] - main_temp[:,2])**2 + (cluster_pm[1] - main_temp[:,3])**2)
    main = main_temp    
    print ' Created pm_tot ' + str(main_temp.shape)
    print ''
    
    print 'Adding identifier column to create fld file'
    main_temp = np.zeros((main.shape[0],main.shape[1]+1))
    main_temp[:,:-1] = main
    main_temp[:,-1] = range(main_temp.shape[0])
    main = main_temp    
    print ' Created main_dist ' + str(main_temp.shape)
    print ''
    
    print 'Adding distance from centre column'
    main_dist = np.zeros((main.shape[0],main.shape[1]+1))
    main_dist[:,:-1] = main
    print ' Created main_dist ' + str(main_dist.shape)
    main_dist[:,-1] = np.sqrt((cluster_centre[0] - main_dist[:,0])**2 + (cluster_centre[1] - main_dist[:,1])**2)
    print ' Added distances ' + str(main_dist[:,-1].shape)
    print ''

    
    #Reduce
    print 'Creating reduced list by tidal radius (Harris 1997)'
    print ' Using tidal_radius=' + str(tidal_rad) + ' deg'
    main_tidal = main_dist[main_dist[:,-1] < tidal_rad]
    print ' ' + str(main_tidal.shape[0]) + ' stars within ' + str(tidal_rad) + ' deg of cluster centre('+str(cluster_centre[0]) + ',' + str(cluster_centre[1]) + ')'
    print ''
#     plotter_47_tuc(main_tidal, main)
    
    print 'Reducing list by proper motion'
    print ' Cluster\'s proper motion ' + str(cluster_pm)
    print ' Using proper motion cutoff=' + str(pp_cutoff) +' (mas/yr)' 
    main_tidal_pp = main_tidal[main_tidal[:,11] < pp_cutoff]
    print ' ' + str(main_tidal_pp.shape[0]) + ' stars with proper motion below ' + str(pp_cutoff)
    print ' '
#     plotter_47_tuc(main_tidal_pp, main)
    
#     print 'Creating density curve'
#     step = 0.05 #steps to average density arcmins
#     all_steps = range(int(1 /step)+1)
#     density = np.zeros(len(all_steps))    
#     for i in all_steps:
#         this_donut = main_tidal_pp[main_tidal_pp[:,-1] > step*i]
#         if len(this_donut)>0:
#             this_donut = this_donut[this_donut[:,-1] < step*(i+1)]
#             density[i] = this_donut.shape[0]/(math.pi*(step*(i+1))**2 - math.pi*(step*i)**2)
#             print 'Between ' + str(step*i) + ' and '+ str(step*(i+1)) + ' deg from cluster centre'
#             print 'area ' +str (i) + ' = '+ str(math.pi*(step*(i+1))**2 - math.pi*(step*i)**2) + ' deg**2'
#             print 'density = ' + str(density[i]) + ' stars per deg**2'
#             print 
    
    print 'Reducing list by magnitude range'
    print ' Using ' + str(mag_max) + ' > V_mag (imag field) > ' +  str(mag_min) 
    main_tidal_pp_mag = main_tidal_pp[(main_tidal_pp[:,10] < mag_max)]
    main_tidal_pp_mag = main_tidal_pp_mag[(main_tidal_pp_mag[:,10] > mag_min)]
    print ' ' + str(main_tidal_pp_mag.shape[0]) + ' stars with ' + str(mag_max) + ' > V_mag_2mass > ' +  str(mag_min)    
    plotter_47_tuc(main_tidal_pp_mag, main)
    
    
    ###Prepare Sky
    sky_out = main
    
    
    ###Create panda tables
    #These lines convert the np arrays into panda tables.
    #This is the final tabel before export to fld
    main_df = pd.DataFrame({ 'No' : main_tidal_pp_mag[:,12],
                            'RA' : main_tidal_pp_mag[:,0]/15 ,
                            'RA1' : dec2sex(main_tidal_pp_mag[:,0]/15)[0],
                            'RA2' : dec2sex(main_tidal_pp_mag[:,0]/15)[1],
                            'RA3' : dec2sex(main_tidal_pp_mag[:,0]/15)[2],
                            'Dec' : main_tidal_pp_mag[:,1],
                            'Dec1' : dec2sex(main_tidal_pp_mag[:,1])[0],
                            'Dec2' : dec2sex(main_tidal_pp_mag[:,1])[1],
                            'Dec3' : dec2sex(main_tidal_pp_mag[:,1])[2],
                            'V' : main_tidal_pp_mag[:,10]})
    
    keep = np.asarray(main_df['No'])
    
    #Sky
    sky_df = pd.DataFrame({ 'No' : sky_out[:,12],
                            'RA' : sky_out[:,0]/15 ,
                            'RA1' : dec2sex(sky_out[:,0]/15)[0],
                            'RA2' : dec2sex(sky_out[:,0]/15)[1],
                            'RA3' : dec2sex(sky_out[:,0]/15)[2],
                            'Dec' : sky_out[:,1],
                            'Dec1' : dec2sex(sky_out[:,1])[0],
                            'Dec2' : dec2sex(sky_out[:,1])[1],
                            'Dec3' : dec2sex(sky_out[:,1])[2],
                            'V' : sky_out[:,10]})
    sky_df = np.array([])
#     sky_df = sky_df[~sky_df['No'].isin(keep)]
    
    #Fiducials
    main_fid = main[main[:,10]>5]
    main_fid = main_fid[main_fid[:,10]<10.5]
    
    fiducials_df = pd.DataFrame({ 'No' : main_fid[:,12],
                            'RA' : main_fid[:,0]/15,
                            'RA1' : dec2sex(main_fid[:,0]/15)[0],
                            'RA2' : dec2sex(main_fid[:,0]/15)[1],
                            'RA3' : dec2sex(main_fid[:,0]/15)[2],
                            'Dec' : main_fid[:,1],
                            'Dec1' : dec2sex(main_fid[:,1])[0],
                            'Dec2' : dec2sex(main_fid[:,1])[1],
                            'Dec3' : dec2sex(main_fid[:,1])[2],
                            'V' : main_fid[:,10]})
      
    fiducials_df = fiducials_df.sort(['V']).tail(80)
    
    #remove target duplicates
    fiducials_df = fiducials_df[~fiducials_df['No'].isin(keep)]
    
    write_fld(main_df, sky_df, fiducials_df, DATA_DIR + 'output.fld', cluster_centre_sex)
        
def plotter_47_tuc(main_table, orig_table):
    
    x_pp = main_table[:,2]
    y_pp = main_table[:,3]
    fig = p.figure(facecolor='grey') 
    ax = p.subplot(111, axisbg='black', title = "Proper Motions")
    p.scatter(x_pp, y_pp, s=5, lw = 0, color = 'blue', label='All stars')
    p.xlabel('pmRA*cos(Dec) (mas/yr)')
    p.ylabel('pmDec (mas/yr)')
    ax.legend()
#     p.show()
 
    x_RA = main_table[:,0]
    y_Dec = main_table[:,1]
    fig = p.figure(facecolor='grey') 
    ax = p.subplot(111, axisbg='black', title = "RA vs Dec")
    p.scatter(x_RA, y_Dec, s=5, lw = 0, color = 'blue', label='All stars')
    p.xlabel('RA (deg)')
    p.ylabel('Dec (deg)')
    ax.legend()
#     p.show()
    
    x_all = orig_table[:,4]-orig_table[:,6]
    y_all = - orig_table[:,4] 
    x_tidal = main_table[:,4]-main_table[:,6]
    y_tidal = - main_table[:,4] 
    fig = p.figure(facecolor='grey') 
    ax = p.subplot(111, axisbg='black', title = "J vs. J-K")
    p.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
    p.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
    p.xlabel('J-K (mag)')
    p.ylabel('J (mag)')
    ax.legend()
#     p.show()
    
    x_all = orig_table[:,4]-orig_table[:,6]
    y_all = - orig_table[:,10] 
    x_tidal = main_table[:,4]-main_table[:,6]
    y_tidal = - main_table[:,10] 
    fig = p.figure(facecolor='grey') 
    ax = p.subplot(111, axisbg='black', title = "V vs. J-K")
    p.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
    p.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
    p.xlabel('J-K (mag)')
    p.ylabel('V (mag)')
    ax.legend()
    p.show()

def read_external_file(fileName, skipRows = 0):
        tableData = np.loadtxt(fileName, unpack=True, skiprows = skipRows)
        
        return tableData


def build_DataFrame(d, tableData):
        
    df_dict = {}
    
    for i in range(len(d)):
        df_dict[d[i]] = np.array(tableData.field(d[i]))
        
    df = pd.DataFrame(df_dict)

    return df


#     
#     main_df = pd.DataFrame({ 'No' : main_tidal_pp_mag[:,12],
#                             'RA' : main_tidal_pp_mag[:,0]/15 ,
#                             'RA1' : dec2sex(main_tidal_pp_mag[:,0]/15)[0],
#                             'RA2' : dec2sex(main_tidal_pp_mag[:,0]/15)[1],
#                             'RA3' : dec2sex(main_tidal_pp_mag[:,0]/15)[2],
#                             'Dec' : main_tidal_pp_mag[:,1],
#                             'Dec1' : dec2sex(main_tidal_pp_mag[:,1])[0],
#                             'Dec2' : dec2sex(main_tidal_pp_mag[:,1])[1],
#                             'Dec3' : dec2sex(main_tidal_pp_mag[:,1])[2],
#                             'V' : main_tidal_pp_mag[:,10]})
    
    
    
    
#     #Creates the main merged table
#     main_df = main_df.merge(prob_df, on='No')
#     main_df = main_df.merge(UBV_df, on='No')
#     
#     
#     #todo. Set constraints her

def plotter():
    p.ion()
    
    DATA_DIR = ''
    
    dataFile = 'NGC2477_XYpos.txt'
    
    a = np.loadtxt(DATA_DIR + dataFile, unpack=True)
    
    starID = a[0]
    mag = a[1]
    x = a[2]
    y = a[3]
    
    
    mag_norm = abs(mag/max(mag))
    mag_norm[mag_norm>=1]=1
    rgb = np.array(mag_norm) * np.ones((3, len(mag_norm)))
    colormap = rgb.transpose()
    #colormap = colormap**(10)
    
    
    fig = p.figure(facecolor='grey') 
    p.draw()
    ax = p.subplot(111, axisbg='black')
    p.scatter(x, y, s=2*mag_norm, lw = 0, color = colormap)
    p.draw()
    #fig = figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    #scatter(x, y, )
    
    p.show()
    
    print 'end'
    


