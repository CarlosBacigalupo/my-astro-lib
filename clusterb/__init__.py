import numpy as np
import pylab as plt
import pandas as pd
import math
import sys

from params import *
import toolbox
import TableBrowser as TB

class webda():
    
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/ngc2477/webda/'
    
    def load_tables(self):
        
        main = np.loadtxt(self.base_dir + 'ad2000.coo', unpack=True, skiprows = 2)
        prob = np.loadtxt(self.base_dir + 'prob.mu', unpack=True, skiprows = 2)
        UBV = np.loadtxt(self.base_dir + 'ubv.ccd', unpack=True, skiprows = 2, usecols = [0,1,2,3])
    
        #These lines convert the np arrays into panda tables.
        #Select the format before cross reference here
        main_df = pd.DataFrame({ 'No' : main[0],
                                'RA' : main[2] + main[3]/60 + main[4]/3600,
                                'RA1' : main[2],
                                'RA2' : main[3],
                                'RA3' : main[4],
                                'Dec' : main[5] + main[6]/60 + main[7]/3600,
                                'Dec1' : main[5],
                                'Dec2' : main[6],
                                'Dec3' : main[7]})
    
        prob_df = pd.DataFrame({ 'No' : prob[0],
                                'Prob' : prob[2]})
    
        UBV_df = pd.DataFrame({ 'No' : UBV[0],
                                'V' : UBV[2],
                                'BV' : UBV[3]})
    
        #Creates the main merged table
#         main_df = main_df.merge(prob_df, on='No')
        main_df = main_df.merge(UBV_df, on='No')
        
        self.df = main_df
        
        #todo. Set constraints here.
        
        #Create .fld file
#         write_fld(main_df)
        
        print main_df.describe()
        a = main_df.describe()
        
        Dec_centre_dec = a['Dec']['mean']
        Dec_centre = toolbox.dec2sex(Dec_centre_dec)
        
        RA_centre_dec = a['RA']['mean']
        RA_centre = toolbox.dec2sex(RA_centre_dec)
        
        print Dec_centre_dec, RA_centre_dec
        print Dec_centre, RA_centre
        print 'end'

    
class vizier():
    
    base_dir = ''
    fileName = ''
    
    #stellar params
    tidal_rad = toolbox.sex2dec(0,60,00)
    cluster_centre = (0.,0.) #(toolbox.sex2dec(8,51,18) * 15, toolbox.sex2dec(11,48,0)  )  #Harris 97 (RA deg, Dec deg)
    cluster_centre_sex = ('08 51 18.0', '11 48 00.00') #Harris 97 (RA deg, Dec deg)
#     cluster_pm = ( 3.69, -0.39) #median yadav 2008 data pmRA * cos(Dec), pmDec 
#     cluster_pm = ( -2.4, -5.1) #median cone
    cluster_pm = (-9.6, -3.7) # pmRA * cos(Dec), pmDec (bellini+ 2010)
    
    #reduction params
    pp_cutoff = 10000 # based on observation (mas/yr)  
    mag_max = 17
    mag_min = 10

    def add_column(self, data, name):
        
        
        print 'Adding ' + name + ' column'   
        self.df[name] = data
        print ' Created ' + name + ' column. Unreduced df shape = ' + str(self.df.shape)
        print ''
        
        
    def add_Vmag_galah(self):
        #Tomaz's  galah eq.
        data = 2.0 * ( (self.df['Jmag'] - self.df['Kmag']) + 0.14) + 0.382 * np.exp( ((self.df['Jmag'] - self.df['Kmag']) -  0.2) /0.50) + self.df['Kmag'] 
        self.add_column(data, 'Vmag')


    def add_BVmag_ppmxl(self):

        data = self.df['b2mag'] - self.df['Vmag']
        self.add_column(data, 'BVmag')


    def add_BVmag_yadav(self):

        data = self.df['Bmag'] - self.df['Vmag']
        self.add_column(data, 'BVmag')


            
            
            
            
    def load_data(self):
        '''
        Loads a tab separated values output form vizier. project specific but most stuff can be reused. 
        '''

        main = np.genfromtxt(self.base_dir + self.fileName, dtype = str, delimiter = '\t', filling_values = '0')
        columns = main[0]
        main = main[3:]
        
        #dodgy cleanup of blank strings, 0 is the wrong value.
        main[main=='        '] = np.nan      
        main[main=='     '] = np.nan      
        main[main=='      '] = np.nan      
        main[main=='  '] = np.nan
        main[main==' '] = np.nan
        main[main==''] = np.nan
        
            
            
        self.df = TB.build_DataFrame(main.transpose(), columns)
        
        for i in self.df.columns:
            try:
                self.df[i] = self.df[i].astype(float)
                print 'Converted ', i
            except:
                print 'Failed ', i
                pass
        
        print str(main.shape[0]) + ' stars imported. Unreduced dataframe shape = ' + str(main.shape)
        print ''
                
        self.cluster_pm_median = (np.median(self.df['pmRA']), np.median(self.df['pmDE']))

    
        #Add Columns      
        data = self.df['pmRA'] - self.cluster_pm[0]
        self.add_column(data, 'delta_pmRA')
       
        data = self.df['pmDE'] - self.cluster_pm[1]
        self.add_column(data, 'delta_pmDE')
       
        data = np.sqrt((self.df['delta_pmRA'])**2 + (self.df['delta_pmDE'])**2)
        self.add_column(data, 'pm_tot')
       
        data = range(self.df.shape[0])
        self.add_column(data, 'ID')
                
        data = np.sqrt((self.cluster_centre[0] - self.df['RAJ2000'])**2 + (self.cluster_centre[1] - self.df['DEJ2000'])**2)
        self.add_column(data, 'd')
        
            
    def clear_reduced(self):
        
        self.reduced = []
    
    def reduce_list(self):
                
        self.reduced = self.df
        
        #Reduce
        print 'Removing NaNs based on BVmag.'
        self.reduced = self.reduced[-np.isnan(self.reduced['BVmag'])]
        print ' ' + str(self.reduced.shape[0]) + ' stars with non NaN BVmag.'
        print ''

        print 'Reducing list by magnitude range.'
        self.reduced = self.reduced[(self.reduced['Vmag'] < self.mag_max)]
        self.reduced = self.reduced[(self.reduced['Vmag'] > self.mag_min)]
        print ' ' + str(self.reduced.shape[0]) + ' stars with ' + str(self.mag_max) + ' > V_mag > ' +  str(self.mag_min)    

        print 'Reducing list by proper motion.'
        print ' Cluster\'s proper motion ' + str(self.cluster_pm) + ', median=' +str(self.cluster_pm_median)
        self.reduced = self.reduced[self.reduced['pm_tot'] < self.pp_cutoff]
        print ' ' + str(self.reduced.shape[0]) + ' stars with proper motion below ' + str(self.pp_cutoff)
        print ' '
        
        print 'Creating reduced list by tidal radius (Harris 1997)'
        print ' Using tidal_radius=' + str(self.tidal_rad) + ' deg'
        self.reduced = self.reduced[self.reduced['d'] < self.tidal_rad]
        print ' ' + str(self.reduced.shape[0]) + ' stars within ' + str(self.tidal_rad) + ' deg of cluster centre('+str(self.cluster_centre[0]) + ',' + str(self.cluster_centre[1]) + ')'
        print ''
        
        
        


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
        
       
    def plot_pm(self):
        x_pp = self.df['delta_pmDE']
        y_pp = self.df['delta_pmRA']
        x_pp_red = self.reduced['delta_pmDE']
        y_pp_red = self.reduced['delta_pmRA']
        x_pp_abs = self.df['pmDE']
        y_pp_abs = self.df['pmRA']
        plt.subplot(121, axisbg='black', title = "Cluster Corrected Proper Motions")
        plt.scatter(x_pp, y_pp, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_pp_red, y_pp_red, s=5, lw = 0, color = 'red', label='Selection')
        plt.xlabel('pmRA*cos(Dec) (mas/yr)')
        plt.ylabel('pmDec (mas/yr)')
        plt.legend()
        plt.subplot(122, axisbg='black', title = "Cluster Corrected Proper Motions")
        plt.scatter(x_pp_abs, y_pp_abs, s=5, lw = 0, color = 'blue', label='All stars')
        plt.show()


    def plot_HR_JK(self):
              
        x_mag = self.df['Jmag'] - self.df['Kmag']
        y_mag = self.df['Jmag']
        x_mag_red = self.reduced['Jmag'] - self.reduced['Kmag']
        y_mag_red = self.reduced['Jmag']
        
        plt.subplot(111, axisbg='black', title = "J-K vs J")
#         plt.scatter(x_mag, y_mag, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_mag_red, y_mag_red, s=5, lw = 0, color = 'white', label='Selection')
        
        plt.xlabel('J - K')
        plt.ylabel('J')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()

    def plot_HR(self):
              
        x_mag = self.df['BVmag']
        y_mag = self.df['Vmag']
        x_mag_red = self.reduced['BVmag']
        y_mag_red = self.reduced['Vmag']
        
        plt.subplot(111, axisbg='black', title = "B-V vs V")
#         plt.scatter(x_mag, y_mag, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_mag_red, y_mag_red, s=5, lw = 0, color = 'white', label='Selection')
        
        plt.xlabel('B - V')
        plt.ylabel('V')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()

    
    def plot_HR_pos_pm(self):
        
        
        x_mag = self.df['BVmag']
        y_mag = self.df['Vmag']
        x_mag_red = self.reduced['BVmag']
        y_mag_red = self.reduced['Vmag']
        plt.subplot(131, axisbg='black', title = "B-V vs V")
        plt.scatter(x_mag, y_mag, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_mag_red, y_mag_red, s=5, lw = 0, color = 'red', label='Selection')
        plt.xlabel('B - V')
        plt.ylabel('V')
        plt.legend()
        plt.gca().invert_yaxis()

        x_RA = self.df['RAJ2000']
        y_Dec = self.df['DEJ2000']
        x_RA_red = self.reduced['RAJ2000']
        y_Dec_red = self.reduced['DEJ2000']
        plt.subplot(132, axisbg='black', title = "RA vs Dec")
        plt.scatter(x_RA, y_Dec, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_RA_red, y_Dec_red, s=5, lw = 0, color = 'red', label='Selection')
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        plt.legend()

        x_pp = self.df['delta_pmDE']
        y_pp = self.df['delta_pmRA']
        x_pp_red = self.reduced['delta_pmDE']
        y_pp_red = self.reduced['delta_pmRA']
        plt.subplot(133, axisbg='black', title = "Cluster Corrected Proper Motions")
        plt.scatter(x_pp, y_pp, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_pp_red, y_pp_red, s=5, lw = 0, color = 'red', label='Selection')
        plt.xlabel('pmRA*cos(Dec) (mas/yr)')
        plt.ylabel('pmDec (mas/yr)')
        plt.legend()        
        plt.show()

    def plot(self, main_table, orig_table):
        
        
        x_all = orig_table[:,4]-orig_table[:,6]
        y_all = - orig_table[:,4] 
        x_tidal = main_table[:,4]-main_table[:,6]
        y_tidal = - main_table[:,4] 
        fig =plt.figure(facecolor='grey') 
        ax =plt.subplot(111, axisbg='black', title = "J vs. J-K")
        plt.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
        plt.xlabel('J-K (mag)')
        plt.ylabel('J (mag)')
        ax.legend()
        plt.show()
        
        x_all = orig_table[:,4]-orig_table[:,6]
        y_all = - orig_table[:,10] 
        x_tidal = main_table[:,4]-main_table[:,6]
        y_tidal = - main_table[:,10] 
        fig =plt.figure(facecolor='grey') 
        ax =plt.subplot(111, axisbg='black', title = "V vs. J-K")
        plt.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
        plt.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
        plt.xlabel('J-K (mag)')
        plt.ylabel('V (mag)')
        ax.legend()
        plt.show()

        
    def create_fld(self):   
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
    main_temp[:,-1] = np.sqrt((self.cluster_pm[0] - main_temp[:,2])**2 + (self.cluster_pm[1] - main_temp[:,3])**2)
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
    fig =plt.figure(facecolor='grey') 
    ax =plt.subplot(111, axisbg='black', title = "Proper Motions")
    plt.scatter(x_pp, y_pp, s=5, lw = 0, color = 'blue', label='All stars')
    plt.xlabel('pmRA*cos(Dec) (mas/yr)')
    plt.ylabel('pmDec (mas/yr)')
    ax.legend()
    #    plt.show()
    
    x_RA = main_table[:,0]
    y_Dec = main_table[:,1]
    fig =plt.figure(facecolor='grey') 
    ax =plt.subplot(111, axisbg='black', title = "RA vs Dec")
    plt.scatter(x_RA, y_Dec, s=5, lw = 0, color = 'blue', label='All stars')
    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    ax.legend()
    #    plt.show()
    
    x_all = orig_table[:,4]-orig_table[:,6]
    y_all = - orig_table[:,4] 
    x_tidal = main_table[:,4]-main_table[:,6]
    y_tidal = - main_table[:,4] 
    fig =plt.figure(facecolor='grey') 
    ax =plt.subplot(111, axisbg='black', title = "J vs. J-K")
    plt.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
    plt.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
    plt.xlabel('J-K (mag)')
    plt.ylabel('J (mag)')
    ax.legend()
    #    plt.show()
    
    x_all = orig_table[:,4]-orig_table[:,6]
    y_all = - orig_table[:,10] 
    x_tidal = main_table[:,4]-main_table[:,6]
    y_tidal = - main_table[:,10] 
    fig =plt.figure(facecolor='grey') 
    ax =plt.subplot(111, axisbg='black', title = "V vs. J-K")
    plt.scatter(x_all, y_all, s=5, lw = 0, color = 'blue', label='All stars')
    plt.scatter(x_tidal, y_tidal, s=5, lw = 0, color = 'red', label='Stars selected so far')
    plt.xlabel('J-K (mag)')
    plt.ylabel('V (mag)')
    ax.legend()
    plt.show()

def load_tables():
    global DATA_DIR
    
    main = np.loadtxt(DATA_DIR + 'ad2000.coo', unpack=True, skiprows = 2)
    prob = np.loadtxt(DATA_DIR + 'prob.mu', unpack=True, skiprows = 2)
    UBV = np.loadtxt(DATA_DIR + 'ascc.mes', unpack=True, skiprows = 2)

    #These lines convert the np arrays into panda tables.
    #Select the format before cross reference here
    main_df = pd.DataFrame({ 'No' : main[0],
                            'RA' : main[2] + main[3]/60 + main[4]/3600,
                            'RA1' : main[2],
                            'RA2' : main[3],
                            'RA3' : main[4],
                            'Dec' : main[5] + main[6]/60 + main[7]/3600,
                            'Dec1' : main[5],
                            'Dec2' : main[6],
                            'Dec3' : main[7]})

    prob_df = pd.DataFrame({ 'No' : prob[0],
                            'Prob' : prob[2]})

    UBV_df = pd.DataFrame({ 'No' : UBV[0],
                            'V' : UBV[2],
                            'BV' : UBV[3]})

    #Creates the main merged table
    main_df = main_df.merge(prob_df, on='No')
    main_df = main_df.merge(UBV_df, on='No')
    
    
    #todo. Set constraints here.
    
    #Create .fld file
    write_fld(main_df)
    
    print main_df.describe()
    a = main_df.describe()
    
    Dec_centre_dec = a['Dec']['mean']
    Dec_centre = dec2sex(Dec_centre_dec)
    
    RA_centre_dec = a['RA']['mean']
    RA_centre = dec2sex(RA_centre_dec)
    
    print Dec_centre_dec, RA_centre_dec
    print Dec_centre, RA_centre
    print 'end'

def plotter():
    plt.ion()
    
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
    
    
    fig =plt.figure(facecolor='grey') 
    plt.draw()
    ax =plt.subplot(111, axisbg='black')
    plt.scatter(x, y, s=2*mag_norm, lw = 0, color = colormap)
    plt.draw()
    #fig = figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    #scatter(x, y, )
    
    plt.show()
    
    print 'end'
    


