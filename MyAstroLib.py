import math
import numpy as np
import scipy.constants as const
from scipy import signal
import pylab as plt
import matplotlib.cm as cm    

import RVSimulator as RVS
import tools

def calculate_dRV_HERMES():
#Calculates the photon noise limited RV precision based on a given spectrum
#Fibre type: F-guide, N-broken, dead or no fibre, 
#            P-program (science), S - Sky, U-unallocated or unused
    
    
    fibreDict = {'F':'guide', 'N':'broken, dead or no fibre', 'P' :'program (science)', 'S': 'Sky', 'U':'unallocated or unused'}
    c = const.c  #[m/s]
    
    
    fileName1='24oct1_combined.fits'
    fileName2='24oct2_combined.fits'
    fileName3='24oct3_combined.fits'
    fileName4='24oct4_combined.fits'
    
    fibreTable = read_fibre_table(fileName1)
    
    fileData1 = RVS.extract_all_from_2dfdr(fileName1)
    fileData2 = RVS.extract_all_from_2dfdr(fileName2)
    fileData3 = RVS.extract_all_from_2dfdr(fileName3)
    fileData4 = RVS.extract_all_from_2dfdr(fileName4)
    
    Lambda1 = RVS.extract_HERMES_wavelength(fileName1)
    Lambda2 = RVS.extract_HERMES_wavelength(fileName2)
    Lambda3 = RVS.extract_HERMES_wavelength(fileName3)
    Lambda4 = RVS.extract_HERMES_wavelength(fileName4)
    
    for i in range(399):
#         print 'Reading fibre ' + str(i+1)
        A01 = fileData1[i]
        A02 = fileData2[i]
        A03 = fileData3[i]
        A04 = fileData4[i]

        A01[np.isnan(A01)]=0
        A02[np.isnan(A02)]=0
        A03[np.isnan(A03)]=0
        A04[np.isnan(A04)]=0   
        
        Lambda = np.hstack((Lambda1[A01>0], Lambda2[A02>0], Lambda3[A03>0], Lambda4[A04>0]))
        
        A01 = A01[A01>0]
        A02 = A02[A02>0]
        A03 = A03[A03>0]
        A04 = A04[A04>0]
        
        A0 = np.hstack((A01, A02, A03, A04))

        if np.sum(A01)>0:        
            W1 = RVS.W(Lambda1[A01>0], A01)
            W1[np.isnan(W1)]=0           
            
        if np.sum(A02)>0:        
            W2 = RVS.W(Lambda2[A02>0], A02)
            W2[np.isnan(W2)]=0           

        if np.sum(A03)>0:        
            W3 = RVS.W(Lambda3[A03>0], A03)
            W3[np.isnan(W3)]=0           

        if np.sum(A04)>0:        
            W4 = RVS.W(Lambda4[A04>0], A04)
            W4[np.isnan(W4)]=0           

        if np.sum(A0)>0:        
            W = RVS.W(Lambda, A0)
            W[np.isnan(W)]=0           
           
            if np.sum(W)!=0:            
                Q1 = RVS.Q(W1, A01)
                Q2 = RVS.Q(W2, A02)
                Q3 = RVS.Q(W3, A03)
                Q4 = RVS.Q(W4, A04)
                Q = RVS.Q(W, A0)
                
                dRV1 = c/math.sqrt(np.sum(W1))
                dRV2 = c/math.sqrt(np.sum(W2))
                dRV3 = c/math.sqrt(np.sum(W3))
                dRV4 = c/math.sqrt(np.sum(W4))
                dRV = c/math.sqrt(np.sum(W))
    
                print i+1, fibreDict[str(fibreTable[i]['TYPE']).strip()], str(fibreTable[i]['NAME']).strip(),str(fibreTable[i]['MAGNITUDE']).strip() 
                print '   ' +   str(Q1), str(dRV1)
                print '   ' +   str(Q2), str(dRV2)
                print '   ' +   str(Q3), str(dRV3)
                print '   ' +   str(Q4), str(dRV4)
                print '   ' +   str(Q), str(dRV)
                print ''               
    # plt.plot(Lambda, A0)
    # plt.show()            

def create_fibre_table_DataFrame():
    
    import numpy as np
    import RVSimulator as RVS
    import TableBrowser as TB
    
    fileName='24oct3_combined.fits'
    
    d = np.array(['NAME', 'RA', 'DEC', 'X', 'Y', 'XERR', 'YERR', 'THETA', 'TYPE', 'PIVOT', 'MAGNITUDE', 'PID', 'COMMENT', 'RETRACTOR', 'WLEN', 'PMRA', 'PMDEC'])
    
    tableData = RVS.load_2dfdr_fibre_table(fileName)
    
    df = TB.build_DataFrame(d, tableData)

    return df




# find_PSF_arc()
# find_PSF_flat()
#read_psf_data_spatial()
# read_psf_data_spectral()

# do_all_plots()

def calculate_RV_shift_HERMES(fibres, files):
    import pyfits as pf
    
    c = const.c
    
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/RV/'    
    JD = np.zeros(len(files))
    Lambda = []
    fileData = []
    fibreTable = []
    
    for i in range(len(files)):    
        a = pf.open(WORKING_DIR + files[i])
        JD[i] = a[0].header['UTMJD']
        
        temp = RVS.extract_HERMES_wavelength(WORKING_DIR + files[i])
        if (len(fileData)==0):
            Lambda = np.array([temp])
        else:
            Lambda = np.append(Lambda, [temp], 0)

        temp = RVS.read_fibre_table(WORKING_DIR + files[i])
        if (len(fileData)==0):
            fibreTable = np.array([temp])
        else:
            fibreTable = np.append(fibreTable, [temp], 0)

        temp = RVS.extract_all_from_2dfdr(WORKING_DIR + files[i])
        if (len(fileData)==0):
            fileData = np.array([temp])
        else:
            fileData = np.append(fileData, [temp],0)
    
    print ' Read ' + str(len(files)) + ' files'
    
    for j in range(399):     
        for i in range(fileData.shape[0]-1):
            fibre1 = fibreTable[i,j]
            fibre2 = fibreTable[i+1,j]
            # fibre types: F-guide, N-broken, dead or no fibre, P-program (science), S - Sky, U-unallocated or unused
            if ((fibre1['NAME'].strip()==fibre2['NAME'].strip()) & (fibre1['TYPE']=='P') & (fibre1['NAME'].strip().lower()!='parked')):
                time1 = JD[i]
                time2 = JD[i+1]
                lambda1 = Lambda[i]
                lambda2 = Lambda[i+1]
                flux1 = fileData[i,j,:]
                flux2 = fileData[i+1,j,:]
                newLambda1, flux1Neat = RVS.clean_flux(flux1, 10, lambda1 )
                newLambda2, flux2Neat = RVS.clean_flux(flux2, 10, lambda2)                   
                
                
                a =  np.correlate(flux1Neat, flux2Neat, 'same')
                b =  np.correlate(flux1, flux2, 'same')
                print ''
                print 'Time Difference: ' + str(abs(time1 - time2)*24*60) + ' min'
                deltaLambda = newLambda2[len(newLambda2)/2] - newLambda2[np.where(a==max(a))[0][0]]
                print 'Wavelength difference: ' + str(deltaLambda) + ' Ang'
                print 'deltaRV: ' + str(c/newLambda2[len(newLambda2)/2] * deltaLambda) + ' m/s'
                
                #plotsssss
#                 plt.plot(newLambda1, flux1Neat/max(flux1Neat))
#                 plt.plot(newLambda2, flux2Neat/max(flux2Neat))
                plt.plot(newLambda2, a/max(a))
#                 plt.plot(lambda1, flux1/max(flux1))
#                 plt.plot(lambda2, flux2/max(flux2))
#                 plt.plot(lambda1, b/max(b))
                plt.title(fibre1['NAME'])
                plt.show()
            else:
                print fibre1[0].strip()
calculate_RV_shift_HERMES([1], ['24oct10025red.fits', '24oct10026red.fits', '24oct10027red.fits', '24oct10030red.fits', '24oct10031red.fits', '24oct10032red.fits'])


