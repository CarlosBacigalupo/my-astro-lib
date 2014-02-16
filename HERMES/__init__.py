import pyfits as pf
import scipy.constants as const
import scipy.optimize as opt
import numpy as np
import pylab as plt
import os
import shutil
import subprocess 
import math
from pymodelfit.builtins import VoigtModel
from pymodelfit.builtins import GaussianModel
from pymodelfit import get_model_instance
import time

import RVSimulator as RVS
import TableBrowser as TB
import toolbox

# calculate_RV_shift_HERMES([1], ['24oct10025red.fits', '24oct10026red.fits', '24oct10027red.fits', '24oct10030red.fits', '24oct10031red.fits', '24oct10032red.fits'])

class dr2df():
    
    date = ''
    base_dir = ''
    target_dir = ''
    dr_dir = ''
    file_ix = []
    

    def create_file_list(self):
        
        self.files1 =  [self.date +'1' + str(name).zfill(4)+ '.fits' for name in self.file_ix]
        self.files2 =  [self.date +'2' + str(name).zfill(4)+ '.fits' for name in self.file_ix]
        self.files3 =  [self.date +'3' + str(name).zfill(4)+ '.fits' for name in self.file_ix]
        self.files4 =  [self.date +'4' + str(name).zfill(4)+ '.fits' for name in self.file_ix]
        
        
    def create_folders(self):
        result = True 
        try:
            os.rmdir(self.target_dir)
        except OSError as ex:
            if ((ex.errno == 66) or (ex.errno == 17)):
                print "Target folder not empty."
                return result
            
        os.mkdir(self.target_dir)
        os.mkdir(self.target_dir+'1/')
        os.mkdir(self.target_dir+'2/')
        os.mkdir(self.target_dir+'3/')
        os.mkdir(self.target_dir+'4/')
                
        
        return result
    
            
    def auto_reduce(self):
        
        if self.create_folders():
            self.create_file_list()
        
            #copy files
            cam = 0
            for camList in [self.files1, self.files2, self.files3, self.files4]:
                cam += 1
                for i in camList:
                    src = self.base_dir  + 'ccd_' + str(cam) + '/' + i
                    dst = self.target_dir + '' + str(cam) + '/' + i
                    if not os.path.exists(dst):
                        shutil.copy(src, dst)
            
            # environment vars for os call
            env = {'PATH': self.dr_dir,
                   'DISPLAY': '/tmp/launch-Sqjh9Y/org.macosforge.xquartz:0',
                   }
            
                        
            #start reduction
            cam = 0  
            out = 0      
            for j in [self.files1, self.files2, self.files3, self.files4]:
                cam += 1
                os.chdir(self.target_dir + str(cam) + '/')
                
                #flat
                result = True 
                try:
                    os.rmdir(self.target_dir + str(cam) + '/' + j[self.flat][:-5]+'_outdir')
                except OSError as ex:
                    if ex.errno == 66:
                        print "Target folder not empty."
                        return False
                os.mkdir (j[self.flat][:-5]+'_outdir')                   
                os_command =  'drcontrol'
                os_command += ' reduce_fflat ' + j[self.flat]
                os_command += ' -idxfile hermes.idx'
                os_command += ' -OUT_DIRNAME '  + j[self.flat][:-5]+'_outdir'
                out = subprocess.call(os_command, env = env, shell = True)

                if out==0: #arc
                    result = True 
                     
                    try:
                        os.rmdir(self.target_dir + str(cam) + '/' + j[self.arc][:-5]+'_outdir')
                    except OSError as ex:
                        if ex.errno == 66:
                            print "Target folder not empty."
                            return False
                    os.mkdir(self.target_dir + str(cam) + '/' + j[self.arc][:-5]+'_outdir')                   
                    os_command =  'drcontrol'
                    os_command += ' reduce_arc ' + self.target_dir + str(cam) + '/' + j[self.arc]
                    os_command += ' -idxfile hermes.idx'
                    os_command += ' -TLMAP_FILENAME ' + self.target_dir + str(cam) + '/' + j[self.flat][:-5] + 'tlm.fits'
                    os_command += ' -OUT_DIRNAME ' + self.target_dir + str(cam) + '/' + j[self.arc][:-5]+'_outdir'
                    os_command =  'drcontrol'
                    os_command += ' reduce_arc '  + j[self.arc]
                    os_command += ' -idxfile hermes.idx'
                    os_command += ' -TLMAP_FILENAME ' + j[self.flat][:-5] + 'tlm.fits'
                    os_command += ' -OUT_DIRNAME ' + j[self.arc][:-5]+'_outdir'
                    out = subprocess.call(os_command, env = env, shell = True)
                     
                    if out==0: #scrunch flat
                        os_command =  'drcontrol'
                        os_command += ' reduce_fflat ' + j[self.flat]
                        os_command += ' -idxfile hermes.idx'
                        os_command += ' -WAVEL_FILENAME ' + j[self.arc][:-5] + 'red.fits'
                        os_command += ' -OUT_DIRNAME ' + j[self.flat][:-5]+'_outdir'
                        out = subprocess.call(os_command, env = env, shell = True)
                        
                        obj_files = j[:]
                        del obj_files[self.arc]
                        del obj_files[self.flat-1]   #-1 as i just removed 1 element....                       
                        for obj in obj_files:
                            result = True             
                            try:
                                os.rmdir(obj[:-5]+'_outdir')
                            except OSError as ex:
                                if ex.errno == 66:
                                    print "Target folder not empty."
                                    return False
                            os.mkdir(obj[:-5]+'_outdir')                   
                            if out==0: 
                                os_command =  'drcontrol'
                                os_command += ' reduce_object ' + obj
                                os_command += ' -idxfile hermes.idx'
                                os_command += ' -WAVEL_FILENAME ' + j[self.arc][:-5] + 'red.fits'
                                os_command += ' -TLMAP_FILENAME ' + j[self.flat][:-5] + 'tlm.fits'
                                os_command += ' -FFLAT_FILENAME ' + j[self.flat][:-5] + 'red.fits'
                                os_command += ' -OUT_DIRNAME ' + obj[:-5]+'_outdir'
#                                 os_command += ' -TPMETH OFFSKY'
                                os.system('killall drcontrol')
                                os.system('killall drexec')
                                out = subprocess.call(os_command, env = env, shell = True)
                                shutil.copyfile(obj[:-5]+'red.fits', self.red_dir + obj[:-5]+'red.fits')
                                
        
class RV():

    files = [] #science reduced files to analyse
    base_dir = ''    
    JS = 1.1574074074074073e-05  #julian second
    JD = []
    Lambda = []
    fileData = []
    fibreTable = []
    header = []
    c = const.c
    
    
    def plot_timeline(self):
        for i in range(len(self.files)):
            plt.plot(np.arange(self.JD[i], self.JD[i]+self.exposed[i],self.JS),np.ones(int(self.exposed[i]/self.JS)+1))
        plt.show()

    def read_fits_files(self):
        
        
        self.JD = np.zeros(len(self.files))
        self.exposed = np.zeros(len(self.files))

        for i in range(len(self.files)):    
            a = pf.open(self.base_dir + self.files[i])
            self.JD[i] = a[0].header['UTMJD']
            self.exposed[i] = a[0].header['EXPOSED']/24/3600
            
            temp = a[0].header
            self.header = self.header + [temp]
#             if (len(self.fileData)==0):
#                 self.header = np.array([temp])
#             else:
#                 self.header = np.append(self.header, [temp], 0)
            
            
            temp = RVS.extract_HERMES_wavelength(self.base_dir + self.files[i])
            if (len(self.fileData)==0):
                self.Lambda = np.array([temp])
            else:
                self.Lambda = np.append(self.Lambda, [temp], 0)
    
            temp = self.load_fibre_table(self.base_dir + self.files[i])
            if (len(self.fileData)==0):
                self.fibreTable = np.array([temp])
            else:
                self.fibreTable = np.append(self.fibreTable, [temp], 0)
    
            temp = RVS.extract_all_from_2dfdr(self.base_dir + self.files[i])
            if (len(self.fileData)==0):
                self.fileData = np.array([temp])
            else:
                self.fileData = np.append(self.fileData, [temp],0)
        
        print ' Read ' + str(len(self.files)) + ' files'
        
    def clear_by_index(self, index):
        
        if index<len(self.files):
            self.files.pop(index)
            self.header.pop(index)
            self.Lambda = np.delete(self.Lambda, index,0)
            self.fibreTable = np.delete(self.fibreTable, index,0)
            self.fileData = np.delete(self.fileData, index,0)

    def check_single_target_alignment_by_fibre(self, fibre):
        for i in range(len(self.fibreTable)-1):
             if np.all(self.fibreTable[i]['target'][fibre-1] == self.fibreTable[i+1]['target'][fibre-1])==False:
                 print 'DIFFERENT target in files ', self.files[i], self.files[i+1] 
                 print self.fibreTable[i]['target'][fibre-1], self.fibreTable[i+1]['target'][fibre-1]
                 print ''
             else:
                 print 'SAME target in files ', self.files[i], self.files[i+1] 
                 print self.fibreTable[i]['target'][fibre-1], self.fibreTable[i+1]['target'][fibre-1]
                 print ''
        return True


    def check_all_target_alignment(self):
        for i in range(len(self.fibreTable)-1):
             if np.all(self.fibreTable[i]['target'] == self.fibreTable[i+1]['target'])==False:
                 print 'DIFFERENT fibre table in files ', self.files[i], self.files[i+1] 
             else:
                 print 'SAME fibre table in files ', self.files[i], self.files[i+1] 
        return True

    
    def single_spectrum_by_target(self):
        print 'empty'
        
        
    def old_calculate_RV_shift(self):
        
        for j in range(399):     
            for i in range(fileData.shape[0]-1):
                fibre1 = fibreTable[i,j]
                fibre2 = fibreTable[i+1,j]
                # fibre types: F-guide, N-broken, dead or no fibre, P-program (science), S - Sky, U-unallocated or unused
                if ((fibre1['NAME'].strip()==fibre2['NAME'].strip()) & (fibre1['TYPE']=='P') & (fibre1['NAME'].strip().lower()!='parked')):
                    time1 = self.JD[i]
                    time2 = self.JD[i+1]
                    lambda1 = self.Lambda[i]
                    lambda2 = self.Lambda[i+1]
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
        
    
    def calculate_RV_shifts(self):
        
        self.RVShifts = np.zeros((len(self.files),400))
        
        for j in range(399):     
            fibre1 = self.fibreTable[0].ix[j]       
            # fibre types: F-guide, N-broken, dead or no fibre, P-program (science), S - Sky, U-unallocated or unused
            if ((fibre1['type']=='P') & (fibre1['target'].strip().lower()!='parked')):
                time1 = self.JD[0]
                lambda1 = self.Lambda[0]
                flux1 = self.fileData[0,j,:]
                newLambda1, flux1Neat = RVS.clean_flux(flux1, 10, lambda1 )
                
                for i in range(0,self.fileData.shape[0]):
    
                    fibre2 = self.fibreTable[i].ix[j]
                    if (fibre1['target'].strip()==fibre2['target'].strip()):
                        time2 = self.JD[i]
                        lambda2 = self.Lambda[i]
                        flux2 = self.fileData[i,j,:]
                        newLambda2, flux2Neat = RVS.clean_flux(flux2, 10, lambda2)                   
                        
                        a =  np.correlate(flux1Neat, flux2Neat, 'same')
#                         b =  np.correlate(flux1, flux2, 'same')
#                         print ''
#                         print fibre1['target']
#                         print 'Time Difference: ' + str(abs(time1 - time2)*24*60) + ' min'
                        deltaLambda = newLambda2[len(newLambda2)/2] - newLambda2[np.where(a==max(a))[0][0]]
#                         print 'Wavelength difference: ' + str(deltaLambda) + ' Ang'
#                         print 'deltaRV: ' + str(self.c/newLambda2[len(newLambda2)/2] * deltaLambda) + ' m/s'
                        RV = self.c/newLambda2[len(newLambda2)/2] * deltaLambda
                        self.RVShifts[i,j] = RV
                        print i,j, RV
                        #plotsssss
        #                 plt.plot(newLambda1, flux1Neat/max(flux1Neat))
        #                 plt.plot(newLambda2, flux2Neat/max(flux2Neat))
#                         plt.plot(newLambda2, a/max(a))
#         #                 plt.plot(lambda1, flux1/max(flux1))
#         #                 plt.plot(lambda2, flux2/max(flux2))
#         #                 plt.plot(lambda1, b/max(b))
#                         title = fibre1['target']+ ' ('+fibre1['mag']+')'
#                         plt.title(title)
#                         plt.show()
            else:
                pass
#                 print 'Skipping ', fibre1['target'].strip()
    

    def load_fibre_table(self, fileName):
        
        self.source = fileName.split('.')[-1]
        
        if self.source == 'fld':
            a = np.loadtxt(fileName, skiprows = 9, dtype=str)
            b = a.transpose()

            self.pivot = np.zeros(len(b[0].astype(str)))
            self.target = b[0].astype(str)
            self.RA_h = b[1].astype(int)
            self.RA_min = b[2].astype(int)
            self.RA_sec = b[3].astype(float)
            self.Dec_deg = b[4].astype(int)
            self.Dec_min = b[5].astype(int)
            self.Dec_sec = b[6].astype(float)
            self.type = b[7]
            self.mag = b[9].astype(float)

            self.RA_dec = 15*(self.RA_h + self.RA_min /60. + self.RA_sec /3600.)     
            self.Dec_dec = self.Dec_deg + self.Dec_min /60. + self.Dec_sec /3600.
                        
        elif self.source == 'lis':
            a = np.loadtxt(fileName, skiprows = 9, dtype=str)
            b = a.transpose()[1:]

            self.pivot = b[0].astype('int')
            self.target = b[1].astype('str')
            self.RA_h = b[2].astype('int')
            self.RA_min = b[3].astype('int')
            self.RA_sec = b[4].astype('float')
            self.Dec_deg = b[5].astype('int')
            self.Dec_min = b[6].astype('int')
            self.Dec_sec = b[7].astype('float')
            self.type = b[8]
            self.mag = b[10].astype('float')
        
            self.mag[self.type=='1'] = 0.
            self.type[self.type=='1'] = 0.

            self.RA_dec = 15*(self.RA_h + self.RA_min /60. + self.RA_sec /3600.)     
            self.Dec_dec = self.Dec_deg + self.Dec_min /60. + self.Dec_sec /3600.
            
        elif ((self.source.lower() == 'fits') or (self.source.lower() == 'fit')):
            a = pf.open(fileName)    
            for i in range(1,len(a)):
                fibreType = a[i].header['EXTNAME']
                if fibreType=='FIBRES' : b = a[i].data
            
            self.fibre = np.arange(1,401)
            self.pivot = b.field('PIVOT')
            self.target = b.field('NAME').strip()
            self.RA_dec = b.field('RA')
            self.Dec_dec = b.field('DEC')
            self.mag = b.field('MAGNITUDE')
            self.type = b.field('TYPE')

            self.RA_h, self.RA_min, self.RA_sec = toolbox.dec2sex(self.RA_dec/15)   
            self.Dec_deg, self.Dec_min, self.Dec_sec = toolbox.dec2sex(self.Dec_dec)
            
        return self.create_dataframe()
    
    def create_dataframe(self):
        
        d = [ 'pivot','target', 'RA_h', 'RA_min', 'RA_sec', 'RA_dec' , 'Dec_deg', 'Dec_min',  'Dec_sec', 'Dec_dec', 'type', 'mag']
       
        tableFrame = np.array([self.pivot, 
                               self.target, 
                               self.RA_h, 
                               self.RA_min, 
                               self.RA_sec, 
                               self.RA_dec, 
                               self.Dec_deg, 
                               self.Dec_min, 
                               self.Dec_sec, 
                               self.Dec_dec, 
                               self.type, 
                               self.mag])
        
        return TB.build_DataFrame(tableFrame, d)


class PSF():
    
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    sexParamFile = base_dir + 'HERMES.sex'
    outputFileName = base_dir + 'out.txt'
    scienceFile = base_dir + '10nov40045.fits'
    biasFile = base_dir + 'BIAScombined4.fits'
    flatFile = base_dir + '10nov10044.fits'
    tramLinefile = base_dir + '10nov10044tlm.fits'
    nFibres = 10

    def __init__(self):
        self.nFibres = 400
        self.nBundles = 40
        self.pShort = []
        
    def bias_subtract_from_overscan(self, image, range):
        
        biasArrayInit = image[:,range]
        
        if len(range) > 1:          
            biasLine = np.median(biasArrayInit,1)
        else:
            biasLine = biasArrayInit
            
        biasArray = np.tile(biasLine,(self.imWidth,1)).transpose()
            
        return image - biasArray      
        
    def bias_subtract(self):
        self.scienceIm_b = self.scienceIm - self.biasIm       
        self.flatIm_b = self.flatIm - self.biasIm        
            
    def open_files(self):
        
         #Read, calibrate and create im objects
        self.bias = pf.open(self.biasFile)
        self.tlm = pf.open(self.tramLinefile)
        self.flat = pf.open(self.flatFile)
        self.science = pf.open(self.scienceFile)

        self.imWidth = self.science[0].header['NAXIS1']
        self.imHeight = self.science[0].header['NAXIS2']
        self.scienceIm = self.science[0].data     
        self.tlmIm = self.tlm[0].data
        self.flatIm = self.flat[0].data
        self.biasIm = self.bias[0].data

    def fit_10f(self, flatCurve):
        #findes best fit for 10f config (5 x 2f) for len(p)=40
        
        self.write_p10(flatCurve)
        
        factorTry = 100
        diagTry = np.ones(len(self.p10)) * 10
        diagTry = np.arange(1,len(self.p10)+1)
        maxfev = int(1e4) # of max tries
        # diagonals!=1 to follow:
#         diagTry[0] = 10    
#         diagTry[self.nFibres] = 0.1  #sigma
#         diagTry[self.nFibres+1] = 100  #gamma
#         diagTry[-self.nFibres:] = np.ones(self.nFibres) * 10
        
        self.reps=0
        p = opt.leastsq(self.find_residuals, self.p10, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)
        
        finalCurve, onePSF = self.find_p10_curve(p[0],flatCurve, output = False)
        
#         plt.plot(flatCurve, label = 'image')
#         plt.plot(finalCurve,label = 'model')
#         plt.title(self.profile)
#         plt.legend()
#         plt.show()

        return p[0][10],p[0][11], p[2]['fvec'], np.sum(p[2]['fvec'])/np.sum(flatCurve)
    
    def fit_pShort(self, flatCurve):
       
        self.write_pShort(flatCurve)
  
        factorTry = 1
        diagTry = np.ones(len(self.pShort))
#         diagTry[0] = 200
#         diagTry[1] = 200
#         diagTry[2] = 0.1
#         diagTry[3] = 1
#         diagTry[4] = 30
#         diagTry[5] = 1
#         diagTry[6] = 1
#         diagTry[7] = 1
#         diagTry[8] = 0.1
#         diagTry[9] = 0.1
        maxfev = int(1e3)
        self.reps=0
        p = opt.leastsq(self.find_residuals, self.pShort, full_output=True, args = flatCurve)#, factor = factorTry, diag = diagTry, maxfev = maxfev)



#         plt.plot(self.find_curve(p[0]))
        finalCurve, onePSF = self.find_pShort_curve(p[0],flatCurve, output = False)
        plt.plot(flatCurve, label = 'flat')
        plt.plot(finalCurve,label = 'model')
#         plt.plot(onePSF/max(onePSF),label = 'Single PSF')
        plt.legend()
        plt.show()
        return p
    
    def fit_p400_A_mu(self, flatCurve):
  
        factorTry = 100
        diagTry = np.ones(len(self.p400_A_mu))
        maxfev = int(1e3)
        self.reps=0

        p = opt.leastsq(self.find_residuals, self.p400_A_mu, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)
        print p[0]
        finalCurve = self.find_p400_A_mu_curve(p[0],flatCurve, output = False)
        
#         #comarisson with pymodelfit
#         self.reps=0
#         self.profile = 'gaussian_fit'
#         p = opt.leastsq(self.find_residuals, self.p400_A_mu, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)
#         finalCurve = self.find_p400_A_mu_curve(p[0],flatCurve, output = False)

        
        self.plot(flatCurve, finalCurve)


        return p
    
    def fit_p400_sigma_gamma(self, flatCurve):
  
        factorTry = 1
        diagTry = np.ones(len(self.p400_sigma_gamma))
#         diagTry[0] = 200
#         diagTry[1] = 200
#         diagTry[2] = 0.1
#         diagTry[3] = 1
#         diagTry[4] = 30
#         diagTry[5] = 1
#         diagTry[6] = 1
#         diagTry[7] = 1
#         diagTry[8] = 0.1
#         diagTry[9] = 0.1
        maxfev = int(1e1)
        self.reps=0
        p = opt.leastsq(self.find_residuals, self.p400_sigma_gamma, full_output=True, args = flatCurve)#, factor = factorTry, diag = diagTry, maxfev = maxfev)

# #         plt.plot(self.find_curve(p[0]))
#         finalCurve, onePSF = self.find_p400_sigma_gamma_curve(p[0],flatCurve, output = False)
#         plt.plot(flatCurve, label = 'flat')
#         plt.plot(finalCurve,label = 'model')
# #         plt.plot(onePSF/max(onePSF),label = 'Single PSF')
#         plt.legend()
#         plt.show()
        return p
   
    def write_p10(self, flatCurve):
        
        firstGap = 24.5
        lGap = 17.
        sGap = 8.0001   

        sigma = 1.9
        gamma = 2.24741646e-01
        
#         sigma = np.ones(10) * sigma
#         gamma = np.ones(10) * gamma
        
        A = np.zeros(10)
#         A[0] = 4.30763396
#         A[1] = 4.79011575
#         A[2] = 3.17765607
#         A[3] = 3.19855864
#         A[4] = 3.91087569
#         A[5] = 3.83988602
#         A[6] = 3.44450670
#         A[7] = 4.28972608
#         A[8] = 3.39712514
#         A[9] = 3.48422336
        A = np.ones(self.nFibres) * np.max(flatCurve) * np.sqrt(2* math.pi * sigma**2) *2

        mu = np.zeros(10)
        mu[0] = 660
        mu[1] = mu[0] + sGap
        mu[2] = 1888
        mu[3] = mu[2] + sGap
        mu[4] = 2628
        mu[5] = mu[4] + sGap
        mu[6] = 3043
        mu[7] = mu[6] + sGap
        mu[8] = 4066
        mu[9] = mu[8] + sGap
        
        self.p10 = np.hstack((A, sigma, gamma, mu))
        
    def write_pShort(self, flatCurve):

        sigma = 1.6775996937959061
        gamma = 5.9997168955962937e-06        
        A = np.ones(self.nFibres) * np.max(flatCurve) * np.sqrt(2* math.pi * sigma**2)
        mu = np.zeros(400)
        self.pShort = np.hstack((A, sigma, gamma, mu))        
        
    def write_p400_A_mu(self, flatCurve):  
             
        sigma = np.tile(1.,400)
        gamma = 0.2
#         A = np.sqrt(2*np.pi*sigma**2)*flatCurve[self.mu.astype(int)].astype(float)
        A = flatCurve[self.mu.astype(int)].astype(float)
        self.p400_A_mu = np.hstack((A, sigma, gamma))        

    def write_p400_sigma_gamma(self, flatCurve):
         
        sigma = 1.9
        gamma = 0.2
        sigma = np.ones(400) * sigma
        gamma = np.ones(400) * gamma
        self.p400_sigma_gamma = np.hstack((sigma, gamma))        
    
    def write_mu_from_tlm(self, column):
        
        self.mu = self.tlmIm[:,column]-0.5
        
    def find_residuals(self, p, flatCurve, output=True):
        self.reps += 1

        if len(p)==len(self.pShort):
            model = self.find_pShort_curve(p, flatCurve)[0]
            
        elif len(p)==len(self.p400_A_mu):
            model = self.find_p400_A_mu_curve(p, flatCurve)

        elif len(p)==len(self.p400_sigma_gamma):
            model = self.find_p400_sigma_gamma_curve(p, flatCurve)[0]

        elif len(p)==len(self.p10):
            model = self.find_p10_curve(p, flatCurve, output=output)[0]

        diff = flatCurve-model
#         plt.plot(flatCurve)
#         plt.plot(model)
#         plt.show() 
        
        #just output to see progress
        if (((self.reps in [1,2,3,4,5,6,7,8,10,20,50,100, 500, 1000]) or (self.reps in range(1100,10000,100))) and (output==True)):
            print self.reps,
        if ((self.reps==1) and (output==True)):
            self.plot(flatCurve, model)
            
        return diff
 
    def find_p10_curve(self, p, flatCurve, output = False, psfFibre = []):
         
        A = p[0:self.nFibres]
        #Sigma, gamma for 10 params
        #sigma = p[self.nFibres:2*self.nFibres]
        #gamma = p[2*self.nFibres:3*self.nFibres]
        #Sigma, gamma for 1 param
        sigma = p[self.nFibres]
        gamma = p[self.nFibres+1]
        mu = p[-self.nFibres:]
                
        if self.profile =='voigt':
            a = get_model_instance('Voigt')
        elif self.profile =='gaussian':
            a = get_model_instance('Gaussian')
            
        model = np.zeros(len(flatCurve))
        onePSF = np.zeros(len(flatCurve))
        
        for i in range(self.nFibres):
            thisA = A[i]
            thisSigma = sigma
            thisGamma = math.fabs(gamma)
#             thisGamma = gamma
            thisMu = mu[i]
    
            if self.profile =='voigt':
                thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisGamma,thisMu)
            elif self.profile =='gaussian':
                thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisMu)
            
            if ((np.sum(thisCurve) >0) and (output == True)):
                print 'Input parameters for fibre ' + str(i) + ':', thisA, thisSigma, thisGamma, thisMu
            
            model += thisCurve
#         plt.plot(flatCurve)
#         plt.plot(model)
#         plt.show()
        
        return model, onePSF
    
    def find_pShort_curve(self, p, flatCurve, output = False):

        A = p[0:self.nFibres]
        sigma = p[self.nFibres]
        gamma = p[self.nFibres+1]
        mu = p[-self.nFibres:]
                
        if self.profile =='voigt':
            a = get_model_instance('Voigt')
        elif self.profile =='gaussian':
            a = get_model_instance('Gaussian')
            
        model = np.zeros(len(flatCurve))
        onePSF = np.zeros(len(flatCurve))
        
        for i in range(self.nFibres):
            thisA = A[i]
            thisSigma = sigma
#             thisGamma = math.fabs(gamma)
            thisGamma = gamma
            thisMu = mu[i]
            if i+1 in self.deadFibres:
                thisCurve = np.zeros(len(flatCurve))
            else:
                if self.profile =='voigt':
                    thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisGamma,thisMu)
                elif self.profile =='gaussian':
                    thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisMu)
            
            if ((np.sum(thisCurve) >0) and (output == True)):
                print 'Input parameters for fibre ' + str(i) + ':', thisA, thisSigma, thisGamma, thisMu
            
            model += thisCurve
            
#         plt.plot(flatCurve)
#         plt.plot(model)
#         plt.show()                    
        return model, onePSF
    
    def find_p400_A_mu_curve(self, p, flatCurve, output = False):

        A = np.array(p[0:self.nFibres])
        sigma = np.array(p[self.nFibres:2*self.nFibres])
        gamma = np.array(p[-1])

#         if self.profile =='voigt':
#             a = get_model_instance('Voigt')
#         elif self.profile =='gaussian':
#             a = get_model_instance('Gaussian')
#             
#         
# 
#         for i in range(self.nFibres):
#             thisA = A[i]
#             thisSigma = sigma
#             thisGamma = gamma
#             thisMu = self.mu[i]
#             
#             if i+1 in self.deadFibres:
#                 thisCurve = np.zeros(len(flatCurve))
#             else:
#                 if self.profile =='voigt':
#                     thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisGamma,thisMu)
#                 elif self.profile =='gaussian':
#                     thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisMu)
#             
#             if ((np.sum(thisCurve) >0) and (output == True)):
#                 print 'Input parameters for fibre ' + str(i) + ':', thisA, thisSigma, thisGamma, thisMu
#             
#             model += thisCurve

#         start_time=time.time()
        thisA = A
        thisSigma = sigma
        thisGamma = np.array([gamma])
        thisMu = self.mu
        
        #set amplitude of deadfibres to 0
        thisA[np.array(self.deadFibres[self.cur_camera])-1]=0
        
        if self.profile =='voigt':
            thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisGamma,thisMu)
        elif self.profile =='gaussian':
            model = toolbox.gaussian(np.arange(len(flatCurve)), sigma = thisSigma, mu = thisMu, A = thisA)                    
        elif self.profile =='gaussian_fit':
            model = np.zeros(len(flatCurve))
            a = get_model_instance('Gaussian')
            for i in range(self.nFibres):
                thisA = A[i]
                thisSigma = sigma[i]
                thisGamma = gamma
                thisMu = self.mu[i]
                
                if i+1 in self.deadFibres:
                    thisCurve = np.zeros(len(flatCurve))
                else:
                    thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisMu)
                
            
                model += thisCurve
                
        return model
    
    def find_p400_sigma_gamma_curve(self, p, flatCurve, output = False):

        sigma = p[:self.nFibres]
        gamma = p[self.nFibres:]
                
        if self.profile =='voigt':
            a = get_model_instance('Voigt')
        elif self.profile =='gaussian':
            a = get_model_instance('Gaussian')
            
        model = np.zeros(len(flatCurve))
        onePSF = np.zeros(len(flatCurve))
        
        for i in range(self.nFibres):
            thisA = self.A[i]
            thisSigma = sigma[i]
#             thisGamma = math.fabs(gamma[i])
            thisGamma = gamma[i]
            thisMu = self.mu[i]
            if i+1 in self.deadFibres:
                thisCurve = np.zeros(len(flatCurve))
            else:
                if self.profile =='voigt':
                    thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisGamma,thisMu)
                elif self.profile =='gaussian':
                    thisCurve = a.f(np.arange(len(flatCurve)), thisA,thisSigma,thisMu)
            
            if ((np.sum(thisCurve) >0) and (output == True)):
                print 'Input parameters for fibre ' + str(i) + ':', thisA, thisSigma, thisGamma, thisMu
            
            model += thisCurve 
            
#         plt.plot(flatCurve)
#         plt.plot(model)
#         plt.show()                    
        return model, onePSF        
    
    def plot(self, flatCurve = [], model=[], extraCurve =[]):
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        if flatCurve !=[]:ax1.plot(flatCurve, label = 'flat')
        if model !=[]:ax1.plot(model,label = 'model')
        if extraCurve !=[]:ax1.plot(extraCurve,label = 'extra')
        
        ax2 = ax1.twiny()
        ax2.set_xlabel("Fibre")
        ax2.set_xlim(ax1.axis()[0:2])
        ax2.set_xticks(self.mu)
        ax2.set_xticklabels(range(1,401))
        ax1.legend()
        plt.show()    
        
    def read_full_image_spatial(self, profile, Columns):

        self.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
        self.scienceFile = self.base_dir + '10nov' + str(self.cur_camera) + '0044.fits'
        self.biasFile = self.base_dir + 'BIAScombined2.fits'
        self.profile = profile
        self.out_dir = self.base_dir + ''
        self.flatFile = self.base_dir + '10nov' + str(self.cur_camera) + '0044.fits'
        self.tramLinefile = self.base_dir + '10nov' + str(self.cur_camera) + '0044tlm.fits'
        self.deadFibres = [[],[],[],[], []]
        self.deadFibres[1] = [41,66, 91, 141,191, 241, 250, 294, 304, 341, 391]
        self.deadFibres[2] = [41,66, 91, 141,191,241, 294, 307, 341, 391]
        self.deadFibres[3] = [41,66, 91, 141,191,241, 250, 294, 307, 341, 391]
        self.deadFibres[4] = [41,66, 91, 141,191,241, 250, 294, 307, 341, 391]
        self.nFibres = 400
        self.open_files()
        self.flatIm_b = self.bias_subtract_from_overscan(self.flatIm, range(-45,-5)) #bias subtract using the last 40 columns(overscan region)

        if profile=='voigt':
            fileOutSpa = open(self.out_dir + 'spatialV'+str(self.cur_camera)+'.txt','w')
        elif profile=='gaussian':
            fileOutSpa = open(self.out_dir + 'spatialG'+str(self.cur_camera)+'.txt','w')
        fileOutSpa.write('fibre, y_cent, A, sigma, gamma, mu \n')
        
        for C in Columns:
            flatCurve = self.flatIm_b[:,C]       
            
            start_time=time.time()
            #find A        
            self.write_mu_from_tlm(C)
            self.write_p400_A_mu(flatCurve)
            p = self.fit_p400_A_mu(flatCurve)
            print p
            print C, 'column fit took ', time.time() - start_time, "seconds" 
            

            A = p[0][:self.nFibres]           
            sigma = p[0][self.nFibres:2*self.nFibres]
            gamma = p[0][-1]
            
            for F in range(400):
                outString = (str(F+1) + ', ' + 
                            str(C) + ', ' + 
                            str(A[F]) + ', ' + 
                            str(sigma[F]) + ', ' + 
                            str(gamma) + ', ' + 
                            str(self.mu[F]) + ', ' + 
                            '\n' )
                fileOutSpa.write(outString)
            
    def read_psf_data(self, profile, referenceFile, camera):
        import matplotlib
        import matplotlib.cm as cm
        import matplotlib.mlab as mlab
        import matplotlib.pyplot as plt
        
    #     from scipy.interpolate import interp1d
    #     from scipy.interpolate import interp2d
    #     import pyfits as pf
        if camera==1:
            cameraName = 'Blue Camera'
        elif camera==2:
            cameraName = 'Green Camera'
        elif camera==3:
            cameraName = 'Red Camera'
        elif camera==4:
            cameraName = 'IR Camera'
            
        
        self.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
        self.out_dir = self.base_dir + 'output/'
        if profile=='voigt':
            dataFile = open(self.out_dir + 'spatialV'+str(camera)+'.txt','r')
        elif profile=='gaussian':
            dataFile = open(self.out_dir + 'spatialG'+str(camera)+'.txt','r')

        a = np.loadtxt(dataFile, skiprows=1, delimiter = ',' , usecols=[0,1,2,3,4,5]).transpose()
    
        Lambda = RVS.extract_HERMES_wavelength(self.base_dir + referenceFile)

        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        
        X, Y = np.meshgrid(np.unique(a[1]), np.unique(a[0]))
        Z = a[3].reshape(X.shape[1],X.shape[0])
        Z = Z.transpose()

        Z[Z<0.5] = np.average(Z)
        Z[Z>4] = np.average(Z)
        
        if profile=='voigt':
            Z = 2 * np.sqrt(2*math.log(2)) * Z
        elif profile=='gaussian':
            Z = 2 * np.sqrt(2*math.log(2)) * Z

        # Or you can use a colormap to specify the colors; the default
        # colormap will be used for the contour lines
        plt.figure()
        im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                        cmap=cm.gray, extent=(0,400,1,400))
#         plt.xticks(range(40), Lambda.astype(int)[::len(Lambda)/40], size='small')
        plt.xticks(range(0,400,100),Lambda.astype(int)[::len(Lambda)/100])
        levels = np.arange(2.5, 7., 0.5)
        CS = plt.contour(Z, levels,
                         origin='lower',
                         cmap=cm.gray, 
                         linewidths=1,
                         extent=(0,400,1,400))

        #Thicken the zero contour.
#         zc = CS.collections[6]
#         plt.setp(zc, linewidth=4)
         
        plt.clabel(CS, levels,  # levels[1::2]  to label every second level
                   inline=0,
                   fmt='%1.2f',
                   fontsize=12)
        
        # make a colorbar for the contour lines
#         CB = plt.colorbar(CS, shrink=0.8, extend='both')
        
        plt.title('Spatial FWHM - ' + cameraName)
        plt.gray()  # Now change the colormap for the contour lines and colorbar
#         plt.flag()
        
        # We can still add a colorbar for the image, too.
        CBI = plt.colorbar(im, orientation='vertical', shrink=1)
        
        # This makes the original colorbar look a bit out of place,
        # so let's improve its position.
        
#         l,b,w,h = plt.gca().get_position().bounds
#         ll,bb,ww,hh = CB.ax.get_position().bounds
#         CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
        plt.xlabel('Wavelength [Ang]')
        plt.ylabel('Fibre')
        plt.savefig(self.out_dir + 'contour_' + str(camera) + 'G.png')
        plt.show()

        
        
        #X-binned  plot
        fibreFlux = np.zeros(40)
        fibreFluxStd = np.zeros(40)
        x = np.zeros(40)
        for fib in [0,99,199,299,399]:
            for i in range(0,37,4):
                fibreFlux[i] = np.average(Z[fib,i:i+10])
                fibreFluxStd[i] = Z[fib,i:i+10].std()
            mask = [fibreFlux!=0]
            plt.errorbar(Lambda[mask], fibreFlux[mask], yerr=fibreFluxStd[mask], label = 'Fibre ' + str(fib+1))
        plt.title('FWHM - ' + cameraName)
        plt.xlabel('Wavelength [Ang]')
        plt.ylabel('PSF [px]')
        plt.legend()
        plt.savefig(self.out_dir + 'lambda_' + str(camera) + 'G.png')
        plt.show()
        
        #y-binned  plot
        fibreFlux = np.zeros(400)
        fibreFluxStd = np.zeros(400)
        x = np.zeros(400)
        for col in [0,100,200,300,400]:
            for i in range(0,359,40):
                fibreFlux[i] = np.average(Z[i:i+40,col])
                fibreFluxStd[i] = Z[i:i+40,col].std()
            mask = [fibreFlux!=0]
            plt.errorbar(np.array(range(len(fibreFlux)))[mask], fibreFlux[mask], yerr=fibreFluxStd[mask], label = 'column ' + str(col*10))
        plt.title('FWHM - '+cameraName)
        plt.xlabel('Fibre')
        plt.ylabel('PSF [px]')
        plt.legend()
        plt.savefig(self.out_dir + 'fibre_' + str(camera) + 'G.png')
        plt.show()
        
        
        
#         rows = 4112
#         cols = 4146
#         bins = 4
#         #bin rows
#         for i in range(bins):
#             binLimits = (rows/bins*i,rows/bins*(i+1))
#             mapBin = ((a[1]>binLimits[0]) & (a[1]<=binLimits[1]))
#             b = a[0,mapBin], a[5,mapBin]
#             c = np.bincount(b[0].astype(int), weights=b[1])
#             d = np.bincount(b[0].astype(int))
#             px = np.arange(len(d))
#             px = px[d>0]
#             c = c[d>0]
#             d = d[d>0]
#             c = c/d
#             plt.plot(px,c, label= ' Range (y) = ' + str(binLimits))
#             plt.title('Spatial PSF - y-binned - ' + method)
#             plt.xlabel('X[Px]')
#             plt.ylabel('Std. Dev. (px)')
#             plt.legend()
# #         plt.savefig(base_dir + 'spa_' + str(cam) + '_y_' + method +'.png')
#         plt.show()
#     
#         #bin cols
#         for i in range(bins):
#             binLimits = (cols/bins*i,cols/bins*(i+1))
#             mapBin = ((a[0]>binLimits[0]) & (a[0]<=binLimits[1]))
#             b = a[1,mapBin], a[5,mapBin]
#             c = np.bincount(b[0].astype(int), weights=b[1])
#             d = np.bincount(b[0].astype(int))
#             px = np.arange(len(d))
#             px = px[d>0]
#             c = c[d>0]
#             d = d[d>0]
#             c = c/d
#             plt.plot(px,c, label= ' Range (x) = ' + str(binLimits))
#             plt.title('Spatial PSF - x-binned - ' + method)
#             plt.xlabel('Y-Pixels')
#             plt.ylabel('Std. Dev. (px)')
#             plt.legend()
#         plt.savefig(base_dir + 'spa_' + str(cam) + '_x_' + method +'.png')
#         plt.show()
    
    
    #     grid = np.zeros((4112,4146))
    #     a[2] = a[2]/max(a[2])*65535
    #     pointSize = 10
    #     for i in range(len(a[0])):
    #         grid[int(a[1][i]-pointSize):int(a[1][i]+pointSize),int(a[0][i]-pointSize):int(a[0][i]+pointSize)] = a[2][i]
    #     plt.imshow(grid, origin='lower')
    #     plt.set_cmap(cm.Greys_r)                                                                                    
    #     plt.show()
    #     pf.writeto(base_dir + 'a.fits', grid)
    
    
        #interpolate function
    #     f = interp2d(a[0][range(0,len(a[0]),1)],a[1][range(0,len(a[0]),1)],a[5][range(0,len(a[0]),1)])#, kind='cubic')
    # #     grid = np.zeros((4112,4146))
    # #     grid = f(a[0].astype('int'), a[1].astype('int'))
    #     grid = f(range(4112),range(4146))
    # #     grid = f(range(100),range(100))
    #     plt.imshow(grid, origin='lower')
    #     plt.set_cmap(cm.Greys_r)                                                                                    
    #     plt.show()
    #     pf.writeto(base_dir + 'a.fits',grid, clobber= True)

    