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
        self.p400_A_mu = []
        
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
  
        factorTry = 1
        diagTry = np.ones(len(self.p400_A_mu))
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
        p = opt.leastsq(self.find_residuals, self.p400_A_mu, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)



#         plt.plot(self.find_curve(p[0]))
#         finalCurve, onePSF = self.find_p400_A_mu_curve(p[0],flatCurve, output = False)
#         plt.plot(flatCurve, label = 'flat')
#         plt.plot(finalCurve,label = 'model')
# #         plt.plot(onePSF/max(onePSF),label = 'Single PSF')
#         plt.legend()
#         plt.show()
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
        p = opt.leastsq(self.find_residuals, self.p400_sigma_gamma, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)

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
        
#         firstGap = 24.5
#         lGap = 17.
#         sGap = 8.0001   

        sigma = 1.6775996937959061
        gamma = 5.9997168955962937e-06
        
#         sigma = np.ones(10) * sigma
#         gamma = np.ones(10) * gamma
        
        A = np.ones(self.nFibres) * np.max(flatCurve) * np.sqrt(2* math.pi * sigma**2)

        mu = np.zeros(400)
        mu[0]=24.5
        mu[1]=32.501
        mu[2]=40.502
        mu[3]=48.503
        mu[4]=56.504
        mu[5]=64.505
        mu[6]=72.506
        mu[7]=80.507
        mu[8]=88.508
        mu[9]=96.509
        mu[10]=120
        mu[11]=128.1
        mu[12]=136.2
        mu[13]=144.3
        mu[14]=152.4
        mu[15]=160.5
        mu[16]=168.6
        mu[17]=176.7
        mu[18]=184.8
        mu[19]=192.9
        mu[20]=214
        mu[21]=222.1
        mu[22]=230.2
        mu[23]=238.3
        mu[24]=246.4
        mu[25]=254.5
        mu[26]=262.6
        mu[27]=270.7
        mu[28]=278.8
        mu[29]=286.9
        mu[30]=312
        mu[31]=320.05
        mu[32]=328.1
        mu[33]=336.15
        mu[34]=344.2
        mu[35]=352.25
        mu[36]=360.3
        mu[37]=368.35
        mu[38]=376.4
        mu[39]=384.45
        mu[40]=409
        mu[41]=417.4
        mu[42]=425.8
        mu[43]=434.2
        mu[44]=442.6
        mu[45]=451
        mu[46]=459.4
        mu[47]=467.8
        mu[48]=476.2
        mu[49]=484.6
        mu[50]=508
        mu[51]=516.4
        mu[52]=524.8
        mu[53]=533.2
        mu[54]=541.6
        mu[55]=550
        mu[56]=558.4
        mu[57]=566.8
        mu[58]=575.2
        mu[59]=583.6
        mu[60]=608
        mu[61]=616.4
        mu[62]=624.8
        mu[63]=633.2
        mu[64]=641.6
        mu[65]=650
        mu[66]=658.4
        mu[67]=666.8
        mu[68]=675.2
        mu[69]=683.6
        mu[70]=707
        mu[71]=715.4
        mu[72]=723.8
        mu[73]=732.2
        mu[74]=740.6
        mu[75]=749
        mu[76]=757.4
        mu[77]=765.8
        mu[78]=774.2
        mu[79]=782.6
        mu[80]=810
        mu[81]=818.4
        mu[82]=826.8
        mu[83]=835.2
        mu[84]=843.6
        mu[85]=852
        mu[86]=860.4
        mu[87]=868.8
        mu[88]=877.2
        mu[89]=885.6
        mu[90]=910
        mu[91]=918.4
        mu[92]=927.3
        mu[93]=936.2
        mu[94]=945.1
        mu[95]=954
        mu[96]=962.9
        mu[97]=971.8
        mu[98]=980.7
        mu[99]=989.6
        mu[100]=1013
        mu[101]=1021.9
        mu[102]=1030.8
        mu[103]=1039.7
        mu[104]=1048.6
        mu[105]=1057.5
        mu[106]=1066.4
        mu[107]=1075.3
        mu[108]=1084.2
        mu[109]=1093.1
        mu[110]=1117
        mu[111]=1125.9
        mu[112]=1134.8
        mu[113]=1143.7
        mu[114]=1152.6
        mu[115]=1161.5
        mu[116]=1170.4
        mu[117]=1179.3
        mu[118]=1188.2
        mu[119]=1197.1
        mu[120]=1219
        mu[121]=1227.9
        mu[122]=1236.8
        mu[123]=1245.7
        mu[124]=1254.6
        mu[125]=1263.5
        mu[126]=1272.4
        mu[127]=1281.3
        mu[128]=1290.2
        mu[129]=1299.1
        mu[130]=1324
        mu[131]=1332.9
        mu[132]=1341.8
        mu[133]=1350.7
        mu[134]=1359.6
        mu[135]=1368.5
        mu[136]=1377.4
        mu[137]=1386.3
        mu[138]=1395.2
        mu[139]=1404.1
        mu[140]=1428
        mu[141]=1436.9
        mu[142]=1445.8
        mu[143]=1454.7
        mu[144]=1463.6
        mu[145]=1472.5
        mu[146]=1481.4
        mu[147]=1490.3
        mu[148]=1499.2
        mu[149]=1508.1
        mu[150]=1533
        mu[151]=1541.9
        mu[152]=1550.8
        mu[153]=1559.7
        mu[154]=1568.6
        mu[155]=1577.5
        mu[156]=1586.4
        mu[157]=1595.3
        mu[158]=1604.2
        mu[159]=1613.1
        mu[160]=1637
        mu[161]=1645.9
        mu[162]=1654.8
        mu[163]=1663.7
        mu[164]=1672.6
        mu[165]=1681.5
        mu[166]=1690.4
        mu[167]=1699.3
        mu[168]=1708.2
        mu[169]=1717.1
        mu[170]=1744
        mu[171]=1752.9
        mu[172]=1761.8
        mu[173]=1770.7
        mu[174]=1779.6
        mu[175]=1788.5
        mu[176]=1797.4
        mu[177]=1806.3
        mu[178]=1815.2
        mu[179]=1824.1
        mu[180]=1849
        mu[181]=1857.9
        mu[182]=1866.8
        mu[183]=1875.7
        mu[184]=1884.6
        mu[185]=1893.5
        mu[186]=1902.4
        mu[187]=1911.3
        mu[188]=1920.2
        mu[189]=1929.1
        mu[190]=1954
        mu[191]=1962.9
        mu[192]=1971.8
        mu[193]=1980.7
        mu[194]=1989.6
        mu[195]=1998.5
        mu[196]=2007.4
        mu[197]=2016.3
        mu[198]=2025.2
        mu[199]=2034.1
        mu[200]=2060
        mu[201]=2068.9
        mu[202]=2077.8
        mu[203]=2086.7
        mu[204]=2095.6
        mu[205]=2104.5
        mu[206]=2113.4
        mu[207]=2122.3
        mu[208]=2131.2
        mu[209]=2140.1
        mu[210]=2166
        mu[211]=2174.9
        mu[212]=2183.8
        mu[213]=2192.7
        mu[214]=2201.6
        mu[215]=2210.5
        mu[216]=2219.4
        mu[217]=2228.3
        mu[218]=2237.2
        mu[219]=2246.1
        mu[220]=2271
        mu[221]=2279.9
        mu[222]=2288.8
        mu[223]=2297.7
        mu[224]=2306.6
        mu[225]=2315.5
        mu[226]=2324.4
        mu[227]=2333.3
        mu[228]=2342.2
        mu[229]=2351.1
        mu[230]=2377
        mu[231]=2385.9
        mu[232]=2394.8
        mu[233]=2403.7
        mu[234]=2412.6
        mu[235]=2421.5
        mu[236]=2430.4
        mu[237]=2439.3
        mu[238]=2448.2
        mu[239]=2457.1
        mu[240]=2481
        mu[241]=2489.9
        mu[242]=2498.8
        mu[243]=2507.7
        mu[244]=2516.6
        mu[245]=2525.5
        mu[246]=2534.4
        mu[247]=2543.3
        mu[248]=2552.2
        mu[249]=2561.1
        mu[250]=2589
        mu[251]=2597.9
        mu[252]=2606.8
        mu[253]=2615.7
        mu[254]=2624.6
        mu[255]=2633.5
        mu[256]=2642.4
        mu[257]=2651.3
        mu[258]=2660.2
        mu[259]=2669.1
        mu[260]=2693
        mu[261]=2701.9
        mu[262]=2710.8
        mu[263]=2719.7
        mu[264]=2728.6
        mu[265]=2737.5
        mu[266]=2746.4
        mu[267]=2755.3
        mu[268]=2764.2
        mu[269]=2773.1
        mu[270]=2797
        mu[271]=2805.9
        mu[272]=2814.8
        mu[273]=2823.7
        mu[274]=2832.6
        mu[275]=2841.5
        mu[276]=2850.4
        mu[277]=2859.3
        mu[278]=2868.2
        mu[279]=2877.1
        mu[280]=2902
        mu[281]=2910.9
        mu[282]=2919.8
        mu[283]=2928.7
        mu[284]=2937.6
        mu[285]=2946.5
        mu[286]=2955.4
        mu[287]=2964.3
        mu[288]=2973.2
        mu[289]=2982.1
        mu[290]=3005
        mu[291]=3013.9
        mu[292]=3022.8
        mu[293]=3031.7
        mu[294]=3040.6
        mu[295]=3049.5
        mu[296]=3058.4
        mu[297]=3067.3
        mu[298]=3076.2
        mu[299]=3085.1
        mu[300]=3110
        mu[301]=3118.9
        mu[302]=3127.8
        mu[303]=3136.7
        mu[304]=3145.6
        mu[305]=3154.5
        mu[306]=3163.4
        mu[307]=3172.3
        mu[308]=3181.2
        mu[309]=3190.1
        mu[310]=3213
        mu[311]=3221.9
        mu[312]=3230.8
        mu[313]=3239.7
        mu[314]=3248.6
        mu[315]=3257.5
        mu[316]=3266.4
        mu[317]=3275.3
        mu[318]=3284.2
        mu[319]=3293.1
        mu[320]=3315
        mu[321]=3323.9
        mu[322]=3332.8
        mu[323]=3341.7
        mu[324]=3350.6
        mu[325]=3359.5
        mu[326]=3368.4
        mu[327]=3377.3
        mu[328]=3386.2
        mu[329]=3395.1
        mu[330]=3415
        mu[331]=3423.9
        mu[332]=3432.8
        mu[333]=3441.7
        mu[334]=3450.6
        mu[335]=3459.5
        mu[336]=3468.4
        mu[337]=3477.3
        mu[338]=3486.2
        mu[339]=3495.1
        mu[340]=3516
        mu[341]=3524.9
        mu[342]=3533.8
        mu[343]=3542.7
        mu[344]=3551.6
        mu[345]=3560.5
        mu[346]=3569.4
        mu[347]=3578.3
        mu[348]=3587.2
        mu[349]=3596.1
        mu[350]=3617
        mu[351]=3625.4
        mu[352]=3633.8
        mu[353]=3642.2
        mu[354]=3650.6
        mu[355]=3659
        mu[356]=3667.4
        mu[357]=3675.8
        mu[358]=3684.2
        mu[359]=3692.6
        mu[360]=3715
        mu[361]=3723.4
        mu[362]=3731.8
        mu[363]=3740.2
        mu[364]=3748.6
        mu[365]=3757
        mu[366]=3765.4
        mu[367]=3773.8
        mu[368]=3782.2
        mu[369]=3790.6
        mu[370]=3813
        mu[371]=3821.4
        mu[372]=3829.8
        mu[373]=3838.2
        mu[374]=3846.6
        mu[375]=3855
        mu[376]=3863.4
        mu[377]=3871.8
        mu[378]=3880.2
        mu[379]=3888.6
        mu[380]=3910
        mu[381]=3918.4
        mu[382]=3926.8
        mu[383]=3935.2
        mu[384]=3943.6
        mu[385]=3952
        mu[386]=3960.4
        mu[387]=3968.8
        mu[388]=3977.2
        mu[389]=3985.6
        mu[390]=4007
        mu[391]=4015.4
        mu[392]=4023.8
        mu[393]=4032.2
        mu[394]=4040.6
        mu[395]=4049
        mu[396]=4057.4
        mu[397]=4065.8
        mu[398]=4074.2
        mu[399]=4082.6
        
        
        
        
        
        
        
        self.pShort = np.hstack((A, sigma, gamma, mu))        
        
        
#         # p = AA, AB, AC, muA, muB, muC, sigmaA, sigmaB, sigmaC, gamma, sGap, lGapA, lGapB, lGapC
#         lGapA = -1e-5
#         lGapB = 1e-2
#         lGapC = 15.
#         sGapA = 0.#-3.8e-5
#         sGapB = 0.#1e-2
#         sGapC = 8
#         AA  = 0.
#         AB = 0.
#         AC = 14000.
#         muA = 0.
#         muB = 0.
#         muC = 0.1
#         sigmaA = 0.
#         sigmaB = 0.
#         sigmaC = 2
#         firstGap = 25
#         self.pShort = [AA,AB,AC,muA,muB,muC,sigmaA,sigmaB,sigmaC,0,sGapA,sGapB,sGapC,lGapA,lGapB,lGapC, firstGap]
#         self.pShort = [3.82731509e+00,   
#                   2.52457594e+00, 
#                   2.10977540e+04,   
#                   8.20869367e+02,
#                   1.34588314e+00,  
#                   -3.69510633e+02,
#                   1.61574258e-02,
#                   -1.60114662e-01,
#                   2.05803744e+00,
#                   2.24741646e-01,   
#                   7.21478242e-06,  
#                   -1.64175390e+03,
#                   -2.45577948e+03,  
#                   -1.00000000e-05,   
#                   1.00000000e-02,   
#                   1.50000000e+01,
#                   -4.28664936e+02]

    def write_p400_A_mu(self, flatCurve):
         

        sigma = 1.9
        gamma = 0.2 #5.9997168955962937e-06
        
#         sigma = np.ones(400) * sigma
#         gamma = np.ones(400) * gamma
        
        A = np.ones(self.nFibres) * np.max(flatCurve) * np.sqrt(2* math.pi * sigma**2)

        mu = np.zeros(400)
        mu[0]=24.5
        mu[1]=32.501
        mu[2]=40.502
        mu[3]=48.503
        mu[4]=56.504
        mu[5]=64.505
        mu[6]=72.506
        mu[7]=80.507
        mu[8]=88.508
        mu[9]=96.509
        mu[10]=120
        mu[11]=128.1
        mu[12]=136.2
        mu[13]=144.3
        mu[14]=152.4
        mu[15]=160.5
        mu[16]=168.6
        mu[17]=176.7
        mu[18]=184.8
        mu[19]=192.9
        mu[20]=214
        mu[21]=222.1
        mu[22]=230.2
        mu[23]=238.3
        mu[24]=246.4
        mu[25]=254.5
        mu[26]=262.6
        mu[27]=270.7
        mu[28]=278.8
        mu[29]=286.9
        mu[30]=312
        mu[31]=320.05
        mu[32]=328.1
        mu[33]=336.15
        mu[34]=344.2
        mu[35]=352.25
        mu[36]=360.3
        mu[37]=368.35
        mu[38]=376.4
        mu[39]=384.45
        mu[40]=409
        mu[41]=417.4
        mu[42]=425.8
        mu[43]=434.2
        mu[44]=442.6
        mu[45]=451
        mu[46]=459.4
        mu[47]=467.8
        mu[48]=476.2
        mu[49]=484.6
        mu[50]=508
        mu[51]=516.4
        mu[52]=524.8
        mu[53]=533.2
        mu[54]=541.6
        mu[55]=550
        mu[56]=558.4
        mu[57]=566.8
        mu[58]=575.2
        mu[59]=583.6
        mu[60]=608
        mu[61]=616.4
        mu[62]=624.8
        mu[63]=633.2
        mu[64]=641.6
        mu[65]=650
        mu[66]=658.4
        mu[67]=666.8
        mu[68]=675.2
        mu[69]=683.6
        mu[70]=707
        mu[71]=715.4
        mu[72]=723.8
        mu[73]=732.2
        mu[74]=740.6
        mu[75]=749
        mu[76]=757.4
        mu[77]=765.8
        mu[78]=774.2
        mu[79]=782.6
        mu[80]=810
        mu[81]=818.4
        mu[82]=826.8
        mu[83]=835.2
        mu[84]=843.6
        mu[85]=852
        mu[86]=860.4
        mu[87]=868.8
        mu[88]=877.2
        mu[89]=885.6
        mu[90]=910
        mu[91]=918.4
        mu[92]=927.3
        mu[93]=936.2
        mu[94]=945.1
        mu[95]=954
        mu[96]=962.9
        mu[97]=971.8
        mu[98]=980.7
        mu[99]=989.6
        mu[100]=1013
        mu[101]=1021.9
        mu[102]=1030.8
        mu[103]=1039.7
        mu[104]=1048.6
        mu[105]=1057.5
        mu[106]=1066.4
        mu[107]=1075.3
        mu[108]=1084.2
        mu[109]=1093.1
        mu[110]=1117
        mu[111]=1125.9
        mu[112]=1134.8
        mu[113]=1143.7
        mu[114]=1152.6
        mu[115]=1161.5
        mu[116]=1170.4
        mu[117]=1179.3
        mu[118]=1188.2
        mu[119]=1197.1
        mu[120]=1219
        mu[121]=1227.9
        mu[122]=1236.8
        mu[123]=1245.7
        mu[124]=1254.6
        mu[125]=1263.5
        mu[126]=1272.4
        mu[127]=1281.3
        mu[128]=1290.2
        mu[129]=1299.1
        mu[130]=1324
        mu[131]=1332.9
        mu[132]=1341.8
        mu[133]=1350.7
        mu[134]=1359.6
        mu[135]=1368.5
        mu[136]=1377.4
        mu[137]=1386.3
        mu[138]=1395.2
        mu[139]=1404.1
        mu[140]=1428
        mu[141]=1436.9
        mu[142]=1445.8
        mu[143]=1454.7
        mu[144]=1463.6
        mu[145]=1472.5
        mu[146]=1481.4
        mu[147]=1490.3
        mu[148]=1499.2
        mu[149]=1508.1
        mu[150]=1533
        mu[151]=1541.9
        mu[152]=1550.8
        mu[153]=1559.7
        mu[154]=1568.6
        mu[155]=1577.5
        mu[156]=1586.4
        mu[157]=1595.3
        mu[158]=1604.2
        mu[159]=1613.1
        mu[160]=1637
        mu[161]=1645.9
        mu[162]=1654.8
        mu[163]=1663.7
        mu[164]=1672.6
        mu[165]=1681.5
        mu[166]=1690.4
        mu[167]=1699.3
        mu[168]=1708.2
        mu[169]=1717.1
        mu[170]=1744
        mu[171]=1752.9
        mu[172]=1761.8
        mu[173]=1770.7
        mu[174]=1779.6
        mu[175]=1788.5
        mu[176]=1797.4
        mu[177]=1806.3
        mu[178]=1815.2
        mu[179]=1824.1
        mu[180]=1849
        mu[181]=1857.9
        mu[182]=1866.8
        mu[183]=1875.7
        mu[184]=1884.6
        mu[185]=1893.5
        mu[186]=1902.4
        mu[187]=1911.3
        mu[188]=1920.2
        mu[189]=1929.1
        mu[190]=1954
        mu[191]=1962.9
        mu[192]=1971.8
        mu[193]=1980.7
        mu[194]=1989.6
        mu[195]=1998.5
        mu[196]=2007.4
        mu[197]=2016.3
        mu[198]=2025.2
        mu[199]=2034.1
        mu[200]=2060
        mu[201]=2068.9
        mu[202]=2077.8
        mu[203]=2086.7
        mu[204]=2095.6
        mu[205]=2104.5
        mu[206]=2113.4
        mu[207]=2122.3
        mu[208]=2131.2
        mu[209]=2140.1
        mu[210]=2166
        mu[211]=2174.9
        mu[212]=2183.8
        mu[213]=2192.7
        mu[214]=2201.6
        mu[215]=2210.5
        mu[216]=2219.4
        mu[217]=2228.3
        mu[218]=2237.2
        mu[219]=2246.1
        mu[220]=2271
        mu[221]=2279.9
        mu[222]=2288.8
        mu[223]=2297.7
        mu[224]=2306.6
        mu[225]=2315.5
        mu[226]=2324.4
        mu[227]=2333.3
        mu[228]=2342.2
        mu[229]=2351.1
        mu[230]=2377
        mu[231]=2385.9
        mu[232]=2394.8
        mu[233]=2403.7
        mu[234]=2412.6
        mu[235]=2421.5
        mu[236]=2430.4
        mu[237]=2439.3
        mu[238]=2448.2
        mu[239]=2457.1
        mu[240]=2481
        mu[241]=2489.9
        mu[242]=2498.8
        mu[243]=2507.7
        mu[244]=2516.6
        mu[245]=2525.5
        mu[246]=2534.4
        mu[247]=2543.3
        mu[248]=2552.2
        mu[249]=2561.1
        mu[250]=2589
        mu[251]=2597.9
        mu[252]=2606.8
        mu[253]=2615.7
        mu[254]=2624.6
        mu[255]=2633.5
        mu[256]=2642.4
        mu[257]=2651.3
        mu[258]=2660.2
        mu[259]=2669.1
        mu[260]=2693
        mu[261]=2701.9
        mu[262]=2710.8
        mu[263]=2719.7
        mu[264]=2728.6
        mu[265]=2737.5
        mu[266]=2746.4
        mu[267]=2755.3
        mu[268]=2764.2
        mu[269]=2773.1
        mu[270]=2797
        mu[271]=2805.9
        mu[272]=2814.8
        mu[273]=2823.7
        mu[274]=2832.6
        mu[275]=2841.5
        mu[276]=2850.4
        mu[277]=2859.3
        mu[278]=2868.2
        mu[279]=2877.1
        mu[280]=2902
        mu[281]=2910.9
        mu[282]=2919.8
        mu[283]=2928.7
        mu[284]=2937.6
        mu[285]=2946.5
        mu[286]=2955.4
        mu[287]=2964.3
        mu[288]=2973.2
        mu[289]=2982.1
        mu[290]=3005
        mu[291]=3013.9
        mu[292]=3022.8
        mu[293]=3031.7
        mu[294]=3040.6
        mu[295]=3049.5
        mu[296]=3058.4
        mu[297]=3067.3
        mu[298]=3076.2
        mu[299]=3085.1
        mu[300]=3110
        mu[301]=3118.9
        mu[302]=3127.8
        mu[303]=3136.7
        mu[304]=3145.6
        mu[305]=3154.5
        mu[306]=3163.4
        mu[307]=3172.3
        mu[308]=3181.2
        mu[309]=3190.1
        mu[310]=3213
        mu[311]=3221.9
        mu[312]=3230.8
        mu[313]=3239.7
        mu[314]=3248.6
        mu[315]=3257.5
        mu[316]=3266.4
        mu[317]=3275.3
        mu[318]=3284.2
        mu[319]=3293.1
        mu[320]=3315
        mu[321]=3323.9
        mu[322]=3332.8
        mu[323]=3341.7
        mu[324]=3350.6
        mu[325]=3359.5
        mu[326]=3368.4
        mu[327]=3377.3
        mu[328]=3386.2
        mu[329]=3395.1
        mu[330]=3415
        mu[331]=3423.9
        mu[332]=3432.8
        mu[333]=3441.7
        mu[334]=3450.6
        mu[335]=3459.5
        mu[336]=3468.4
        mu[337]=3477.3
        mu[338]=3486.2
        mu[339]=3495.1
        mu[340]=3516
        mu[341]=3524.9
        mu[342]=3533.8
        mu[343]=3542.7
        mu[344]=3551.6
        mu[345]=3560.5
        mu[346]=3569.4
        mu[347]=3578.3
        mu[348]=3587.2
        mu[349]=3596.1
        mu[350]=3617
        mu[351]=3625.4
        mu[352]=3633.8
        mu[353]=3642.2
        mu[354]=3650.6
        mu[355]=3659
        mu[356]=3667.4
        mu[357]=3675.8
        mu[358]=3684.2
        mu[359]=3692.6
        mu[360]=3715
        mu[361]=3723.4
        mu[362]=3731.8
        mu[363]=3740.2
        mu[364]=3748.6
        mu[365]=3757
        mu[366]=3765.4
        mu[367]=3773.8
        mu[368]=3782.2
        mu[369]=3790.6
        mu[370]=3813
        mu[371]=3821.4
        mu[372]=3829.8
        mu[373]=3838.2
        mu[374]=3846.6
        mu[375]=3855
        mu[376]=3863.4
        mu[377]=3871.8
        mu[378]=3880.2
        mu[379]=3888.6
        mu[380]=3910
        mu[381]=3918.4
        mu[382]=3926.8
        mu[383]=3935.2
        mu[384]=3943.6
        mu[385]=3952
        mu[386]=3960.4
        mu[387]=3968.8
        mu[388]=3977.2
        mu[389]=3985.6
        mu[390]=4007
        mu[391]=4015.4
        mu[392]=4023.8
        mu[393]=4032.2
        mu[394]=4040.6
        mu[395]=4049
        mu[396]=4057.4
        mu[397]=4065.8
        mu[398]=4074.2
        mu[399]=4082.6

        
        self.p400_A_mu = np.hstack((A, sigma, gamma, mu))        

    def write_p400_sigma_gamma(self, flatCurve):
         

        sigma = 1.9
        gamma = 0.2 #5.9997168955962937e-06
        
        sigma = np.ones(400) * sigma
        gamma = np.ones(400) * gamma
        
        self.p400_sigma_gamma = np.hstack((sigma, gamma))        
    
    def find_residuals(self, p, flatCurve, output=True):
        self.reps += 1
        if len(p)==len(self.pShort):
            model = self.find_pShort_curve(p, flatCurve)[0]
            
        elif len(p)==len(self.p400_A_mu):
            model = self.find_p400_A_mu_curve(p, flatCurve)[0]

        elif len(p)==len(self.p400_sigma_gamma):
            model = self.find_p400_sigma_gamma_curve(p, flatCurve)[0]

        elif len(p)==len(self.p10):
            model = self.find_p10_curve(p, flatCurve, output=output)[0]

        diff = flatCurve-model

        
        #just output to see progress
        if ((self.reps in [0,1,2,3,4,5,6,7,8,10,20,50,100, 500, 1000]) and (output==True)):
            print 'total diff = ', np.sum(diff), 'iterations: ', self.reps
#             print p[:5], 'A:', np.average(p[-2*self.nFibres:-self.nFibres]), ' mu:', np.average(p[-self.nFibres:]) 
            print ' '
            
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

        A = p[0:self.nFibres]
        sigma = p[self.nFibres]
        gamma = p[self.nFibres+1]
#         mu = p[-self.nFibres:]
                
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
    
    def read_full_image_spatial(self, profile, nColumns):

#         self.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
        self.base_dir = ''
        self.scienceFile = self.base_dir + '10nov20044.fits'
        self.biasFile = self.base_dir + 'BIAScombined2.fits'
        self.profile = profile
        self.out_dir = self.base_dir + ''
        self.flatFile = self.base_dir + '10nov20044.fits'
        self.tramLinefile = self.base_dir + '10nov20044tlm.fits'
        self.deadFibres = [41,66, 91, 141,191,241, 250, 294, 304, 341, 391]
        self.nFibres = 400
        self.open_files()
        self.flatIm_b = self.bias_subtract_from_overscan(self.flatIm, range(-45,-5)) #bias subtract using the last 40 columns(overscan region)

        if profile=='voigt':
            fileOutSpa = open(self.out_dir + 'spatialV2.txt','w')
        elif profile=='gaussian':
            fileOutSpa = open(self.out_dir + 'spatialG2.txt','w')
        fileOutSpa.write('fibre, y_cent, A, sigma, gamma, mu \n')
        
        for C in range(0, self.imWidth,self.imWidth/nColumns):
            print time.asctime()
            flatCurve = self.flatIm_b[:,C]       
            
            #find A        
            self.write_p400_A_mu(flatCurve)
            self.write_mu_from_tlm(self.p400_A_mu,C)
            p = self.fit_p400_A_mu(flatCurve)
            self.A = p[0][:self.nFibres]
            
            #find sigma, gamma
            self.write_p400_sigma_gamma(flatCurve)
            p = self.fit_p400_sigma_gamma(flatCurve)
            
            sigma = p[0][:self.nFibres]
            gamma = p[0][self.nFibres:]
            for F in range(400):
                outString = (str(F+1) + ', ' + 
                            str(C) + ', ' + 
                            str(self.A[F]) + ', ' + 
                            str(sigma[F]) + ', ' + 
                            str(gamma[F]) + ', ' + 
                            str(self.mu[F]) + ', ' + 
                            '\n' )
                fileOutSpa.write(outString)

    def write_mu_from_tlm(self, p, column):
        
        self.p400_A_mu = p[:-self.nFibres]
        self.mu = self.tlmIm[:,column]-0.5
        
    def read_psf_data(self, profile):
        import matplotlib
        import matplotlib.cm as cm
        import matplotlib.mlab as mlab
        import matplotlib.pyplot as plt
        
    #     from scipy.interpolate import interp1d
    #     from scipy.interpolate import interp2d
    #     import pyfits as pf
        
        self.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
        self.out_dir = self.base_dir + 'output/'
        if profile=='voigt':
            dataFile = open(self.out_dir + 'spatialV1.txt','r')
        elif profile=='gaussian':
            dataFile = open(self.out_dir + 'spatialG1.txt','r')

        a = np.loadtxt(dataFile, skiprows=1, delimiter = ',' , usecols=[0,1,2,3,4,5]).transpose()
    

        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        
        X, Y = np.meshgrid(np.unique(a[1]), np.unique(a[0]))
        Z = a[3].reshape(X.shape[1],X.shape[0])
        Z = Z.transpose()



        # Or you can use a colormap to specify the colors; the default
        # colormap will be used for the contour lines
        plt.figure()
        im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                        cmap=cm.gray, extent=(0,400,0,400))
        levels = np.arange(1.4, 2.2, 0.1)
        CS = plt.contour(Z, levels,
                         origin='lower',
                         cmap=cm.jet, 
                         linewidths=2,
                         extent=(0,400,0,400))
        
        #Thicken the zero contour.
#         zc = CS.collections[6]
#         plt.setp(zc, linewidth=4)
        
#         plt.clabel(CS, levels[1::2],  # label every second level
#                    inline=1,
#                    fmt='%1.1f',
#                    fontsize=14)
        
        # make a colorbar for the contour lines
#         CB = plt.colorbar(CS, shrink=0.8, extend='both')
        
        plt.title('Lines with colorbar')
        #plt.hot()  # Now change the colormap for the contour lines and colorbar
#         plt.flag()
        
        # We can still add a colorbar for the image, too.
#         CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)
        
        # This makes the original colorbar look a bit out of place,
        # so let's improve its position.
        
#         l,b,w,h = plt.gca().get_position().bounds
#         ll,bb,ww,hh = CB.ax.get_position().bounds
#         CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
        
        
        plt.show()






    
        #scatter plot
    #     plt.scatter(a[0],a[1],s = a[2])                                                                                      
    #     plt.show()
        

        
        
        #binned plot
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
#             plt.xlabel('X-Pixels')
#             plt.ylabel('Std. Dev. (px)')
#             plt.legend()
#         plt.savefig(base_dir + 'spa_' + str(cam) + '_y_' + method +'.png')
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

    