import pyfits as pf
import scipy.constants as const
import numpy as np
import pylab as plt
import os
import shutil
import subprocess 

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

    
    
    