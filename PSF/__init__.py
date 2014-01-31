################################################
#deprecated
#part of the HERMES module now


import numpy as np
import pylab as plt
import pyfits as pf
import scipy.optimize as opt

# import wsm
from pymodelfit.builtins import VoigtModel
from pymodelfit.builtins import GaussianModel
from pymodelfit import get_model_instance
import math 


class HERMES():
    
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    sexParamFile = base_dir + 'HERMES.sex'
    outputFileName = base_dir + 'out.txt'
    scienceFile = base_dir + '10nov40045.fits'
    biasFile = base_dir + 'BIAScombined4.fits'
    flatFile = base_dir + '10nov10044.fits'
    tramLinefile = base_dir + '10nov10044tlm.fits'
    nFibres = 10

    def __init__(self):
        self.open_files()
        self.deadFibres = [41,66, 91, 141,191,241, 294, 307, 341, 391]
        self.nFibres = 400
        self.nBundles = 40
        self.pShort = []
        self.p400 = []
        
    def open_files(self):
        
         #Read, calibrate and create im object
        self.bias = pf.open(self.biasFile)
        self.tlm = pf.open(self.tramLinefile)
        self.flat = pf.open(self.flatFile)
        self.science = pf.open(self.scienceFile)

        self.imWidth = self.science[0].header['NAXIS1']
        self.imHeight = self.science[0].header['NAXIS2']
        self.im = self.science[0].data     
        biasIm = self.bias[0].data
        self.flatIm = self.flat[0].data
        self.im -= biasIm       
        self.flatIm -= biasIm

    def fit_multi_PSF(self):
       
        flatCurve = self.flatIm[:, 10]
#         flatCurve = self.flatIm[:110, 100]
        self.multi_fit(flatCurve)

    def fit_A_PSF(self):
        self.reps = 0
        flatCurve = self.flatIm[:, 10]
        self.A_fit(flatCurve)

    def A_fit(self, flatCurve):
        
        sigma = 2.
        gamma = 2.24741646e-01
        firstGap = 24.5
        lGap = 17.
        sGap = 8.0001
        mu = np.ones(self.nFibres) * 0.
        A = np.ones(self.nFibres) * np.average(flatCurve) * math.sqrt(2* math.pi * sigma**2)
        print A[0]
        self.p400 = np.hstack(( sigma, gamma, firstGap, lGap, sGap, A, mu))
        
        factorTry = 1
        diagTry = np.ones(len(self.p400))
        maxfev = int(1e3) # of max tries
        # diagonals <> 1 to follow:
        diagTry[5] = 0.1
        diagTry[-self.nFibres:] = np.ones(self.nFibres) * 0.1
        p = opt.leastsq(self.find_residuals, self.p400, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)

        finalCurve, onePSF = self.find_A_curve(p[0],flatCurve, output = False)
        
        plt.plot(flatCurve, label = 'flat')
        plt.plot(finalCurve,label = 'model')
        plt.plot(onePSF,label = 'Single PSF')
        plt.legend()
        plt.show()
        
        plt.plot(flatCurve/max(flatCurve), label = 'flat')
        plt.plot(finalCurve/max(finalCurve),label = 'model')
        plt.plot(onePSF/max(onePSF),label = 'Single PSF')
        plt.legend()
        plt.show()
              
    def multi_fit(self, flatCurve):
        
        # p = AA, AB, AC, muA, muB, muC, sigmaA, sigmaB, sigmaC, gamma, sGap, lGapA, lGapB, lGapC
        lGapA = -1e-5
        lGapB = 1e-2
        lGapC = 15.
        sGapA = 0.#-3.8e-5
        sGapB = 0.#1e-2
        sGapC = 8
        AA  = 0.
        AB = 0.
        AC = 14000.
        muA = 0.
        muB = 0.
        muC = 0.1
        sigmaA = 0.
        sigmaB = 0.
        sigmaC = 2
        firstGap = 25
        self.pShort = [AA,AB,AC,muA,muB,muC,sigmaA,sigmaB,sigmaC,0,sGapA,sGapB,sGapC,lGapA,lGapB,lGapC, firstGap]
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
        factorTry = 1
        diagTry = np.ones(len(self.pShort))
        diagTry[0] = 200
        diagTry[1] = 200
        diagTry[2] = 0.1
        diagTry[3] = 1
        diagTry[4] = 30
        diagTry[5] = 1
        diagTry[6] = 1
        diagTry[7] = 1
        diagTry[8] = 0.1
        diagTry[9] = 0.1
        
        maxfev = int(1e3)
        
        p = opt.leastsq(self.find_residuals, self.pShort, full_output=True, args = flatCurve, factor = factorTry, diag = diagTry, maxfev = maxfev)
#         print p
#         print p[0]
         
        plt.subplot(111)
#         plt.plot(flatCurve)
#         plt.plot(self.find_curve(p[0]))
        print "Final info:"
        finalCurve, onePSF = self.find_curve(p[0],flatCurve, output = True)
        plt.plot(flatCurve/max(flatCurve), label = 'flat')
        plt.plot(finalCurve/max(finalCurve),label = 'model')
        plt.plot(onePSF/max(onePSF),label = 'Single PSF')
        plt.legend()
        plt.show()

    def find_residuals(self, p, flatCurve):
        self.reps += 1
        if len(p)==len(self.pShort):
            model = self.find_curve(p, flatCurve)[0]
            
        elif len(p)==len(self.p400):
            model = self.find_A_curve(p, flatCurve)[0]

        diff = flatCurve-model
        #diff_norm = flatCurve/max(flatCurve)-model/max(model)

        if self.reps in [1,2,3,4,5,6,7,8,10,20,50,100, 500, 1000]:
            print 'total diff=', np.sum(diff), 'iterations: ', self.reps
            print p[:5], 'A:', np.average(p[-2*self.nFibres:-self.nFibres]), ' mu:', np.average(p[-self.nFibres:]) 
            print ' '
        return diff
 
    def find_A_curve(self, p, flatCurve, output = False, psfFibre = 5):
         
        A = p[-2*self.nFibres:-self.nFibres]
        mu = p[-self.nFibres:]
        sigma = p[0]
        gamma = p[1]
        firstGap = p[2]
        lGap = p[3]
        sGap = p[4]

        a = get_model_instance('Voigt')
        model = np.zeros(len(flatCurve))
        onePSF = np.zeros(len(flatCurve))
        
        for i in range(self.nBundles):
            lGapTot = lGap * i
            for j in range(10):
                fibre = i*10 + j + 1
                if fibre not in self.deadFibres: 
                    sGapTot = sGap * (fibre - 1)
                    thisA = A[fibre-1]
                    thisMu = mu[fibre-1]
                    gap = firstGap + sGapTot + lGapTot
                    thisCurve = a.f(np.arange(len(flatCurve))-gap, thisA,sigma,gamma,thisMu)
                    if ((np.sum(thisCurve) >0) and (output == True)):
                        print thisA,sigma, gamma, mu
                    if fibre==psfFibre:
                        onePSF = thisCurve
                    model += thisCurve
        
#         plt.plot(flatCurve)
#         plt.plot(model)
#         plt.show()
        return model, onePSF
  
    def find_curve(self, p, flatCurve, output = False):

        AA = p[0]
        AB = p[1]
        AC = p[2]
        muA = p[3]
        muB = p[4]
        muC = p[5]
        sigmaA = p[6]
        sigmaB = p[7]
        sigmaC = p[8]
        gamma = p[9]
        sGapA = p[10]
        sGapB = p[11]
        sGapC = p[12]
        lGapA = p[13]
        lGapB = p[14]
        lGapC = p[15]
        firstGap = p[16]
#         print AA, AB, AC
#         print muA, muB, muC
#         print sigmaA, sigmaB, sigmaC, gamma
#         print lGapA, lGapB, lGapC 

        a = get_model_instance('Voigt')
        model = np.zeros(len(flatCurve))

        for i in range(self.nBundles):
            lGapTotList = np.arange(i)*10+1
            lGapTot = np.sum(lGapA*lGapTotList**2 + lGapB*lGapTotList + lGapC)
            for j in range(10):
                fibre = i*10 + j + 1
                if fibre not in self.deadFibres: 
                    sGapTotList = np.arange(fibre-1)
                    sGapTot = np.sum(sGapA*sGapTotList**2 + sGapB*sGapTotList + sGapC)
                    A = AA*fibre**2 + AB*fibre + AC
                    sigma = sigmaA*fibre**2 + sigmaB*fibre + sigmaC
                    mu = muA*fibre**2 + muB*fibre + muC
#                     lGap = lGapA*fibre**2 + lGapB*fibre + lGapC
#                     sGap = sGapA*fibre**2 + sGapB*fibre + sGapC
                    gap = firstGap + sGapTot + lGapTot
#                     thisCurve = a.f(np.arange(self.imHeight)-gap, A,sigma,gamma,mu)
                    thisCurve = a.f(np.arange(len(flatCurve))-gap, A,sigma,gamma,mu)
                    if ((np.sum(thisCurve) >0) and (output == True)):
                            print A,sigma, gamma, mu
                    if j==5:
                        onePSF = thisCurve
                    model += thisCurve
                    
        return model, onePSF
        
#         if False==True:
#             x = np.arange(1,400)
#             plt_A = AA*x**2 + AB*x + AC
#             plt_sigma = sigmaA*x**2 + sigmaB*x + sigmaC
#             plt_mu = muA*x**2 + muB*x + muC
#             plt.subplot(221)
#             plt.plot(x,plt_A, label = 'A')
#             plt.legend()
#             plt.subplot(222)
#             plt.plot(x,plt_sigma, label = 'sigma')
#             plt.plot(x,plt_mu, label = 'mu')
#             plt.legend()
#             plt.show()

#                     plt.plot(thisCurve)
#         plt.plot(model)
#         plt.show()
    
        
#         self.find_residuals(flatCurve,N,0,0,1,0,self.imHeight/N,self.imHeight/N,0,0,1,1)
        
        
#         for j in range(0,self.imWidth/10,10): #range(areaRange):
#             for i in self.tlm[0].data[:,j]:
#             
#                 area = self.im[i-sigTest:i+sigTest, j:j+10]
# 
#                 if (plotPoints==True):    
#                     plt.subplot(121)
#                     plt.imshow(area, origin='lower')
#                     plt.set_cmap(cm.Greys_r) 
#         
#                 #Spatial direction
#                 line = im[i-sigTest:i+sigTest, j]
#                 line[np.isnan(line)] = 0
#                 b = np.arange(int(i-sigTest), int(i+sigTest))
#                 if method=='voigt':
#                     a = get_model_instance('Voigt')
#                     a.parvals = (max(line),1,1,np.average(b))
#                     a.parvals = VoigtModel.fitData(a, b, line)
#                 elif method=='gaussian':
#                     a = get_model_instance('gaussian')
#                     a.parvals = (max(line),1,np.average(b))
#                     a.parvals = GaussianModel.fitData(a, b, line)
#                 c = np.arange(min(b),max(b),0.01)
#     
#                 #output and plot
#                 outString = (str(j*10) + ', ' + 
#                             str(i) + ', ' + 
#                             str(a.FWHM) + ', ' + 
#                             str(a.pardict['A']) + ', ' + 
#                             str(a.pardict['mu']) + ', ' + 
#                             str(a.pardict['sig']) + ', ' + 
#                             '\n' )
#                 if (a.chi2Data()[1]<100): fileOutSpa.write(outString)
#                 print 'Spatial: ' + str(a.FWHM), a.chi2Data()
#                 if (plotPoints==True):
#                     plt.subplot(122)
#                     plt.title('Spatial ' + str(a.FWHM) + ', (' + str(a.chi2Data()[1]) + ')') # + str(a.chi2Data()[1]) + ')')
#                     plt.plot(b, line)
#                     plt.plot(c, a(c))
#                     plt.show()      

    def read_psf_data_spatial(method, cam):
    #     from scipy.interpolate import interp1d
    #     from scipy.interpolate import interp2d
    #     import pyfits as pf
        
        base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
        dataFile = base_dir + 'psf_HERMES1_spa_' + method + str(cam) +'.txt'
        
        a = np.loadtxt(dataFile,skiprows=1, delimiter = ',' , usecols = [0,1,2,3,4,5]).transpose()
        a = a.transpose()[a[5]<7].transpose()
        a = a.transpose()[a[5]>1].transpose()
    
        #scatter plot
    #     plt.scatter(a[0],a[1],s = a[2])                                                                                      
    #     plt.show()
        
        #binned plot
        rows = 4112
        cols = 4146
        bins = 4
        #bin rows
        for i in range(bins):
            binLimits = (rows/bins*i,rows/bins*(i+1))
            mapBin = ((a[1]>binLimits[0]) & (a[1]<=binLimits[1]))
            b = a[0,mapBin], a[5,mapBin]
            c = np.bincount(b[0].astype(int), weights=b[1])
            d = np.bincount(b[0].astype(int))
            px = np.arange(len(d))
            px = px[d>0]
            c = c[d>0]
            d = d[d>0]
            c = c/d
            plt.plot(px,c, label= ' Range (y) = ' + str(binLimits))
            plt.title('Spatial PSF - y-binned - ' + method)
            plt.xlabel('X-Pixels')
            plt.ylabel('Std. Dev. (px)')
            plt.legend()
        plt.savefig(base_dir + 'spa_' + str(cam) + '_y_' + method +'.png')
        plt.show()
    
        #bin cols
        for i in range(bins):
            binLimits = (cols/bins*i,cols/bins*(i+1))
            mapBin = ((a[0]>binLimits[0]) & (a[0]<=binLimits[1]))
            b = a[1,mapBin], a[5,mapBin]
            c = np.bincount(b[0].astype(int), weights=b[1])
            d = np.bincount(b[0].astype(int))
            px = np.arange(len(d))
            px = px[d>0]
            c = c[d>0]
            d = d[d>0]
            c = c/d
            plt.plot(px,c, label= ' Range (x) = ' + str(binLimits))
            plt.title('Spatial PSF - x-binned - ' + method)
            plt.xlabel('Y-Pixels')
            plt.ylabel('Std. Dev. (px)')
            plt.legend()
        plt.savefig(base_dir + 'spa_' + str(cam) + '_x_' + method +'.png')
        plt.show()
    
    
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

    
    
def do_all_plots():
    
    for cam in [1,2,3,4]:
        for method in ['Gauss', 'Voigt']:
            read_psf_data_spatial(method, cam)
            
    for cam in [1,2,3,4]:
        for method in ['Gauss', 'Voigt']:
            read_psf_data_spectral(method, cam)
 

def find_PSF_arc():
    import pyfits as pf
    import wsm
    from scipy import interpolate
    from pymodelfit.builtins import VoigtModel
    from pymodelfit.builtins import GaussianModel
    from pymodelfit import get_model_instance
    from scipy.stats import chisquare
    
    #Variable declaration
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    sexParamFile = base_dir + 'HERMES.sex'
    outputFileName = base_dir + 'out.txt'
    scienceFile = base_dir + '10nov40045.fits'
    biasFile = base_dir + 'BIAScombined4.fits'
    tramLinefile = base_dir + '10nov10044tlm.fits'
    method = 'gaussian'
    e = 2.718281
    areaRange = 9 
        
    #booleans
    plotMap = 0
    plotPoints = 0
    createPoints = 0
    calculateSpatial = 1
    calculateSpectral = 1
    
    #Read, calibrate and create im object
    bias = pf.open(biasFile)
    tlm = pf.open(tramLinefile)
    biasIm = bias[0].data
    hdulist = pf.open(scienceFile)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    im = hdulist[0].data     
    im = im-biasIm
    

    
    if createPoints==True: wsm.wt.ia.analyse_image_sex(scienceFile, sexParamFile, outputFileName)
 
    #Loads from calibration output file
    imageMapX, imageMapY, image_map_sigx, image_map_sigy   = wsm.wt.ia.load_image_map_sex(outputFileName)
    if imageMapX==[]: return
    print 'Points found: ' + str(len(imageMapX))
    imageMapX = imageMapX[image_map_sigx<8]
    imageMapY = imageMapY[image_map_sigx<8]
    image_map_sigy = image_map_sigy[image_map_sigx<8]
    image_map_sigx = image_map_sigx[image_map_sigx<8]
    imageMapX = imageMapX[image_map_sigx>1]
    imageMapY = imageMapY[image_map_sigx>1]
    image_map_sigy = image_map_sigy[image_map_sigx>1]
    image_map_sigx = image_map_sigx[image_map_sigx>1]
    imageMapX = imageMapX[image_map_sigy<8]
    imageMapY = imageMapY[image_map_sigy<8]
    image_map_sigx = image_map_sigx[image_map_sigy<8]
    image_map_sigy = image_map_sigy[image_map_sigy<8]
    imageMapX = imageMapX[image_map_sigy>1]
    imageMapY = imageMapY[image_map_sigy>1]
    image_map_sigx = image_map_sigx[image_map_sigy>1]
    image_map_sigy = image_map_sigy[image_map_sigy>1]
    print 'Points found after cleanup: ' + str(len(imageMapX))

 
    if (plotMap==True):    
        im2=im
#         im2 = np.log10(im)
    
#         im2 = wsm.wt.ic.normalise_image(im2) 
        im2 = np.sqrt(im2)
#         im2 = np.sqrt(im2)
     
         
        plt.imshow(im2,origin='lower')
        plt.set_cmap(cm.Greys_r)
         
        plt.scatter(imageMapX, imageMapY ,s=10, color="green" , marker='o', alpha = 1)
         
        plt.show()
        
    #Open output files
#     fileOutSpa = open(base_dir + 'psf_HERMES1_spa_Voigt_arc.txt','w')
    if method=='voigt':
        fileOutSpe = open(base_dir + 'psf_HERMES1_spe_Voigt4.txt','w')
    elif method=='gaussian':
        fileOutSpe = open(base_dir + 'psf_HERMES1_spe_Gauss4.txt','w')
#     fileOutSpa.write('x_cent, y_cent, FWHM, A, mu, sigma, gamma \n')
    fileOutSpe.write('x_cent, y_cent, FWHM, A, mu, sigma, gamma \n')
#     
#     imageMapY = np.arange(25,97,8) 
#     imageMapY = np.tile(imageMapY,40)
#     offset1 = np.arange(imHeight-25,96)
#     offsetAll = np.repeat(offset1, 10)
#     imageMapY = imageMapY + offsetAll
    
    for i in range(len(imageMapX)):
        #initial parameters and area plot
        regionHalfWidth = image_map_sigx[i]
        regionhalfHeight = image_map_sigy[i] * 1.1
        rangeX = np.arange(int(imageMapX[i]-regionHalfWidth),int(imageMapX[i]+regionHalfWidth))
        rangeX = rangeX[rangeX>=0]
        if len(rangeX)>3:
            rangeY = range(int(imageMapY[i]-regionhalfHeight),int(imageMapY[i]+regionhalfHeight))
            area = im[int(min(rangeY)):int(max(rangeY)), int(min(rangeX)):int(max(rangeX))]
#             print ''
#             print area.shape
            if (plotPoints==True):     
                plt.subplot(221)
                plt.imshow(area, origin='lower')
                plt.set_cmap(cm.Greys_r) 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
            #Spectral direction
            if (calculateSpectral==True):
                line = im[int(imageMapY[i]), min(rangeX):max(rangeX)+1]
                line[np.isnan(line)] = 0
#                 print line
#                 print rangeX
                if method=='voigt':
                    a = get_model_instance('Voigt')
                    a.parvals = (max(line),1,1,np.average(rangeX))
                elif method=='gaussian':
                    a = get_model_instance('gaussian')
                    a.parvals = (max(line),1,np.average(rangeX))
                    
                b = np.arange(min(rangeX),max(rangeX),0.01)            
                if method=='voigt':
                    a.parvals = VoigtModel.fitData(a, rangeX, line)
                elif method=='gaussian':
                    a.parvals = GaussianModel.fitData(a, rangeX, line)
    
                #Output and plot
                outString = (str(imageMapX[i]) + ', ' + 
                            str(imageMapY[i]) + ', ' + 
                            str(a.FWHM) + ', ' + 
                            str(a.pardict['A']) + ', ' + 
                            str(a.pardict['mu']) + ', ' + 
                            str(a.pardict['sig']) + ', ' + 
    #                             str(a.pardict['gamma']) + 
                            '\n' )
                fileOutSpe.write(outString)
                print 'Spectral: ' + str(a.FWHM), a.chi2Data()[1]#, a.pardict
                if (plotPoints==True): 
                    plt.subplot(223)
                    plt.title('Spectral')
                    plt.plot(rangeX, line)
                    plt.plot(b, a(b))
                    plt.show()
    
#             #Spatial direction
#             line = im[min(rangeY):max(rangeY)+1,int(imageMapX[i])]
#             line[np.isnan(line)] = 0
#             a = get_model_instance('Voigt')
#             b = np.arange(min(rangeY),max(rangeY),0.01)
#             a.parvals = (max(line),0.5,0.5,np.average(rangeY))
#             a.parvals = VoigtModel.fitData(a, rangeY, line)
# 
#             #output and plot
#             outString = (str(imageMapX[i]) + ', ' + 
#                         str(imageMapY[i]) + ', ' + 
#                         str(a.FWHM) + ', ' + 
#                         str(a.pardict['A']) + ', ' + 
#                         str(a.pardict['mu']) + ', ' + 
#                         str(a.pardict['sig']) + ', ' + 
#                         str(a.pardict['gamma']) + 
#                         '\n' )
#             if (a.chi2Data()[1]<100): fileOutSpa.write(outString)
# 
#             if (plotPoints==True):
#                 plt.subplot(222)
#                 plt.title('Spatial ' + str(a.FWHM) + ', (' + str(a.chi2Data()[1]) + ')') # + str(a.chi2Data()[1]) + ')')
#                 plt.plot(rangeY, line)
#                 plt.plot(b, a(b))
#                 plt.show()      
         
def find_PSF_flat():
    import pyfits as pf
    import wsm
    from pymodelfit.builtins import VoigtModel
    from pymodelfit.builtins import GaussianModel
    from pymodelfit import get_model_instance
    
    #Variable declaration
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    scienceFile = base_dir + '10nov30044.fits'
    biasFile = base_dir + 'BIAScombined3.fits'
    tramLinefile = base_dir + '10nov30044tlm.fits'
    method = 'gaussian'
    
    #booleans
    plotPoints = 0
    
    #Read, calibrate and create im object
    bias = pf.open(biasFile)
    tlm = pf.open(tramLinefile)
    biasIm = bias[0].data
    hdulist = pf.open(scienceFile)
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']
    im = hdulist[0].data     
    im -= biasIm
        
    #Open output files
    if method=='voigt':
        fileOutSpa = open(base_dir + 'psf_HERMES1_spa_Voigt4.txt','w')
    elif method=='gaussian':
        fileOutSpa = open(base_dir + 'psf_HERMES1_spa_Gauss4.txt','w')
    fileOutSpa.write('x_cent, y_cent, FWHM, A, mu, sigma \n')
    sigTest = 5
#     imageMapY = np.arange(25,97,8) 
#     imageMapY = np.tile(imageMapY,40)
#     offset1 = np.arange(imHeight-25,96)
#     offsetAll = np.repeat(offset1, 10)
#     imageMapY = imageMapY + offsetAll
    
#     print imageMapY
    for j in range(0,imWidth/10,10): #range(areaRange):
        for i in tlm[0].data[:,j]:
        
            area = im[i-sigTest:i+sigTest, j:j+10]
            if (plotPoints==True):    
                plt.subplot(121)
                plt.imshow(area, origin='lower')
                plt.set_cmap(cm.Greys_r) 
    
            #Spatial direction
            line = im[i-sigTest:i+sigTest, j]
            line[np.isnan(line)] = 0
            b = np.arange(int(i-sigTest), int(i+sigTest))
            if method=='voigt':
                a = get_model_instance('Voigt')
                a.parvals = (max(line),1,1,np.average(b))
                a.parvals = VoigtModel.fitData(a, b, line)
            elif method=='gaussian':
                a = get_model_instance('gaussian')
                a.parvals = (max(line),1,np.average(b))
                a.parvals = GaussianModel.fitData(a, b, line)
            c = np.arange(min(b),max(b),0.01)

            #output and plot
            outString = (str(j*10) + ', ' + 
                        str(i) + ', ' + 
                        str(a.FWHM) + ', ' + 
                        str(a.pardict['A']) + ', ' + 
                        str(a.pardict['mu']) + ', ' + 
                        str(a.pardict['sig']) + ', ' + 
                        '\n' )
            if (a.chi2Data()[1]<100): fileOutSpa.write(outString)
            print 'Spatial: ' + str(a.FWHM), a.chi2Data()
            if (plotPoints==True):
                plt.subplot(122)
                plt.title('Spatial ' + str(a.FWHM) + ', (' + str(a.chi2Data()[1]) + ')') # + str(a.chi2Data()[1]) + ')')
                plt.plot(b, line)
                plt.plot(c, a(c))
                plt.show()      
             
def read_psf_data_spatial(method, cam):
#     from scipy.interpolate import interp1d
#     from scipy.interpolate import interp2d
#     import pyfits as pf
    
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    dataFile = base_dir + 'psf_HERMES1_spa_' + method + str(cam) +'.txt'
    
    a = np.loadtxt(dataFile,skiprows=1, delimiter = ',' , usecols = [0,1,2,3,4,5]).transpose()
    a = a.transpose()[a[5]<7].transpose()
    a = a.transpose()[a[5]>1].transpose()

    #scatter plot
#     plt.scatter(a[0],a[1],s = a[2])                                                                                      
#     plt.show()
    
    #binned plot
    rows = 4112
    cols = 4146
    bins = 4
    #bin rows
    for i in range(bins):
        binLimits = (rows/bins*i,rows/bins*(i+1))
        mapBin = ((a[1]>binLimits[0]) & (a[1]<=binLimits[1]))
        b = a[0,mapBin], a[5,mapBin]
        c = np.bincount(b[0].astype(int), weights=b[1])
        d = np.bincount(b[0].astype(int))
        px = np.arange(len(d))
        px = px[d>0]
        c = c[d>0]
        d = d[d>0]
        c = c/d
        plt.plot(px,c, label= ' Range (y) = ' + str(binLimits))
        plt.title('Spatial PSF - y-binned - ' + method)
        plt.xlabel('X-Pixels')
        plt.ylabel('Std. Dev. (px)')
        plt.legend()
    plt.savefig(base_dir + 'spa_' + str(cam) + '_y_' + method +'.png')
    plt.show()

    #bin cols
    for i in range(bins):
        binLimits = (cols/bins*i,cols/bins*(i+1))
        mapBin = ((a[0]>binLimits[0]) & (a[0]<=binLimits[1]))
        b = a[1,mapBin], a[5,mapBin]
        c = np.bincount(b[0].astype(int), weights=b[1])
        d = np.bincount(b[0].astype(int))
        px = np.arange(len(d))
        px = px[d>0]
        c = c[d>0]
        d = d[d>0]
        c = c/d
        plt.plot(px,c, label= ' Range (x) = ' + str(binLimits))
        plt.title('Spatial PSF - x-binned - ' + method)
        plt.xlabel('Y-Pixels')
        plt.ylabel('Std. Dev. (px)')
        plt.legend()
    plt.savefig(base_dir + 'spa_' + str(cam) + '_x_' + method +'.png')
    plt.show()


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
       
def read_psf_data_spectral(method, cam):
    from scipy.interpolate import interp1d
    from scipy.interpolate import interp2d
    import pyfits as pf
    
    base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    dataFile = base_dir + 'psf_HERMES1_spe_' + method + str(cam) +'.txt'
    
    a = np.loadtxt(dataFile,skiprows=1, delimiter = ',' , usecols = [0,1,2,3,4,5]).transpose()
    a = a.transpose()[a[5]<7].transpose()
    a = a.transpose()[a[5]>1].transpose()

    #scatter plot
#     plt.scatter(a[0],a[1],s = a[2]*10)                                                                                      
#     plt.show()
    
    #binned plot
    rows = 4112
    cols = 4146
    bins = 4
    #bin rows
    for i in range(bins):
        binLimits = (rows/bins*i,rows/bins*(i+1))
        mapBin = ((a[1]>binLimits[0]) & (a[1]<=binLimits[1]))
        b = a[0,mapBin], a[5,mapBin]
        c = np.bincount(b[0].astype(int), weights=b[1])
        d = np.bincount(b[0].astype(int))
        px = np.arange(len(d))
        px = px[d>0]
        c = c[d>0]
        d = d[d>0]
        c = c/d
        plt.plot(px, c, label= ' Range (y) = ' + str(binLimits))
        plt.title('Spectral PSF - y-binned - ' + method)
        plt.xlabel('X-Pixels')
        plt.ylabel('Std. Dev. (px)')
        plt.legend()
    plt.savefig(base_dir + 'spec_' + str(cam) + '_y_' + method +'.png')
    plt.show()

    #bin cols
    for i in range(bins):
        binLimits = (cols/bins*i,cols/bins*(i+1))
        mapBin = ((a[0]>binLimits[0]) & (a[0]<=binLimits[1]))
        b = a[1,mapBin], a[5,mapBin]
        c = np.bincount(b[0].astype(int), weights=b[1])
        d = np.bincount(b[0].astype(int))
        px = np.arange(len(d))
        px = px[d>0]
        c = c[d>0]
        d = d[d>0]
        c = c/d
        plt.plot(px, c, label= ' Range (x) = ' + str(binLimits))
        plt.title('Spectral PSF - x-binned - ' + method)
        plt.xlabel('Y-Pixels')
        plt.ylabel('Std. Dev. (px)')
        plt.legend()
    plt.savefig(base_dir + 'spec_' + str(cam) + '_x_' + method +'.png')
    plt.show()


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
    