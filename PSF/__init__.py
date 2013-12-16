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
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    sexParamFile = WORKING_DIR + 'HERMES.sex'
    outputFileName = WORKING_DIR + 'out.txt'
    scienceFile = WORKING_DIR + '10nov40045.fits'
    biasFile = WORKING_DIR + 'BIAScombined4.fits'
    tramLinefile = WORKING_DIR + '10nov10044tlm.fits'
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
#     fileOutSpa = open(WORKING_DIR + 'psf_HERMES1_spa_Voigt_arc.txt','w')
    if method=='voigt':
        fileOutSpe = open(WORKING_DIR + 'psf_HERMES1_spe_Voigt4.txt','w')
    elif method=='gaussian':
        fileOutSpe = open(WORKING_DIR + 'psf_HERMES1_spe_Gauss4.txt','w')
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
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    scienceFile = WORKING_DIR + '10nov30044.fits'
    biasFile = WORKING_DIR + 'BIAScombined3.fits'
    tramLinefile = WORKING_DIR + '10nov30044tlm.fits'
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
    im = im-biasIm
        
    #Open output files
    if method=='voigt':
        fileOutSpa = open(WORKING_DIR + 'psf_HERMES1_spa_Voigt4.txt','w')
    elif method=='gaussian':
        fileOutSpa = open(WORKING_DIR + 'psf_HERMES1_spa_Gauss4.txt','w')
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
    
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    dataFile = WORKING_DIR + 'psf_HERMES1_spa_' + method + str(cam) +'.txt'
    
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
    plt.savefig(WORKING_DIR + 'spa_' + str(cam) + '_y_' + method +'.png')
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
    plt.savefig(WORKING_DIR + 'spa_' + str(cam) + '_x_' + method +'.png')
    plt.show()


#     grid = np.zeros((4112,4146))
#     a[2] = a[2]/max(a[2])*65535
#     pointSize = 10
#     for i in range(len(a[0])):
#         grid[int(a[1][i]-pointSize):int(a[1][i]+pointSize),int(a[0][i]-pointSize):int(a[0][i]+pointSize)] = a[2][i]
#     plt.imshow(grid, origin='lower')
#     plt.set_cmap(cm.Greys_r)                                                                                    
#     plt.show()
#     pf.writeto(WORKING_DIR + 'a.fits', grid)


    #interpolate function
#     f = interp2d(a[0][range(0,len(a[0]),1)],a[1][range(0,len(a[0]),1)],a[5][range(0,len(a[0]),1)])#, kind='cubic')
# #     grid = np.zeros((4112,4146))
# #     grid = f(a[0].astype('int'), a[1].astype('int'))
#     grid = f(range(4112),range(4146))
# #     grid = f(range(100),range(100))
#     plt.imshow(grid, origin='lower')
#     plt.set_cmap(cm.Greys_r)                                                                                    
#     plt.show()
#     pf.writeto(WORKING_DIR + 'a.fits',grid, clobber= True)
       
def read_psf_data_spectral(method, cam):
    from scipy.interpolate import interp1d
    from scipy.interpolate import interp2d
    import pyfits as pf
    
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
    dataFile = WORKING_DIR + 'psf_HERMES1_spe_' + method + str(cam) +'.txt'
    
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
    plt.savefig(WORKING_DIR + 'spec_' + str(cam) + '_y_' + method +'.png')
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
    plt.savefig(WORKING_DIR + 'spec_' + str(cam) + '_x_' + method +'.png')
    plt.show()


#     grid = np.zeros((4112,4146))
#     a[2] = a[2]/max(a[2])*65535
#     pointSize = 10
#     for i in range(len(a[0])):
#         grid[int(a[1][i]-pointSize):int(a[1][i]+pointSize),int(a[0][i]-pointSize):int(a[0][i]+pointSize)] = a[2][i]
#     plt.imshow(grid, origin='lower')
#     plt.set_cmap(cm.Greys_r)                                                                                    
#     plt.show()
#     pf.writeto(WORKING_DIR + 'a.fits', grid)


    #interpolate function
#     f = interp2d(a[0][range(0,len(a[0]),1)],a[1][range(0,len(a[0]),1)],a[5][range(0,len(a[0]),1)])#, kind='cubic')
# #     grid = np.zeros((4112,4146))
# #     grid = f(a[0].astype('int'), a[1].astype('int'))
#     grid = f(range(4112),range(4146))
# #     grid = f(range(100),range(100))
#     plt.imshow(grid, origin='lower')
#     plt.set_cmap(cm.Greys_r)                                                                                    
#     plt.show()
#     pf.writeto(WORKING_DIR + 'a.fits',grid, clobber= True)
    