# import math
import numpy as np
# import scipy.constants as const
# from scipy import signal
import pylab as plt
# import matplotlib.cm as cm

# import RVSimulator as RVS
# import TableBrowser as TB
# import clusterb as CB
# import PSF
import toolbox
# import pandas as pd
import HERMES


# x=np.arange(0,100)
# A=np.array([1,1,1])
# mu =np.array([25,50,75])
# sigma = np.array([2,1,0.5])
# plt.plot(toolbox.gaussian(x,mu,sigma,A))
# plt.show()


a = HERMES.PSF()
a.cur_camera = 4
a.imWidth =4100
columns = range(0, a.imWidth,(a.imWidth-1)/1)
print 'Columns ',columns
a.read_full_image_spatial('gaussian', columns)

# a.read_psf_data('gaussian','ref4.fits',4)



# #ccd1 single fibre per bundle fit
# a.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/single/1/'
# a.scienceFile = a.base_dir + '12nov10063.fits'
# a.biasFile = a.base_dir + '11novBIAScomb_p0.fits'
# a.flatFile = a.base_dir + '12nov10063.fits'
# a.tramLinefile = a.base_dir + '12nov10063tlm.fits'
# a.sexParamFile = a.base_dir + 'HERMES.sex'
# a.outputFileName = a.base_dir + 'out.txt'
# a.nFibres = 10
# 
# a.profile = 'gaussian'
# a.base_dir = '/Users/Carlos/Documents/HERMES/reductions/resolution_gayandhi/'
# a.out_dir = a.base_dir + 'output'
# a.flatFile = a.base_dir + '10nov10044.fits'
# a.tramLinefile = a.base_dir + '10nov10044tlm.fits'
# a.deadFibres = [41,66, 91, 141,191,241, 294, 307, 341, 391]
# 
# a.open_files()
# a.flatIm_b = a.bias_subtract_from_overscan(a.flatIm, range(-45,-5)) #bias subtract using the last 40 columns(overscan region)
# # a.bias_subtract()
# a.nFibres = 400
# # print a.fit_pShort(a.flatIm[:, 10])
# 
# column = 1000
# flatCurve = a.flatIm_b[:,column]
# a.write_p400_A_mu(flatCurve)
# a.write_mu_from_tlm(a.p400_A_mu,column)
# 
# p = a.fit_p400_A_mu(flatCurve)
# a.A = p[0][:a.nFibres]
# a.write_p400_sigma_gamma(flatCurve)
# p = a.fit_p400_sigma_gamma(flatCurve)


# 
# print a.fit_10f(a.scienceIm_b[:,10])
 
 
 
 
 
# dr2df = HERMES.dr2df()    
# batch_dir = '/Users/Carlos/Documents/HERMES/reductions/batch_reductions/'
# dr2df.base_dir = '/Users/Carlos/Documents/HERMES/data/140109/'
# dr2df.dr_dir = '/Users/Carlos/Documents/workspace/2dfdr/5.49/src_code/2dfdr-5.49/bin/'
# dr2df.red_dir = '/Users/Carlos/Documents/HERMES/reductions/m67_complete/'
# 
# ix_array = [range(48,53) ] #2 sets of files for m67 short
# date_array = ['09jan'] #prefix for the filenames
# arc_array=[4] #index of the arc
# flat_array=[5] #index of the flat
#  
# for i in range(len(date_array)):
#     dr2df.target_dir = batch_dir + str(i) +'/'
#     dr2df.file_ix = ix_array[i]
#     dr2df.date = date_array[i] 
#     dr2df.arc = arc_array[i] 
#     dr2df.flat = flat_array[i] 
#     dr2df.auto_reduce() 

# RV = HERMES.RV()
# RV.files = ['17dec10039red.fits', '17dec10040red.fits', '17dec10041red.fits']
# RV.base_dir = '/Users/Carlos/Documents/HERMES/reductions/m67_complete/'
# RV.read_fits_files()
# RV.calculate_RV_shifts()
     
def aaaa():
    
    
    
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/ngc2477/'
    FLD_DIR = WORKING_DIR + 'fld/'
    
    main = np.loadtxt(WORKING_DIR + 'giants.txt', unpack=True)
#         prob = np.loadtxt(self.base_dir + 'prob.mu', unpack=True, skiprows = 2)
#         UBV = np.loadtxt(self.base_dir + 'ubv.ccd', unpack=True, skiprows = 2, usecols = [0,1,2,3])
    
    #These lines convert the np arrays into panda tables.
    #Select the format before cross reference here
    main_df = pd.DataFrame({ 'target' : main[0],
                            'RA_h' : main[1],
                            'RA_min' : main[2],
                            'RA_sec' : main[3],
                            'Dec_deg' : main[4],
                            'Dec_min' : main[5],
                            'Dec_sec' : main[6],
                            'P' : np.ones(len(main[6]))*9,
                            'mag' : main[7]})

#         prob_df = pd.DataFrame({ 'No' : prob[0],
#                                 'Prob' : prob[2]})
#     
#         UBV_df = pd.DataFrame({ 'No' : UBV[0],
#                                 'V' : UBV[2],
#                                 'BV' : UBV[3]})
#     
        #Creates the main merged table
#         main_df = main_df.merge(prob_df, on='No')
#         main_df = main_df.merge(UBV_df, on='No')
  
    sky_df = main_df[0:0]
    fiducials_df = main_df[0:0]
    toolbox.write_fld(main_df, sky_df, fiducials_df, FLD_DIR + 'suplement.fld', ('7 52 09.8','-38 31 48.00'),title = 'NGC2477 Short Field')



def ngc2477_report():

    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/ngc2477/'
    fldFile = 'test.fld' 
    booCreateFld = False
    booCone = True
    
    if booCone==True:
        cone = CB.vizier()
        cone.base_dir = WORKING_DIR
        cone.fileName = 'vizier_cone_full.tsv'
        cone.cluster_pm = (-0.83, 1.89) #bonatto+ 2011
        cone.cluster_centre = (toolbox.sex2dec(07,52,09.8)*15,toolbox.sex2dec(-38,31,48))
        cone.tidal_rad = toolbox.sex2dec(2,00,00)
        cone.pp_cutoff = 1
        cone.mag_max = 30
        cone.mag_min = -10
        cone.load_data()
        cone.add_Vmag_galah()
        cone.add_BVmag_ppmxl()
        cone.reduce_list()
    #     cone.plot_pm()
    
        cone.plot_HR_pos_pm()    

def m67_report():
    

    #initial variables
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/m67/'
    fldFile = 'test.fld'
    cluster_centre = ('08:51:18.0', '+11:48:00')
    booCreateFld = False
    booStetson = True
    booPlotStetson = True
    booCone = False
    booYadav = False
    
    if booStetson == True:
        #stetson's positions data frame
        tableData = TB.read_external_file(WORKING_DIR + 'Positions.txt')  
        d = range(len(tableData))
        d[0] = 'RA_dec'
        d[1] = 'Dec_dec'
        d[2] = 'RA_h'
        d[3] = 'RA_min'
        d[4] = 'RA_sec'
        d[5] = 'Dec_deg'
        d[6] = 'Dec_min'
        d[7] = 'Dec_sec'
        d[12] = 'target'
        stetson_df = TB.build_DataFrame(tableData, d)
        
        #adds stetson's photometry data qand B-V column
        tableData = TB.read_external_file(WORKING_DIR + 'Photometry.txt')  
        header = tableData.transpose()[0] 
        d = np.hstack((['target'], header))
        photometry_df = TB.build_DataFrame(tableData, d)
        stetson_df = stetson_df.merge(photometry_df, on='target')
        stetson_df['B'] = stetson_df['B'].astype(float)
        stetson_df['V'] = stetson_df['V'].astype(float)
        stetson_df['B'][stetson_df['B']==99.999] = np.nan
        stetson_df['BV'] = stetson_df['B'] - stetson_df['V']
        
        #bright stars
        dataFileNames = [WORKING_DIR + '17dec10044red.fits', WORKING_DIR + '17dec10045red.fits', WORKING_DIR + '17dec10046red.fits'] 
        M67_bright = TB.field()
        M67_bright.dataFileNames = dataFileNames 
        M67_bright.load_data()
        bright_df = stetson_df.merge(M67_bright.FibreTable[0].df ,on = 'target')
           
        #12V14 stars
        dataFileNames = [WORKING_DIR + '17dec10039red.fits', WORKING_DIR + '17dec10040red.fits', WORKING_DIR + '17dec10041red.fits'] 
        M67_12V14 = TB.field()
        M67_12V14.dataFileNames = dataFileNames 
        M67_12V14.load_data()
        V12V14_df = stetson_df.merge(M67_12V14.FibreTable[0].df ,on = 'target')
       
        #V14V16 (using data from the lis file temporarily)
        fileName = WORKING_DIR + 'M67_14V16_p1.lis'
        M67_14V16 = TB.FibreTable(fileName)
        V14V16_df = stetson_df.merge(M67_14V16.df ,on = 'target')
        
        if booPlotStetson==True:    
            plt.scatter(stetson_df['BV'], stetson_df['V'], color = 'grey', label = 'Complete Stetson catalogue' + str(stetson_df.shape))
            plt.scatter(bright_df['BV'], bright_df['V'], color = 'red', label = 'Bright' + str(bright_df.shape))
            plt.scatter(V12V14_df['BV'], V12V14_df['V'], color = 'green', label = '12V14' + str(V12V14_df.shape))
            plt.scatter(V14V16_df['BV'], V14V16_df['V'], color = 'blue', label = '14V16' + str(V14V16_df.shape))
        
            plt.xlabel('B - V')
            plt.ylabel('V')
            plt.gca().invert_yaxis()
            plt.legend()
            plt.show()

    if booYadav==True:
        yadav = CB.vizier()
        yadav.base_dir = WORKING_DIR
        yadav.fileName = 'vizier_yadav_full.tsv'
        yadav.load_data()
        yadav.add_BVmag_yadav()
        yadav.reduce_list()
    #     yadav.plot_pm()
        yadav.plot_HR_pos_pm()
        
    if booCone==True:
        cone = CB.vizier()
        cone.base_dir = WORKING_DIR
        cone.fileName = 'vizier_cone_full.tsv'
        cone.load_data()
        cone.add_Vmag_galah()
        cone.add_BVmag_ppmxl()
        cone.reduce_list()
    #     cone.plot_pm()
    
        cone.plot_HR_pos_pm()
        
        plt.scatter(cone.reduced['BVmag'], cone.reduced['Vmag'], color ='grey')
        plt.scatter(yadav.reduced['BVmag'], yadav.reduced['Vmag'], color = 'blue')
        plt.gca().invert_yaxis()
        plt.show()
        
    if booCreateFld==True:
        stetson_df['mag'] = stetson_df['V']
        sky_df = stetson_df
        fiducials_df = stetson_df
        toolbox.write_fld(stetson_df, sky_df, fiducials_df, WORKING_DIR + 'test.fld', cluster_centre)
    
    
def m67_analysis():
    
    WORKING_DIR = '/Users/Carlos/Documents/HERMES/reductions/ngc2477/'
    
#     a = CB.webda()
#     a.load_tables()
#     plt.scatter(a.df['BV'], a.df['V'])
    
    
    cone = CB.vizier()
    cone.base_dir = WORKING_DIR
    cone.fileName = 'vizier_cone_UCAC4_full.tsv'
    cone.load_data()
    cone.reduce_list()
    

#### 
#     fileName = WORKING_DIR + 'M67_bright_p1.lis'
#     M67_bright = TB.FibreTable(fileName)
# 
#     fileName = WORKING_DIR + 'M67_12V14_p0.lis'
#     M67_12V14 = TB.FibreTable(fileName)
# 
    fileName = WORKING_DIR + 'M67_14V16_p1.lis'
    M67_14V16 = TB.FibreTable(fileName)
# 
#     fileName = WORKING_DIR + 'M67_bright.fld'
#     M67_bright_fld = TB.FibreTable(fileName)
# 
#     fileName = WORKING_DIR + 'M67_12V14.fld'
#     M67_12V14_fld = TB.FibreTable(fileName)
# 
#     fileName = WORKING_DIR + 'M67_14V16.fld'
#     M67_14V16_fld = TB.FibreTable(fileName)
#
#   
#     plt.subplot(121)
#     plt.scatter(M67_bright.RA_dec[M67_bright.type=='P'], M67_bright.Dec_dec[M67_bright.type=='P'], color = 'red', label = 'Bright Stars')   
#     plt.scatter(M67_12V14.RA_dec[M67_12V14.type=='P'], M67_12V14.Dec_dec[M67_12V14.type=='P'], color = 'green', label = '12V14 stars')   
#     plt.scatter(M67_14V16.RA_dec[M67_14V16.type=='P'], M67_14V16.Dec_dec[M67_14V16.type=='P'], color = 'blue', label = '14V16 stars')   
#     plt.legend()
# 
#     plt.subplot(122)
#     plt.scatter(M67_bright_fld.RA_dec, M67_bright_fld.Dec_dec, color = 'black')   
#     plt.scatter(M67_12V14_fld.RA_dec, M67_12V14_fld.Dec_dec, color = 'black')   
#     plt.scatter(M67_14V16_fld.RA_dec, M67_14V16_fld.Dec_dec, color = 'black')   
#     plt.scatter(M67_bright.RA_dec, M67_bright.Dec_dec, color = 'red', label = 'Bright Stars')   
#     plt.scatter(M67_12V14.RA_dec, M67_12V14.Dec_dec, color = 'green', label = '12V14 stars')   
#     plt.scatter(M67_14V16.RA_dec, M67_14V16.Dec_dec, color = 'blue', label = '14V16 stars')   
#     plt.legend()
# 
#     plt.show()
          
    #positions data frame
    tableData = TB.read_external_file(WORKING_DIR + 'Positions.txt')  
    d = range(len(tableData))
    d[0] = 'RA_dec'
    d[1] = 'Dec_dec'
    d[2] = 'RA_h'
    d[3] = 'RA_min'
    d[4] = 'RA_s'
    d[5] = 'Dec_deg'
    d[6] = 'Dec_min'
    d[7] = 'Dec_s'
    d[12] = 'target'
    stetson_df = TB.build_DataFrame(tableData, d)
    
    #photometry data frame
    tableData = TB.read_external_file(WORKING_DIR + 'Photometry.txt')  
    header = tableData.transpose()[0] 
    d = np.hstack((['target'], header))
    photometry_df = TB.build_DataFrame(tableData, d)
    
    stetson_df = stetson_df.merge(photometry_df, on='target')
    
    mag_B = stetson_df['B'].astype(float)
    mag_V = stetson_df['V'].astype(float)
    mag_BV = mag_B - mag_V
    plt.scatter(mag_BV[mag_BV<10], mag_V[mag_BV<10], color = 'grey', label = 'Complete Stetson catalogue' + stetson_df.shape())

    #bright stars
    dataFileNames = [WORKING_DIR + '17dec10044red.fits', WORKING_DIR + '17dec10045red.fits', WORKING_DIR + '17dec10046red.fits'] 
    M67_bright = TB.field()
    M67_bright.dataFileNames = dataFileNames 
    M67_bright.load_data()
    bright_df = stetson_df.merge(M67_bright.FibreTable[0].df ,on = 'target')
    mag_B = bright_df['B'].astype(float)
    mag_V = bright_df['V'].astype(float)
    mag_BV = mag_B - mag_V
    plt.scatter(mag_BV, mag_V, color = 'yellow', label = 'Bright' + bright_df.shape())
    
    
    #12V14 stars
    dataFileNames = [WORKING_DIR + '17dec10039red.fits', WORKING_DIR + '17dec10040red.fits', WORKING_DIR + '17dec10041red.fits'] 
    M67_12V14 = TB.field()
    M67_12V14.dataFileNames = dataFileNames 
    M67_12V14.load_data()
    V12V14_df = stetson_df.merge(M67_12V14.FibreTable[0].df ,on = 'target')
    mag_B = V12V14_df['B'].astype(float)
    mag_V = V12V14_df['V'].astype(float)
    mag_BV = mag_B - mag_V
    plt.scatter(mag_BV, mag_V, color = 'green', label = '12V14' + V12V14_df.shape())
   
    
    #V14V16 (using data from the lis file temporarily)
    V14V16_df = stetson_df.merge(M67_14V16.df ,on = 'target')
    mag_B = V14V16_df['B'].astype(float)
    mag_V = V14V16_df['V'].astype(float)
    mag_BV = mag_B - mag_V
    plt.scatter(mag_BV, mag_V, color = 'blue', label = '14V16' + V14V16.shape())

    plt.scatter(cone.reduced['BVmag'], cone.reduced['Vmag'], color = 'black', label = 'Reduced cone search')

    plt.xlabel('B - V')
    plt.ylabel('V')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    
    

    
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
    
    #changed by fibretable class - change!
    tableData = RVS.load_2dfdr_fibre_table(fileName) 
    
    df = TB.build_DataFrame(d, tableData)

    return df




# find_PSF_arc()
# find_PSF_flat()
#read_psf_data_spatial()
# read_psf_data_spectral()

# do_all_plots()

# psf_HERMES()
# m67_analysis()
# m67_report()
# ngc2477_report()
# aaaa()  #quick fix to my 2477 data
# run_HERMES_rv()
