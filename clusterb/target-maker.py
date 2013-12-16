from constants import *
import numpy as np

global a

def reduce_targetList(targetList, constraints = {}):
    
    for item in constraints.keys():
        const = constraints[item]
        idx = np.where(np.array(header)==item)[0][0]
        if len(const) == 1:
            targetList.transpose()[targetList[idx] == const].transpose()  
        else:
            targetList = targetList.transpose()[targetList[idx] >= const[0]].transpose()
            if len(targetList)>0:
                targetList = targetList.transpose()[targetList[idx] <= const[1]].transpose() 

    return targetList


    
def add_lumin_to_targetList(targetList):

    parallax = targetList[parallaxCol]
    vMag = targetList[vMagCol]

    d = 1/parallax*1000 #mas/1000
    dMod = 5 * np.log10(d/10)
    aMag = vMag- dMod
        
    targetList[dCol] =  d
    targetList[luminCol] =  np.exp((SOLAR_AMAG - aMag)/2.5)

    return targetList

def add_mass_to_targetList(targetList):

    lumin = targetList[luminCol]   
    mass = lumin**(1/4)
        
    targetList[massCol] =  mass

    return targetList

def add_PMax_to_tagetList(targetList):
            
    mass = targetList[massCol]
    temp = targetList[tempCol]
    lumin = targetList[luminCol]
    
    nuMax = ((mass * ((temp) ** 3.5)) / lumin) * (NU_MAX_SUN) #in microHz
    
    PMax = 1/nuMax*1e6
    
    targetList[nuMaxCol] = nuMax
    targetList[PMaxCol] = PMax
    targetList[PMaxHrsCol] = PMax/3600
    targetList[PMaxDaysCol] = PMax/3600/24
    
    return targetList

def add_temp_to_targetList(targetList):  
#    V = targetList[VCol]
#    B = targetList[BCol]    
    CI = targetList[CICol]
    
    targetList[tempCol] = 8540/(CI + 0.865)/TEMP_EFF_SUN

    return targetList

def load_initial_data(catalogueFileName):
    a = np.genfromtxt(catalogueFileName, unpack=True, delimiter = '|', usecols = (0,3,4,5,6,7,8), missing_values = 0 )

    return a

def add_columns(targetList, cols = 1):
     a = np.vstack((targetList, np.zeros((cols, len(targetList[0])))))

     return  a


def create_targetList():
#    targetList = load_initial_data('hip_short.txt.tar.gz')
    targetList = load_initial_data('hip_short.txt')
#    targetList = load_initial_data('initial_catalogue.txt')
    
    targetList = add_columns(targetList, cols = len(header_add))
    
    targetList = add_temp_to_targetList(targetList)
    
    targetList = add_lumin_to_targetList(targetList)
    
    targetList = add_mass_to_targetList(targetList)
    
    targetList = add_PMax_to_tagetList(targetList)
    

    return targetList


a = create_targetList()
b = reduce_targetList(a,{'PMaxDays':[1,100]})
print b
