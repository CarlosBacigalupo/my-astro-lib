import math
import numpy as np
import string

def gd2jd(date, TZ = 11):
    """gd2jd.py converts a UT Gregorian date to Julian date.

    Usage: gd2jd.py (2009, 02, 25, 01, 59, 59)

    To get the current Julian date:
        import time
        gd2jd(time.gmtime())

    Hours, minutes and/or seconds can be omitted -- if so, they are
    assumed to be zero.

    Year and month are converted to type INT, but all others can be
    type FLOAT (standard practice would suggest only the final element
    of the date should be float)
"""
#     print date

    date = list(date)

    if len(date)<3:
        print "You must enter a date of the form (2009, 02, 25)!"
        return -1
    elif len(date)==3:
        for ii in range(3): date.append(0)
    elif len(date)==4:
        for ii in range(2): date.append(0)
    elif len(date)==5:
        date.append(0)

    yyyy = int(date[0])
    mm = int(date[1])
    dd = float(date[2])
    hh = float(date[3]-TZ)
    min = float(date[4])
    sec = float(date[5])
#     print yyyy,mm,dd,hh,min,sec

    UT=hh+min/60+sec/3600

#     print "UT="+`UT`

    total_seconds=hh*3600+min*60+sec
    fracday=total_seconds/86400

#     print "Fractional day: %f" % fracday
# print dd,mm,yyyy, hh,min,sec, UT

    if (100*yyyy+mm-190002.5)>0:
        sig=1
    else:
        sig=-1

    JD = 367*yyyy - int(7*(yyyy+int((mm+9)/12))/4) + int(275*mm/9) + dd + 1721013.5 + UT/24 - 0.5*sig +0.5

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

#     print "\n"+months[mm-1]+" %i, %i, %i:%i:%i UT = JD %f" % (dd, yyyy, hh, min, sec, JD),

# Now calculate the fractional year. Do we have a leap year?
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    daysum=0
    for y in range(mm-1):
        daysum=daysum+days[y]
    daysum=daysum+dd-1+UT/24

    if days[1]==29:
        fracyear=yyyy+daysum/366
    else:
        fracyear=yyyy+daysum/365
#     print " = " + `fracyear`+"\n"


    return JD


def jd2gd(jd):

    """Task to convert a list of julian dates to gregorian dates
    description at http://mathforum.org/library/drmath/view/51907.html
    Original algorithm in Jean Meeus, "Astronomical Formulae for
    Calculators"

    2009-02-15 13:36 IJC: Converted to importable, callable function
    """
    
    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5:
        mm=E-1

    if E>13.5:
        mm=E-13

    if mm>2.5:
        yyyy=C-4716

    if mm<2.5:
        yyyy=C-4715

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    min=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    hh = 24.0*(dd % 1.0)
    min = 60.0*(hh % 1.0)
    sec = 60.0*(min % 1.0)

    dd =  dd-(dd%1.0)
    hh =  hh-(hh%1.0)
    min =  min-(min%1.0)


    print str(jd)+" = "+str(months[mm-1])+ ',' + str(dd) +',' +str(yyyy)
    print string.zfill(h,2)+":"+string.zfill(min,2)+":"+string.zfill(sec,2)+" UTC"

    print (yyyy, mm, dd, hh, min, sec)

    return (yyyy, mm, dd, hh, min, sec)


def sex2dec(hd, min, sec):
    """ 
    Converts a Sexagesimal number to a Decimal number.
    
    Parameters
    ----------
    hd : int
        hours or degrees
    m : int
        minutes or arcminutes
    s : float
        seconds or arcseconds
    
    Returns
    -------
    hd : float
        A decimal number
    
    """
    return float(hd) + min/60.0 + sec/3600.0
    
def dec2sex(deci):
    """ 
    Converts a Decimal number (in hours or degrees) to Sexagesimal.
    
    Parameters
    ----------
    deci : float
        A decimal number to be converted to Sexagismal.
    
    Returns
    -------
    hd : int
        hours or degrees
    m : int
        minutes or arcminutes
    s : float
        seconds or arcseconds
    
    """
    #original func, now array version
#     (hfrac, hd) = math.modf(deci)
#     (min_frac, m) = math.modf(abs(hfrac) * 60)
#     s = min_frac * 60.
    
    hd = np.floor(np.abs(deci)) * deci/np.abs(deci)
    hfrac = np.abs(deci - hd)
    hfrac *= 60
    
    m = np.floor(hfrac)
    mfrac = hfrac - m
    
    s = mfrac * 60
    
    return (hd, m, s)

def write_fld(main_df, sky_df, fiducials_df, fld_file, cluster_centre_sex):
    
    f = open(fld_file ,'w')
    
    out_string = 'LABEL 47 Tuc - Commissioning field \n'
    f.write(out_string) 
    out_string = 'UTDATE 2013 10 10\n'
    f.write(out_string) 
    out_string = 'CENTRE ' + cluster_centre_sex[0] + ' ' + cluster_centre_sex[1] +  '\n'     
    f.write(out_string) 
    out_string = 'EQUINOX J2000.0\n'
    f.write(out_string) 
    out_string = '*** End of Header\n'
    f.write(out_string) 
    out_string = '*'
    f.write(out_string)
    out_string = '**Program Stars\n'
    f.write(out_string) 
    
    
    ###Main Targets
    stars = 392    
    if main_df.shape[0] < 392:
        stars = main_df.shape[0]
        
    for i in range(stars):
        out_string = ('{:8s}'.format(str(int(main_df['No'][i]))) + '  ' + 
                      '{:02.0f}'.format(int(main_df['RA1'][i])) + '  ' +
                      str(int(main_df['RA2'][i])) + ' ' + 
                      str('{:05.2f}'.format(main_df['RA3'][i])) + '  ' +
                      '{:02.0f}'.format(int(main_df['Dec1'][i])) + ' ' +
                      str(int(main_df['Dec2'][i])) + '  ' + 
                      str('{:05.2f}'.format(main_df['Dec3'][i])) + '  ' +
                      'P 1' + ' ' +
                      str('{:05.2f}'.format(main_df['V'][i])) + '  ' + 
                      '0' + '  ' +
                      'star' + '  ' +
                      '\n')             
        f.write(out_string) 
    
    out_string = '*** End of program stars\n'
    f.write(out_string) 
    out_string = '*'
    f.write(out_string)
    out_string = '**Sky Fibres\n'
    f.write(out_string) 
    
    ####Sky fibres
    sky = 392    
    if sky_df.shape[0] < 392:
        sky = sky_df.shape[0]
        
    for i in range(sky):
        out_string = ('{:8s}'.format(str(int(sky_df['No'][i]))) + '  ' + 
                      '{:02.0f}'.format(int(sky_df['RA1'][i])) + '  ' +
                      str(int(sky_df['RA2'][i])) + ' ' + 
                      str('{:05.2f}'.format(sky_df['RA3'][i])) + '  ' +
                      '{:02.0f}'.format(int(sky_df['Dec1'][i])) + ' ' +
                      str(int(sky_df['Dec2'][i])) + '  ' + 
                      str('{:05.2f}'.format(sky_df['Dec3'][i])) + '  ' +
                      'S 8' + ' ' +
                      str('{:05.2f}'.format(sky_df['V'][i])) + '  ' + 
                      '0' + '  ' +
                      'sky' + '  ' +
                      '\n')             
        f.write(out_string) 
        
        
        
    ###Fiducials 
    out_string = '*** End of Sky fibres\n'
    f.write(out_string) 
    out_string = '*'
    f.write(out_string)
    out_string = '**Guiding Stars\n'
    f.write(out_string) 
    
    fiducials = 20    
    if fiducials_df.shape[0] < 20:
        fiducials = fiducials_df.shape[0]
        
    for i in range(fiducials):
        out_string = ('{:8s}'.format(str(int(fiducials_df['No'][i]))) + '  ' + 
                      '{:02.0f}'.format(int(fiducials_df['RA1'][i])) + '  ' +
                      str(int(fiducials_df['RA2'][i])) + ' ' + 
                      str('{:05.2f}'.format(fiducials_df['RA3'][i])) + '  ' +
                      '{:02.0f}'.format(int(fiducials_df['Dec1'][i])) + ' ' +
                      str(int(fiducials_df['Dec2'][i])) + '  ' + 
                      str('{:05.2f}'.format(fiducials_df['Dec3'][i])) + '  ' +
                      'F 9' + ' ' +
                      str('{:05.2f}'.format(fiducials_df['V'][i])) + '  ' + 
                      '0' + '  ' +
                      'fiducial star' + '  ' +
                      '\n')             
        f.write(out_string) 
    
    
    
    
    
    
    
    
    f.close()
    
