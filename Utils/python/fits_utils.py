# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 10:27:00 2016

@author: matvii
"""
from math import *
from astropy.io import fits
from utils import *
def read_fits_date(fitsfile):
    fi=fits.open(fitsfile)
    date=0
    head=fi[0]
    try:
        date=head.header['MJD-OBS']
    except:
        date=float('nan')
    fi.close()
    return date

def calculate_rotated_frame(fitsfile,E):
    fi=fits.open(fitsfile)
    NOCD=0
    head=fi[0]
    try:
        CD11=float(head.header['CD1_1'])
    except:
        NOCD=1
    try:
        CD12=float(head.header['CD1_2'])
    except:
        NOCD=1
    try:
        CD21=float(head.header['CD2_1'])
    except:
        NOCD=1
    try:
        CD22=float(head.header['CD2_2'])
    except:
        NOCD=1
    fi.close()
    if NOCD==1:
        upr=np.array([0,0.3977,0.9175])
        angle=0;
        return upr
    sgn=np.sign(CD11*CD22-CD12*CD21)
    if sgn>0:
        print('CD sign is positive')
    angle1=atan2(sgn*CD12,CD22)
    angle2=atan2(-CD21,sgn*CD11)
    up=np.array([0,0.3977,0.9175])
    if(abs(angle1-angle2)>0.005):
        print('Cannot understand rotation angle in  '+fitsfile)
    M=Ortho_Proj(E,up)
    Rz=np.array([[cos(angle1),-sin(angle1),0],[sin(angle1),cos(angle1),0],[0,0,1]])
    Mrot=np.dot(np.dot(M.transpose(),Rz.transpose()),M)
    upr=np.dot(Mrot,up.transpose()).transpose()
    return upr
    
        
        
    
    
        
        
            
    