# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:05:09 2016

@author: matvii
"""

def parse_ini(inifile):
    fd=open(inifile,'r')
    while True:
        iline=fd.readline()
        if 'MinTim' in iline:
            MinTim=float(iline.split('=')[1].split(' ')[0])
            break
    while True:
        iline=fd.readline()
        if 'UseAO=' in iline:
            nAO=int(iline.split('=')[1].split(' ')[0])
            break
    Date=np.zeros(nAO)
    PixScale=np.zeros(nAO)
    AOSize=np.zeros(nAO)
    Filename=[]
    tao1='AOFile='
    tao2='Date'
    tao3='PixScale'
    tao4='AOSize'
    tao='AO'+str(1)
    while len(line)!=0:
        line=fd.readline()
        if tao in line:
            break
    for j in range(1,nAO+1):
        tao='[AO'+str(j)+']'
        taon='[AO'+str(j+1)+']'
        Date[j-1]=float('nan')
        while (taon not in line) and len(line)!=0:
            line=fd.readline()
            if len(line)>0 and line[0]=='#':
                line=fd.readline()
            if tao1 in line:
                Filename.append(line.split('=')[1].split(' ')[0])
            if tao2 in line:
                Date[j-1]=float(line.split('=')[1].split(' ')[0])
            if tao3 in line:
                PixScale[j-1]=float(line.split('=')[1].split(' ')[0])
            if tao4 in line:
                AOSize[j-1]=int(line.split('=')[1].split(' ')[0].split(',')[0])
    fd.close()
    fd=open(inifile)
    line=fd.readline()
    while len(line)!=0:
        if 'EphFile' in line:
            EphFile=line.split('=')[1].split(' ')[0]
        if 'ShapeFile' in line:
            ShapeFile=line.split('=')[1].split(' ')[0]
        if 'AnglesFile' in line:
            Anglefile=line.split('=')[1].split(' ')[0]
        line=fd.readline()
    fd.close()
    tlist,vlist=Read_Shape(ShapeFile)
    fd=open(Anglefile)
    line=fd.readline()
    angles=np.array([(90-float(line.split(' ')[0]))*pi/180,float(line.split(' ')[1])*pi/180,24*2*pi/float(line.split(' ')[2])])
    fd.close()
    cao=1.731446326742403e+02
    M=np.loadtxt(EphFile)
    E=np.zeros([nAO,3])
    E0=np.zeros([nAO,3])
    for j in range(0,nAO):
        if isnan(Date[j]):
            Date[j]=read_fits_date(Filename[j])+2400000.5
        if min(abs(M[:,0]-Date[j]))>0.01:
            sys.exit("No valid date in the ephemeris file")
        Ind=np.argmin(abs(M[:,0]-Date[j]))
        E[j,:]=np.array(M[Ind,4:8])
        E0[j,:]=np.array(M[Ind,1:4])
        Date[j]=Date[j]-np.sqrt(np.dot(E[j,:],E[j,:]))/cao-MinTim
    dist=np.zeros(nAO)
    km2arcsec=np.zeros(nAO)
    up=np.zeros([nAO,3])
    
    for j in range(0,nAO):
        E0[j,:]=E0[j,:]/np.sqrt(np.dot(E0[j,:],E0[j,:]))
        dist[j]=np.sqrt(np.dot(E[j,:],E[j,:]))
        E[j,:]=E[j,:]/dist[j]
        km2arcsec[j]=1/(dist[j]*149597871.0)*180/pi*3600.0
        up[j,:]=calculate_rotated_frame(Filename[j],E[j,:])

                
            
    
    
    
    
            