#!/usr/bin/python2
#Read ecliptic geocentric coordinates of the asteroid at the time instant JD0
#No LT correction
#Usage: asteroid number time_file outputfile
#time file contains observation times in JD
#output file format time(JD) E0 E LT_TIME
import sys
from telnetlib import Telnet
import time
import numpy as np
import glob
def get_horizons(tn,asteroid,JD0):
  #JD0=2456575.50
  JD1=JD0+1;
  sJD0=str(JD0);
  sJD1=str(JD1);
  tn.read_until("Horizons>")
  tn.write(str(asteroid)+";"+"\n")
  tn.read_until("<cr>:", 2.)
  tn.write("e\n")
  tn.read_until(":", 2.)
  tn.write("v\n")
  res=tn.read_until(":", 2.)
  if res.find("Use previous"):
    tn.write("\n")
  else:
    tn.write("geo\n")
  tn.read_until(":", 2.)
  tn.write("eclip\n")
  tn.read_until("] : ", 2.)
  tn.write("JD "+sJD0+'\n')
  tn.read_until("] : ", 2.)
  tn.write("JD "+sJD1+'\n')
  tn.read_until(" : ", 2.)
  tn.write("1d\n")
  tn.read_until(" : ", 2.)
  #tn.write("\n")
  tn.write("n\n")
  tn.write("\n")
  tn.write("2\n")
  tn.write("\n")
  tn.write("\n")
  tn.write("\n")
  tn.write("\n")
  tn.write("\n")
  tn.read_until("$$SOE")
  eph=tn.read_until("$$EOE")
  tn.read_until("[R]edisplay, ? :",2.)
  tn.write("N\n")
  wlist=eph.split()
  xc=wlist.index('X')
  if wlist[xc+1]=='=':
      xval=wlist[xc+2]
  else:
          xval=wlist[xc+1].replace('=','')
  yc=wlist.index('Y')
  if wlist[yc+1]=='=':
      yval=wlist[yc+2]
  else:
          yval=wlist[yc+1].replace('=','')
  zc=wlist.index('Z')
  if wlist[zc+1]=='=':
      zval=wlist[zc+2]
  else:
          zval=wlist[zc+1].replace('=','')
  
  
  tc=wlist.index('LT=')+1
  V=np.array([float(xval),float(yval),float(zval),JD0-float(wlist[tc])])
  return V
def get_horizons_sun(tn,asteroid,JD0):
  #JD0=2456575.50
  JD1=JD0+1;
  sJD0=str(JD0);
  sJD1=str(JD1);
  tn.read_until("Horizons>")
  tn.write(str(asteroid)+"\n")
  tn.read_until("<cr>:", 2.)
  tn.write("e\n")
  tn.read_until(":", 2.)
  tn.write("v\n")
  res=tn.read_until(" : ", 2.)
  if res.find("Use previous"):
    tn.write("\n")
  else:
    tn.write("geo\n")
  tn.read_until(" : ", 2.)
  tn.write("eclip\n")
  tn.read_until(" : ", 2.)
  tn.write("JD "+sJD0+'\n')
  tn.read_until(" : ", 2.)
  tn.write("JD "+sJD1+'\n')
  tn.read_until(" : ", 2.)
  tn.write("1d\n")
  tn.read_until(" : ", 2.)
  tn.write("\n")
  tn.read_until("$$SOE")
  eph=tn.read_until("$$EOE")
  tn.read_until("[R]edisplay, ? :",2.)
  tn.write("N\n")
  wlist=eph.split()
  xc=wlist.index('X')
  if wlist[xc+1]=='=':
      xval=wlist[xc+2]
  else:
          xval=wlist[xc+1].replace('=','')
  yc=wlist.index('Y')
  if wlist[yc+1]=='=':
      yval=wlist[yc+2]
  else:
          yval=wlist[yc+1].replace('=','')
  zc=wlist.index('Z')
  if wlist[zc+1]=='=':
      zval=wlist[zc+2]
  else:
          zval=wlist[zc+1].replace('=','')
  
  
  tc=wlist.index('LT=')+1
  V=np.array([float(xval),float(yval),float(zval),JD0-float(wlist[tc])])
  return V
if len(sys.argv) != 4:
    print '3 arguments required: asteroid number,file with observation dates, outputfile'
    sys.exit(-1)

asteroid=int(sys.argv[1])
print "Asteroid: "
print asteroid
print "output filename: "
input=sys.argv[2]
output=sys.argv[3]
outf=open(output,"w")
print output
#file=glob.glob(input)
f=open(input,'r')
print "processing"
print f

tn=Telnet('ssd.jpl.nasa.gov', 6775)
for line in f:
  time=float(line)
  print "Processing time "
  print line
  #fits_open=pyfits.open(f)
  #hdu=fits_open[0]
  #time=hdu.header['MJD-OBS']+2400000.5
  V2=get_horizons(tn,asteroid,time) #Asteroid direction
  LT=V2[3]
  V1=get_horizons_sun(tn,10,LT) #Sun direction from Earth
  V3=V1[0:3]-V2[0:3] #Points from asteroid to Sun
  
  V2=-V2[0:3] #Points from asteroid to Earth
  outf.write(str(time)+" ")
  V3.tofile(outf,sep=" ")
  outf.write(" ")
  V2.tofile(outf,sep=" ")
  outf.write(" ")
  LT.tofile(outf,sep=" ")
  outf.write('\n')
tn.close()
outf.close()
