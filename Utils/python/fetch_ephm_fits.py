#!/usr/bin/python2
#Read ecliptic geocentric coordinates of the asteroid at the time instant JD0
#No LT correction
#Usage: asteroid number 'fits-files' outputfile
#E.g. ./fetch_ephm_fits.py 9 '*.fits' ephm.dat
#output file format time(JD) E0 E
import pyfits
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
  V=np.array([float(wlist[6]),float(wlist[7]),float(wlist[8]),JD0-float(wlist[12])])
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
  V=np.array([float(wlist[6]),float(wlist[7]),float(wlist[8]),JD0-float(wlist[12])])
  return V
asteroid=int(sys.argv[1])
print "Asteroid: "
print asteroid
print "output filename: "
input=sys.argv[2]
output=sys.argv[3]
outf=open(output,"w")
print output
file=glob.glob(input)
print "processing"
print len(file)
print "files"
tn=Telnet('ssd.jpl.nasa.gov', 6775)
for f in file:
  print "Processing file "
  print f
  fits_open=pyfits.open(f)
  hdu=fits_open[0]
  try:
    obstime=hdu.header['MJD-OBS']+2400000.5
  except:
      obstime=hdu.header['JDMEAN']
  V2=get_horizons(tn,asteroid,obstime) #Asteroid direction
  LT=V2[3]
  V1=get_horizons_sun(tn,10,LT) #Sun direction from Earth
  V3=V1[0:3]-V2[0:3] #Points from asteroid to Sun
  
  V2=-V2[0:3] #Points from asteroid to Earth
  outf.write(str(obstime)+" ")
  V3.tofile(outf,sep=" ")
  outf.write(" ")
  V2.tofile(outf,sep=" ")
  outf.write('\n')
tn.close()
outf.close()
