#!/usr/bin/python3
#Read ecliptic geocentric coordinates of the asteroid at the time instant JD0
#No LT correction
#Usage: asteroid number 'fits-files' outputfile
#E.g. ./fetch_ephm_fits.py 9 '*.fits' ephm.dat
#output file format time(JD) E0 E
import pyfits
import sys
import time
import numpy as np
import glob
import urllib.request
import re
def get_horizons(asteroid,JD0):
  #JD0=2456575.50
    url='http://vo.imcce.fr/webservices/miriade/ephemcc_query.php?-name=a:'+str(asteroid)+'&-ep='+str(JD0)+'&-observer=500&-tcoor=2&-rplane=2'
    resp=urllib.request.urlopen(url)
    encoding = resp.headers.get_content_charset('utf-8')
    cont=resp.read()
    htext=cont.decode(encoding)
    res=htext.split('\n')
    for l in res:
        if l.startswith('<vot:TR><vot:TD>'):
            ldata=l
            break

    data=ldata.split("<vot:TD>")
    V=np.zeros(3)

    
    V[0]=float(re.search(r'-*\d+.?\d*',data[3]).group())
    V[1]=float(re.search(r'-*\d+.?\d*',data[4]).group())
    V[2]=float(re.search(r'-*\d+.?\d*',data[5]).group())
    
  
    return V
def get_horizons_sun(JD0):
    url='http://vo.imcce.fr/webservices/miriade/ephemcc_query.php?-name=p:Earth&-ep='+str(JD0)+'&-observer=@sun&-tcoor=2&-rplane=2'
    resp=urllib.request.urlopen(url)
    encoding = resp.headers.get_content_charset('utf-8')
    cont=resp.read()
    htext=cont.decode(encoding)
    res=htext.split('\n')
    for l in res:
        if l.startswith('<vot:TR><vot:TD>'):
            ldata=l
            break

    data=ldata.split("<vot:TD>")
    V=np.zeros(3)

    
    V[0]=-float(re.search(r'-*\d+.?\d*',data[3]).group())
    V[1]=-float(re.search(r'-*\d+.?\d*',data[4]).group())
    V[2]=-float(re.search(r'-*\d+.?\d*',data[5]).group())
    
    return V
if len(sys.argv) != 4:
    print("3 arguments required: asteroid number,fits files, outputfile")
    sys.exit(-1)
asteroid=int(sys.argv[1])
print('Asteroid number: '+str(asteroid))
input=sys.argv[2]
output=sys.argv[3]
print ('output filename: '+output)

outf=open(output,"w")

file=glob.glob(input)
print('processing '+str(len(file))+' files')

for f in file:
  print('Processing file '+ f)
  fits_open=pyfits.open(f)
  hdu=fits_open[0]
  try:
    obstime=hdu.header['MJD-OBS']+2400000.5
  except:
      obstime=hdu.header['JDMEAN']
  print('Object: '+hdu.header['OBJECT'])
  V2=get_horizons(asteroid,obstime) #Asteroid direction
  LT=obstime
  V1=get_horizons_sun(LT) #Sun direction from Earth
  V3=V1[0:3]-V2[0:3] #Points from asteroid to Sun
  
  V2=-V2[0:3] #Points from asteroid to Earth
  outf.write(str(obstime)+" ")
  V3.tofile(outf,sep=" ")
  outf.write(" ")
  V2.tofile(outf,sep=" ")
  outf.write('\n')

outf.close()
