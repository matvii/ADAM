#!/usr/bin/python
import numpy as np
import argparse
from utils import *
parser=argparse.ArgumentParser()
parser.add_argument('-T','--JD0',type=float,required=True)
parser.add_argument('-O','--OCfile',type=str,required=True)
parser.add_argument('-S','--Shapefile',type=str)
parser.add_argument('-A','--Angles',type=str,required=True)
parser.add_argument('-o','--Offset',type=str)
parser.add_argument('-I','--Index',type=int)
parser.add_argument('-E','--Ellipsoid',type=str)
args=parser.parse_args()
T=args.JD0
OCfile=args.OCfile
shapefile=args.Shapefile
index=1
angles=np.asarray([float(i) for i in args.Angles.split(',')])
if args.Offset is not None:
    offset=np.asarray([float(i) for i in args.Offset.split(',')])
else:
    offset=np.array([0.0,0.0])
if args.Index is not None:
    index=args.Index
if args.Ellipsoid is not None:
    ell=np.asarray([float(i) for i in args.Ellipsoid.split(',')])
    tlist,vlist=Generate_Ellipsoid(ell[0],ell[1],ell[2],10)
if args.Shapefile is not None:
    tlist,vlist=Read_Shape(shapefile)
OCCs=Read_OCC(OCfile,T)
Plot_Occ(tlist,vlist,angles,offset,OCCs,index)
