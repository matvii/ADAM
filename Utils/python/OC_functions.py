# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 11:59:10 2016

@author: matvii
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
def Ortho_Proj(view,up):
    view=np.array(view)
    up=np.array(up)
    z=view/np.sqrt((view*view).sum())
    x=np.cross(up,z)
    x=x/np.sqrt((x*x).sum())
    y=np.cross(z,x)
    y=y/np.sqrt((x*x).sum())
    M=np.zeros((3,3))
    M[0,:]=x
    M[1,:]=y
    M[2,:]=z
    return np.matrix(M)

def Rot_Matrix(bet0,lam0,omg0,t,phi0):
    bet=np.pi/180*(90-bet0)
    lam=np.pi/180*lam0
    omg=2*np.pi*24/omg0
    f=omg*t+phi0*np.pi/180
    cf=np.cos(f)
    sf=np.sin(f)
    cb=np.cos(bet)
    sb=np.sin(bet)
    cl=np.cos(lam)
    sl=np.sin(lam)
    fmat=np.matrix([[cf,sf,0],[-sf,cf,0],[0,0,1]])
    blmat=np.matrix([[cb*cl,cb*sl,-sb],[-sl,cl,0],[sb*cl,sb*sl,cb]])
    return fmat*blmat

def Read_OCC(ocfile,JD0):
    f=open(ocfile,'r')
    s=f.readline()
    nOCC=int(s)
    OCCs=[]
    E=np.zeros(3)
    V=np.zeros(3)
    Chords=[]
    Chordtype=[]
    TIME=[]
    ChordErr=[]
    Dist=np.zeros(nOCC)
    c=299792.458*60*60*24
    MinTim=2443846.0;
    for n in range(0,nOCC):
        s=f.readline()
        svert=s.split()
        E=np.array([float(svert[0]),float(svert[1]),float(svert[2])])
        norm=np.sqrt((E*E).sum())
        E=E/norm
        Dist=norm
        s=f.readline()
        svert=s.split()
        V=np.array([float(svert[0]),float(svert[1]),float(svert[2])])
        s=f.readline()
        s=f.readline()
        nchords=int(s)
        ch=np.zeros((nchords,4))
        chtype=np.zeros(nchords,dtype=np.int)
        tim=np.zeros((nchords,2))
        err=np.zeros((nchords,2))
        for j in range(0,nchords):
            s=f.readline()
            svert=s.split()
            [tim[j,0],err[j,0],ch[j,0],ch[j,1],tim[j,1],
             err[j,1],ch[j,2],ch[j,3]]=np.array([float(svert[0])-MinTim,float(svert[1]),
             float(svert[2]),float(svert[3]),float(svert[4])-MinTim,float(svert[5])
             ,float(svert[6]),float(svert[7])])
            chtype[j]=int(svert[8])
        mt=np.mean(tim)
        meanx=np.mean(ch[:,(0,2)])
        meany=np.mean(ch[:,(1,3)])
        ch[:,(0,2)]=ch[:,(0,2)]-meanx
        ch[:,(1,3)]=ch[:,(1,3)]-meany
        tim=tim-Dist/c
        OCCs.append(OCC(E,V,Dist,ch,chtype,tim,err))
    f.close()
    return OCCs
    
def Plot_Occ(tlist,vlist,angles,offset,OCCs,index):
    if index>len(OCCs):
        sys.exit("Index too large for occultation")
    up=np.array([0,0.3977,0.9175])
    M=Ortho_Proj(OCCs[index-1].E,up)
    R=Rot_Matrix(angles[0],angles[1],angles[2],OCCs[index-1].Meantime,0)
    MR=M*R.transpose()
    vlist2=(MR*vlist.transpose()).transpose()
    vlist2=vlist2[:,0:2]
    V=M*np.matrix(OCCs[index-1].V).transpose()
    v=V[0:2]
    vlist2[:,0]=vlist2[:,0]+offset[0]
    vlist2[:,1]=vlist2[:,1]+offset[1]
    #plt.figure(figsize=(50,50))
    plt.gca().set_aspect('equal')
    plt.triplot(np.array(vlist2[:,0]).flatten(),np.array(vlist2[:,1]).flatten(), tlist-1, 'g-',alpha=0.5)
    for j in range(0,len(OCCs[index-1].Chordtype)):
        a=OCCs[index-1].Chords[j,0:2]
        b=OCCs[index-1].Chords[j,2:4]
        ctype=OCCs[index-1].Chordtype[j]
        if ctype>=0:
            plt.plot([a[0],b[0]],[a[1],b[1]],'k-')
        if ctype<0:
            plt.plot([a[0],b[0]],[a[1],b[1]],'k--')
        #plt.savefig('/tmp/test.pdf',format='pdf')



    
 
    

    
    