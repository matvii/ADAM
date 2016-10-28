# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 10:51:21 2016

@author: matvii
"""
def Triangulate_Sphere(nrows):
    from numpy import pi
    nvert=4*nrows**2+2
    nfac=8*nrows**2
    t=np.zeros(nvert)
    f=np.zeros(nvert)
    dth=pi/(2*nrows)
    k=0
    t[0]=0
    f[0]=0
    for i in range(1,nrows+1):
        dph=pi/(2*i)
        for j in range(0,4*i):
            k=k+1
            t[k]=i*dth
            f[k]=j*dph
    for i in range(nrows-1,0,-1):
        dph=pi/(2*i)
        for j in range(0,4*i):
            k=k+1
            t[k]=pi-i*dth
            f[k]=j*dph
    t[k+1]=pi
    f[k+1]=0
    ntri=0
    ifp=np.zeros((8*nrows**2,3),int)
    nod=np.zeros((2*nrows+1,4*nrows+1),int)
    nnod=1
    nod[0,0]=nnod
    for i in range(1,nrows+1):
        for j in range(0,4*i):
            nnod=nnod+1
            nod[i,j]=nnod
            if j==0:
                nod[i,4*i]=nnod
    
    for i in range(nrows-1,0,-1):
        for j in range(0,4*i):
            nnod=nnod+1
            nod[2*nrows-i,j]=nnod
            if j==0:
                nod[2*nrows-i,4*i]=nnod
    nod[2*nrows,0]=nnod+1
    ntri=-1
    for j1 in range(1,nrows+1):
        for j3 in range(1,5):
            j0=(j3-1)*j1
            ntri=ntri+1
            ifp[ntri,0]=nod[j1-1,j0-(j3-1)]
            ifp[ntri,1]=nod[j1,j0]
            ifp[ntri,2]=nod[j1,j0+1]
            for j2 in range(j0+1,j0+j1):
                ntri=ntri+1
                ifp[ntri,0]=nod[j1,j2]
                ifp[ntri,1]=nod[j1-1,j2-(j3-1)]
                ifp[ntri,2]=nod[j1-1,j2-1-(j3-1)]
                ntri=ntri+1
                ifp[ntri,0]=nod[j1-1,j2-(j3-1)]
                ifp[ntri,1]=nod[j1,j2]
                ifp[ntri,2]=nod[j1,j2+1]
    
    
    for j1 in range(nrows+1,2*nrows+1):
        for j3 in range(1,5):
            j0=(j3-1)*(2*nrows-j1)
            ntri=ntri+1
            ifp[ntri,0]=nod[j1,j0]
            ifp[ntri,1]=nod[j1-1,j0+1+(j3-1)]
            ifp[ntri,2]=nod[j1-1,j0+(j3-1)]
            for j2 in range(j0+1,j0+(2*nrows-j1)+1):
                ntri=ntri+1
                ifp[ntri,0]=nod[j1,j2]
                ifp[ntri,1]=nod[j1-1,j2+(j3-1)]
                ifp[ntri,2]=nod[j1,j2-1]
                ntri=ntri+1
                ifp[ntri,0]=nod[j1,j2]
                ifp[ntri,1]=nod[j1-1,j2+1+(j3-1)]
                ifp[ntri,2]=nod[j1-1,j2+(j3-1)]
    return ifp,t,f

def Generate_Ellipsoid(a,b,c,nrows):
    tlist,t,f=Triangulate_Sphere(nrows)
    x=a*np.sin(t)*np.cos(f)
    y=b*np.sin(t)*np.sin(f)
    z=c*np.cos(t)
    vlist=np.zeros((len(x),3))
    vlist[:,0]=x
    vlist[:,1]=y
    vlist[:,2]=z
    tlist=tlist-1
    return tlist,vlist

        