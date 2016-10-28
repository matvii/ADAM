import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
def Read_Shape(shapefile):
    f=open(shapefile,'r')
    s=f.readline()
    nvert=int(s.split()[0])
    nfac=int(s.split()[1])
    vlist=np.zeros((nvert,3))
    tlist=np.zeros((nfac,3),dtype=np.int)
    tlist2=[]
    #Read vertex list
    for i in range(0,nvert):
        s=f.readline()
        svert=s.split()
        vlist[i,]=np.asarray([float(svert[0]),float(svert[1]),float(svert[2])])
    for i in range(0,nfac):
        s=f.readline()
        svert=s.split()
        tlist[i,]=np.asarray([int(svert[0]),int(svert[1]),int(svert[2])],dtype=np.int)
    f.close()
    return tlist,vlist
class OCC:
    def __init__(self,E0,V0,Dist,ch,chtype,tim,err):
        self.E=np.array(E0)
        self.V=np.array(V0)
        self.D=Dist
        self.Chords=np.array(ch)
        self.Chordtype=np.array(chtype)
        self.TIME=np.array(tim)
        self.Err=np.array(err)
        self.Meantime=tim.mean()

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
    return tlist,vlist

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
    return M

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
    fmat=np.array([[cf,sf,0],[-sf,cf,0],[0,0,1]])
    blmat=np.array([[cb*cl,cb*sl,-sb],[-sl,cl,0],[sb*cl,sb*sl,cb]])
    return np.dot(fmat,blmat)

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
    MR=np.dot(M,R.transpose())
    vlist2=np.dot(MR,vlist.transpose()).transpose()
    vlist2=vlist2[:,0:2]
    V=np.dot(M,OCCs[index-1].V.transpose())
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
    plt.show()
        #plt.savefig('/tmp/test.pdf',format='pdf')


def v2fmatrix(tlist,vlist):
    """
    Matrix indexing vertices->facets
    """
    nvert=vlist.shape[0]
    nfac=tlist.shape[0]
    M=np.zeros((nvert,nfac),dtype=np.bool)
    for i in range(0,nfac):
        M[tlist[i,:],i]=1
    return M

def Local_Horizon(tlist,vlist,facindex,v2fmatrix):
    """
    Returns facets over the local horizon of  facet facindex.
    v2fmatrix is matrix from v2fmatrix funcion
    """
    nvert=vlist.shape[0]
    nfac=tlist.shape[0]
    normal=np.cross(vlist[tlist[facindex,1],:]-vlist[tlist[facindex,0],:],vlist[tlist[facindex,2],:]-vlist[tlist[facindex,0],:])
    cent=(vlist[tlist[facindex,0],:]+vlist[tlist[facindex,1],:]+vlist[tlist[facindex,2],:])/3
    nlen=np.sqrt(np.dot(normal,normal))
    normal=normal/nlen
    res=np.dot(vlist-cent,normal)
    res[tlist[facindex,:]]=0
    vh=res>0
    #vec=np.asarray(range(0,nfac))
    return sum(v2fmatrix[vh,:])>0

def is_in_triangle(point,dir,v1,v2,v3):
    """
    Check if the line segment determined by point and vector dir
    intersects the triangle determined by (v1,v2,v3)        
    """
    s1=v1-v2
    s2=v1-v3
    b=v1-point
    det=s1[0]*s2[1]*dir[2]-s1[0]*dir[1]*s2[2]+s2[0]*dir[1]*s1[2]-s2[0]*s1[1]*dir[2]+dir[0]*s1[1]*s2[2]-dir[0]*s2[1]*s1[2]
    if abs(det)<1e-6:
        return 0
    gamma=(dir[2]*(s1[0]*b[1]-b[0]*s1[1])+dir[1]*(b[0]*s1[2]-s1[0]*b[2])+dir[0]*(s1[1]*b[2]-b[1]*s1[2]))/det
    if gamma<0 or gamma>1:
        return 0
    beta=(b[0]*(s2[1]*dir[2]-dir[1]*s2[2])+b[1]*(dir[0]*s2[2]-s2[0]*dir[2])+b[2]*(s2[0]*dir[1]-s2[1]*dir[0]))/det
    if beta<0 or beta>1-gamma:
        return 0
    return 1
