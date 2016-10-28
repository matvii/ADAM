def read_shape(shapefile):
    import numpy as np
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
    
tlist,vlist=read_shape('/tmp/mshape.txt')
from mayavi.mlab import *
triangular_mesh(vlist[:,0],vlist[:,1],vlist[:,2],tlist-1)

# Plot 2d triangle
import matplotlib.pyplot as plt
import matplotlib.tri as tri
plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(vlist[:,0], vlist[:,1], tlist-1, 'g-')
