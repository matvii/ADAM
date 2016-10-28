#!/bin/python2
"""
Displays 3d view of the shape given by shapefile
Usage: ./Display_Shape shapefile
"""
import sys
from utils import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
input=sys.argv[1]
tlist,vlist=Read_Shape(input)
tlist=tlist-1
fig = plt.figure()
ax=fig.gca(projection='3d')
maxv=max(vlist.max(),abs(vlist.min()))
lim=maxv+0.2*maxv
ax.set_xlim3d(-lim,lim)
ax.set_ylim3d(-lim,lim)
ax.set_zlim3d(-lim,lim)
ax.plot_trisurf(vlist[:,0],vlist[:,1],vlist[:,2],triangles=tlist, edgecolor='none')


plt.show()