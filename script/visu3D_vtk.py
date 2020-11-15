#!/usr/bin/env python3

from pylab import *
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from mayavi import mlab
       

#................................ c1 field .................................
fname=sys.argv[1]
file=open(fname, 'r')

reader = vtk.vtkStructuredPointsReader()
reader.SetFileName(fname)
reader.Update()
    
data = reader.GetOutput()

#Grab the field in my vtk file
numpy_array = vtk_to_numpy(data.GetPointData().GetScalars())
Nx, Ny, Nz = data.GetDimensions()
    
f=reshape(numpy_array,(Nx,Ny,Nz),order='C')


fig3d = mlab.figure()           
fmin=f[:,:,:].min()
fmax=f[:,:,:].max()            
    
mlab.contour3d(f[:,:,:], contours=[(fmin+fmax)/2.0])

fig3d.scene.background = (1.,1.,1.)                             #define background and foreground colors
fig3d.scene.foreground = (0.3,0.3,0.3)     

mlab.outline(extent=[0,Nx,0,Ny,0,Nz])                           #plot box boudaries
origin = (float(Nx)/2.0,float(Ny)/2.0,float(Nz)/2.0) 
fig3d.scene.camera.focal_point = origin
    

#mlab.savefig("test.png")
mlab.show()                                                     #draw the prepresentation

