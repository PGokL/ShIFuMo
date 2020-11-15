#!/usr/bin/env python

import sys
import tarfile
from pylab import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy

fname = sys.argv[1]
coupe = sys.argv[2:5]

file=open(fname, 'r')

reader = vtk.vtkStructuredPointsReader()
reader.SetFileName(fname)
reader.Update()
data = reader.GetOutput()

#Grab the field in my vtk file
numpy_array = vtk_to_numpy(data.GetPointData().GetScalars())
Nx, Ny, Nz = data.GetDimensions()
    
nb = data.GetNumberOfScalarComponents()
f=reshape(numpy_array,(Nx,Ny,Nz),order='C')

if(coupe[0]==":"):                              
    plot(f[:,int(coupe[1]),int(coupe[2])],'o-')
if(coupe[1]==":"):
    plot(f[int(coupe[0]),:,int(coupe[2])],'o-')
if(coupe[2]==":"):
    plot(f[int(coupe[0]),int(coupe[1]),:],'o-')
                
legend(loc=0)

show()


