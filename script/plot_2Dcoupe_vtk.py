#!/usr/bin/env python

from pylab import *
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy


fname=sys.argv[1]
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

figure()   
if(coupe[0]!=":"):
    im=imshow(transpose(f[int(coupe[0]),:,:]),origin='lower')
if(coupe[1]!=":"):
    im=imshow(transpose(f[:,int(coupe[1]),:]),origin='lower')
if(coupe[2]!=":"):
    im=imshow(transpose(f[:,:,int(coupe[2])]),origin='lower')
cbar = colorbar(im, orientation='vertical', shrink=0.5, aspect=20)
            
show()
