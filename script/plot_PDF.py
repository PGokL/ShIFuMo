#!/usr/bin/env python

from pylab import *
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy

#ind = sys.argv[1]
#fname="ppts_"+ind.zfill(6)+".out"

fname = sys.argv[1]
file=open(fname, 'r')
flag,posi,posj,posk,Rall = loadtxt(file, usecols=[1,2,3,4,5], unpack=True)

R = list(Rall*flag)
R[:] = (value for value in R if value > 0)
R = array(R)

#print len(R), "ppts"

R_ave = average(R)
R_std = std(R)

hist(R/R_ave,bins=15,density=True,label="$\\langle R \\rangle=$"+str(round(R_ave,2))+";std="+str(round(R_std,2)))

x=arange(0.0,1.5,0.001)

LSW = (3**4*exp(1.)/2**(5./3.)) *x*x* exp(-1./(1-2./3.*x)) / ((x+3)**(7./3.) * (1.5-x)**(11./3.))

plot(x,LSW,color="k",ls='--',label="LSW")

xlabel("$R/\\langle R \\rangle$",fontsize=20)
ylabel("PDF",fontsize=20)

xlim([0.0,2.0])
ylim([0.0,2.5])

legend(loc=0)

show()

#savefig("R_pdf_"+ind.zfill(6)+".png")
#savefig("R_pdf.png")
