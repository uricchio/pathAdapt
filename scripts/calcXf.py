from pathAdapt import pathSims
import numpy as np
import sys

L = 5*10**6
p = 10  # epi per gen in 10**4 pop size

sMin = 0.00001
sMax = 0.01
fac = 1.1

epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.0005,mu=0.95,rec=0.05,L=L,mutRate = 2.5e-8)

xi = float(sys.argv[2])/100.
s = -1.*float(sys.argv[3])
xf = epiSim.calcXf(s,xi,0.00005,1)

print(xi,xf,epiSim.deltaFixProb(s,xi,xf),s)
