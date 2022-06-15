from pathAdapt import pathSims
import numpy as np
import sys

L = 5*10**6
p = 1  # epi per gen in 10**4 pop size

epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.008,mu=0.73,rec=0.23,L=L,mutRate = 2.5e-8)

# mutRatesGam calc
a =epiSim.stochSim(p*(float(sys.argv[1])/10000.), 10000, 0.00001, 0.05, 0.00001, 0.00005,1)
print(a[0],a[1],a[2])
