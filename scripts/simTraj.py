from pathAdapt import pathSims
import numpy as np
import sys

L = 5*10**6
p = 10  # epi per gen in 10**4 pop size

sMin = 0.00001
sMax = 0.01
fac = 1.1

epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.2,mu=0.95,rec=0.05,L=L,mutRate = 2.5e-8)
#epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.01,mu=0.003,rec=0.23,,L=L,mutRate = 2.5e-8)
totBackground = 0
totEpiAdd = 0
totEpiAddMax = 0

a = epiSim.genSingleEvoEpi(-0.001,1./(2.*int(sys.argv[1])),0.00005,1)
for freq in a:
    print(freq)
