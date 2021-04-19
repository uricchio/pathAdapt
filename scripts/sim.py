from pathAdapt import pathSims
import numpy as np

epiSim = pathSims.SimulateEpi()

for j in range(6,61):
    s = -1.*(1.3**-j)
    epiSim.DiscSFSSelNeg(s)
    for i in range(1,2000):
        f = i/2000.
        epiSim.initAllele(s,f)
        epiSim.determSim()

