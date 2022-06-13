from pathAdapt import pathSims
import numpy as np
import sys

L = 5*10**6
p = 10  # epi per gen in 10**4 pop size

epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.1,mu=0.01,rec=0.03,L=L,mutRate = 0.001/(4*10000))
totBackground = 0
totEpiAdd = 0

# number of epi events per gen in rescaled pop
kN = 0
k10000 = 0

# mutRatesGam calc
epiSim.mutRatesGamma(0.0001,0.01,0.0001)

for s in np.arange(0.0001,0.01,0.0001):
    s *= -1
    De = epiSim.fixOverAllFreq(s,0.0005,1)
    
    # total number of expected fixations per gen over all seg. sites in absence of epi
    totBackground += De[1]

    # total number of additional fixations per gen with 1 standing resis variant
    totEpiAdd += De[0]

    print(De[1],De[0],epiSim.mutRatesGam[round(abs(s),15)])
    
    # number of seg sites calc
    kN += epiSim.mutRatesGam[round(abs(s),15)]*epiSim.totSites(s,int(sys.argv[1]))
    k10000 += epiSim.mutRatesGam[round(abs(s),15)]*epiSim.totSites(s,10000)

print((p*(kN/k10000)*totEpiAdd)/(totBackground+p*(kN/k10000)*totEpiAdd),p*(kN/k10000))
