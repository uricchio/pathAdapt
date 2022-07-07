from pathAdapt import pathSims
import numpy as np
import sys

if len(sys.argv) < 5:
    print("usage: python", sys.argv[0], "<populationSize> <deathRate> <sHalf> <scale>")
    exit() 


L = 10**6
p = 10  # epi per gen in 10**4 pop size

mu = float(sys.argv[2])
#rec = float(sys.argv[3])
rec = 1. - mu
sHalf = float(sys.argv[3])
scale = float(sys.argv[4])

# these params determine the grid of selection strengths over which we will compute
sMin = 0.000000001
sMax = 0.005
fac = 1.1

epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.0005,mu=mu,rec=rec,L=L,mutRate = 2.5e-8)
#epiSim = pathSims.SimulateEpi(N=int(sys.argv[1]),beta=0.01,mu=0.003,rec=0.23,,L=L,mutRate = 2.5e-8)
totBackground = 0
totEpiAdd = 0
totEpiAddMax = 0

# number of epi events per gen in rescaled pop
kN = int(sys.argv[1])+0.
k10000 = 10000

# mutRatesGam calc
epiSim.mutRatesGamma(sMin,sMax,fac)
epiSim.relNumMutsGamma(sMin,sMax,fac)

s = -sMin
svals = []
while s > -sMax:
    De = epiSim.fixOverAllFreq(s,sHalf,scale)
    
    # total number of expected fixations per gen over all seg. sites in absence of epi
    totBackground += De[1]

    # total number of additional fixations per gen with 1 standing resis variant
    totEpiAdd += De[0]
    totEpiAddMax += De[2] 
 
    svals.append([s,(((p*(kN/k10000)*totEpiAdd)+totBackground)/epiSim.mutRatesGam[round(abs(s),15)])-(totBackground/epiSim.mutRatesGam[round(abs(s),15)])])
   
    s *= fac

#for i in range(len(svals)):
#    svals[i][1] /= totEpiAdd

for s in svals:
    for thing in s:
        print(thing,end=' ')
    print()

#print((p*(kN/k10000)*totEpiAdd)/(totBackground),(p*(kN/k10000)*totEpiAddMax)/(totBackground+p*(kN/k10000)*totEpiAddMax),sys.argv[1],mu,sHalf)
