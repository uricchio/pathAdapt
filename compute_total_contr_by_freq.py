import sys
import numpy as np

# sys.argv[1] should be an output file from my simulation code in pathSims.py

fh = open(sys.argv[1],'r')

sums = {}
vals = {}
pfix0 = {}
pfixN = {}
for line in fh:
    data = line.strip().split()
    if len(data) < 10:
        break
    if float(data[3]) not in sums:
        sums[float(data[3])] = [float(data[6])*float(data[9])]
        vals[float(data[3])] = [data[7], data[8]]
        pfix0[float(data[3])] = [float(data[4])*float(data[9])]
        pfixN[float(data[3])] = [float(data[5])*float(data[9])]
    else: 
        sums[float(data[3])].append(float(data[6])*float(data[9]))
        pfix0[float(data[3])].append(float(data[4])*float(data[9]))
        pfixN[float(data[3])].append(float(data[5])*float(data[9]))

for s in sums:
    print (s, np.sum(sums[s]), vals[s][0], vals[s][1], np.sum(pfix0[s]), np.sum(pfixN[s]))
    

