rm simData/alpha*txt
for mu in {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.999,1.0}; do python scripts/sim.py 10000 $mu 0.0000005 2 >> simData/alpha.nearlyNeutral.txt ; done
for mu in {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.999,1.0}; do python scripts/sim.py 10000 $mu 0.00005 1 >> simData/alpha.diminishingReturns.txt ; done
for mu in {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.999,1.0}; do python scripts/sim.py 10000 $mu 0.0005 3 >> simData/alpha.severe.txt ; done
