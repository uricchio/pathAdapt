 for i in {1..100}; do for s in {0.01,0.001,0.0001,0.00001,0.0000001}; do python scripts/calcXf.py 10000 $i $s >> simData/simDeltaPfix/sim.deltaFix.txt ; done; done
