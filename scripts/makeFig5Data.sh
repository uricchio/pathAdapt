# for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.75 0.00005 1 >> simData/Fig5Data.dimReturns.txt; done
# for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.75 0.0000005 2 >> simData/Fig5Data.nearNeutral.txt; done
# for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.75 0.0005 3 >> simData/Fig5Data.severe.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.5 0.00005 1 >> simData/Fig5Data.mid.dimReturns.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.5 0.0000005 2 >> simData/Fig5Data.mid.nearNeutral.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.5 0.0005 3 >> simData/Fig5Data.mid.severe.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.25 0.00005 1 >> simData/Fig5Data.low.dimReturns.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.25 0.0000005 2 >> simData/Fig5Data.low.nearNeutral.txt; done
 for N in {100,200,500,1000,2000,5000,10000,20000,50000}; do python scripts/sim.py $N 0.25 0.0005 3 >> simData/Fig5Data.low.severe.txt; done
