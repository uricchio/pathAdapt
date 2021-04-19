import sys
import os
import numpy
import re
from collections import defaultdict

"""
A class for making various input files and directories that are needed for the pipeline
"""

class SimulateEpi(N=500,I=1,R=0,beta=0.1,mu=0.1):

    def __init__():
        self.N = N
        self.I = I
        self.R = R
        self.beta = beta
        self.mu = mu
        self.rec = 1-mu        
        self.raf = 0


    def stochSim():
        while self.I > 0:
            self.I -= numpy.random.binom(self.mu,self.I)   
