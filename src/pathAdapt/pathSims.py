import sys
import os
import numpy as np
import re
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from collections import defaultdict

# to do: fix the beta in homozygous ind, which is reduced too much even when the tradeoff is very weak

"""
A class for simulating single epidemics

The idea is to:
    1.  assume a population at equilibirum
    2.  sample a resistance allele from the deleterious allele spectrum
    3,  suppose this allele is now subject to mortality/fecundity tradeoff, or maybe just additional fecundity cost for getting sick?
    4.  ???
    4.  profit

TO do:
    1. reweight by prob of being polymorphic at the given freqeuncy to get expected excess substituions
    2. try different tradeoff fuctions
    3. try different beta/mu

"""

class SimulateEpi():

    def __init__(self,N=1000,I0=1,I1=0,I2=0,R0=0,R1=0,R2=0,S0=998,S1=1,beta=0.003,mu=0.5,rec=0.1,L=1,al=0.184,be=0.000402):
        self.N = N
        self.I0 = I0
        self.I1 = I1
        self.I2 = I2
        self.S0 = S0
        self.S1 = S1
        self.S2 = N-I0-I1-S0-S1
        self.R0 = R0
        self.R1 = R1
        self.R2 = R2
        self.beta = beta
        self.mu = mu
        self.rec = rec        
        self.raf = 0
        self.muts = np.zeros(L)     # array of allele frequncies
        self.res = np.zeros(L)      # array storing which are resistance alleles
        self.resi = 0
        self.effects = np.zeros(L)  # array effect sizes
        self.L = L
        self.al = al
        self.be = be
        self.sfs = []

    # kimura SFS NOT normalized
    def kimuraSFS(self,s,x):
        if abs(s) > 10e-10:
            return np.float(mpmath.exp(2*self.N*s)*(1-mpmath.exp((-2*self.N*s)*(1-x)))/((mpmath.exp(2*s*self.N)-1)*x*(1-x)))
        return 1./x

    # SFS for gamma distributed selection coeffs
    def fullNeg(self,x):
        return np.float((2.**-self.al)*(self.be**self.al)*(-mpmath.zeta(self.al,x+self.be/2.) + mpmath.zeta(self.al,(2+self.be)/2.))/((-1.+x)*x))
       
    # FULL sfs over all frequencies, for gamma
    def DiscSFSSelGam(self):
        NN = int(round(2*self.N))
        dFunc = np.vectorize(self.fullNeg)
        return dFunc([i/(NN+0.) for i in range(1,NN)])
   
    # FULL sfs over all frequencies, for single s
    def DiscSFSSelNeg(self,s):
        NN = int(round(2*self.N))
        dFunc = np.vectorize(self.kimuraSFS)
        a = dFunc(s,[i/(NN+0.) for i in range(1,NN)])
        b = np.sum(a)
        self.sfs = np.divide(a,b)
        return self.sfs

    def DiscSFSSelNegFreq(self,s,j):
        NN = int(round(2*self.N))
        dFunc = np.vectorize(self.kimuraSFS)
        a = dFunc(s,[i/(NN+0.) for i in range(1,NN)])
        b = np.sum(a)
        return np.divide(a,b)[j]

    # get effect sizes and frequencies for mutations, sampling from SFS at equil
    # this sampling is currently WRONG -- need to condition on being polymorhpic
    def initPop(self):
        # select effects from a gamma distribution for effects
        for i in range(0,self.L):
            self.effects[i] = -1.*np.random.gamma(self.al,1./self.be)/(self.N*2)          
     
        # select frequencies randomly from SFS conditional on being polymorphic    
        for i in range(0,self.L):
            self.muts[i] = np.random.choice(range(1,2*self.N),1,p=self.DiscSFSSelNeg(self.effects[i]))
        #print(self.muts)

    def initAllele(self,s,f):
         for i in range(0,self.L):
            self.effects[i] = s
            self.muts[i] = f*2.*self.N


    # need a function to define how the tradeOff works
    # i.e., transform costs into potential benefits during infection
    def tradeOff(self,s):

        numer = (self.al/(self.be*2*self.N))+(1000000./self.N)*s
        denom = self.al/(self.be*2*self.N)
        
        if (numer/denom)*self.beta > 0:
            return (numer/denom)*self.beta 
        return 0.

    def sampleRes(self):
    
        i = np.random.randint(0,len(self.muts))
        self.res[i] = 1
        self.resi = i
  
    def initInfect(self, effect):
        totalProp = self.S0*self.beta+self.S1*self.tradeOff(effect)+self.S2*0.5*self.tradeOff(effect)
        if totalProp == 0.:
            totalProp = self.S1+self.S2+self.S0
            self.I0 = self.S0/totalProp
            self.I1 = self.S1/totalProp
            self.I2 = self.S2/totalProp
            return
        self.I0 = self.S0*self.beta/totalProp
        self.I1 = self.S1*self.tradeOff(effect)/totalProp
        self.I2 = self.S2*0.5*self.tradeOff(effect)/totalProp

    # simulate an epidemic, given the frequencies of the res alleles
    def determSim(self):
 
        #sample resistance allele
        self.sampleRes()

        # first, get afs of resistance alleles
        afs = []
        for i in range(0,len(self.muts)):
            if self.res[i]:
                afs.append([self.muts[i],self.effects[i]])
                #afs.append([50,-0.01])

        # need to consider how to account for diploid state here -- 3 classes for 0,1,2 copies of derived allele
        p = afs[0][0]/(2.*self.N)
        self.S2 = int(math.floor(p**2*self.N))
        self.S1 = int(math.ceil(2*p*(1-p)*self.N))-self.I1-self.I0 
        self.S0 = max(0,self.N-(self.S1+self.I0+self.I1+self.S2+self.I2))

        effect = afs[0][1]
        self.initInfect(effect)

        #print(self.S0,self.S1,self.S2,self.I0,self.I1,self.I2,self.N)

        #print (self.N, self.I0, self.I1, self.I2, self.R0, self.R1, self.R2, self.S0, self.S1, self.S2, self.I0+self.I1+self.R0+self.R1+self.S0+self.S1+self.S2+self.I2+self.R2)

        # differential eqns
        def SIR(t, z, mu, beta0, beta1, r):
            I0, I1, I2, S0, S1, S2, R0, R1, R2, N = z
            return [-r*I0-mu*I0+S0*beta0*(I0+I1+I2),-r*I1-mu*I1+S1*beta1*(I0+I1+I2), -r*I2-mu*I2+S2*0.5*beta1*(I0+I1+I2),
                     -S0*beta0*(I0+I1+I2), -S1*beta1*(I0+I1+I2), -S2*0.5*beta1*(I0+I1+I2),
                      r*I0, r*I1, r*I2,
                      -mu*(I0+I1+I2)]
     
        radius = 1
        tmax = 25
        while radius > 1e-5:
            tmax *= 2
            sol = solve_ivp(SIR, [0, tmax], [self.I0, self.I1, self.I2, self.S0, self.S1, self.S2, self.R0, self.R1, self.R2, self.N], 
                            args=(self.mu, self.beta, self.tradeOff(afs[0][1]), self.rec), dense_output=True)
            t0 = np.linspace(0, tmax, 1000)
            radius = sum(np.absolute(SIR(tmax,np.transpose(sol.sol(t0))[-1],self.mu, self.beta, self.tradeOff(afs[0][1]), self.rec)))

        RecCop = sol.sol(t0)[7][-1]+2*sol.sol(t0)[8][-1]
        SusCop = sol.sol(t0)[4][-1]+2*sol.sol(t0)[5][-1]
        InfCop = sol.sol(t0)[1][-1]+2*sol.sol(t0)[2][-1]

        newFreq = (InfCop+RecCop+SusCop)/(2*sol.sol(t0)[9][-1])
        origFreq = afs[0][0]/(2.*self.N)
  
        origPFix = self.fixProb(origFreq) 
        newPFix = self.fixProb(newFreq)   

        #print  (RecCop,SusCop,InfCop,sol.sol(t0)[9][-1]) 

        print(afs[0][0]/(2.*self.N),(InfCop+RecCop+SusCop)/(2*sol.sol(t0)[9][-1]),sol.sol(t0)[9][-1],afs[0][1],origPFix,newPFix,newPFix-origPFix,self.beta, self.tradeOff(afs[0][1]),self.sfs[int(round(afs[0][0]))-1])

        return (newPFix-origPFix)

    def fixProb(self,f):
        """
        fixation probability of a selected allele
        """
        return (1-mpmath.exp(-2*self.N*self.effects[self.resi]*f))/(1-mpmath.exp(-2*self.N*self.effects[self.resi]))

    def simSFS(self):
        return

    def sampResAllele(self):
        return

    def runAll(self):
        return

