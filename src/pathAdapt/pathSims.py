import sys
import os
import numpy as np
import re
import math
import mpmath
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.special import lambertw
from scipy.stats import gamma
from collections import defaultdict

# to do: fix the beta in homozygous ind, which is reduced too much even when the tradeoff is very weak

"""
A class for simulating density-dependent epidemics and their effect on differentiation rate 

The idea is to:
    1.  assume a population at equilibirum
    2.  sample a resistance allele from the deleterious allele spectrum
    3.  suppose this allele is now subject to mortality/transmission/recovery tradeoff
    4.  compute change in fixation rate given outcome of epidemic
"""

class SimulateEpi():

    def __init__(self,N=1000,I0=1,I1=0,I2=0,R0=0,R1=0,R2=0,S0=998,S1=1,beta=0.1,mu=0.001,rec=0.003,L=1,al=0.184,be=0.000402,mutRate=1e-8):
        self.mutRate = mutRate
        self.mutRatesGam = {}
        self.relNumMuts = {}
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
        self.gridSize = 0.01
        self.theta = 4*N*mutRate
        self.useLimFreqs = True
        small_freqs = np.arange(1./(2.*self.N),0.01,1./(2.*self.N))
        big_freqs = np.arange(0.01,1,0.01)
        self.freqs = np.concatenate([small_freqs,big_freqs])

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
        freqs = []
        if self.useLimFreqs:
            return dFunc([i/(NN+0.) for i in self.freqs])
        return dFunc([i/(NN+0.) for i in range(1,NN)])
   
    # FULL sfs over all frequencies, for single s
    def DiscSFSSelNeg(self,s):
        NN = int(round(2*self.N))
        dFunc = np.vectorize(self.kimuraSFS)
        a = []
        if self.useLimFreqs:
            a = dFunc(s,self.freqs)
        else:
            a = dFunc(s,[i/(NN+0.) for i in range(1,NN)])
        b = np.sum(a)
        return np.divide(a,b)

    def DiscSFSSelNegNotNorm(self,s,N):
        NN = 2*N
        dFunc = np.vectorize(self.kimuraSFS)
        if self.useLimFreqs:
            a= dFunc(s,self.freqs)
        else:
            a = dFunc(s,[i/(NN+0.) for i in range(1,NN)])
        return a

    def mutRatesGamma(self,sMin,sMax,fac):
        # everything relative to human anc pop size
        prev = gamma.cdf(sMin,a=self.al,scale=1/(2.*10000*self.be))       
        mutRatesGam = {}
        mutRatesGam[sMin] = prev
        
        s = sMin*fac
        while s < sMax:
            nextVal = gamma.cdf(s,a=self.al,scale=1/(2.*10000*self.be))
            mutRatesGam[round(s,15)] = nextVal - prev
            prev = nextVal
            s *= fac
        self.mutRatesGam = mutRatesGam

    def relNumMutsGamma(self,sMin,sMax,fac):
        allVars = 0
        s = sMin
        while s < sMax:
            SFS = self.DiscSFSSelNegNotNorm(-s,self.N)
            tot = np.sum(SFS)
            self.relNumMuts[abs(round(s,15))] = tot
            allVars += tot
            s *= fac
        for val in self.relNumMuts:
            self.relNumMuts[val] /= tot
 
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
  
    # simulate an epidemic, given the frequencies of the res alleles
    def determSim(self,p=0.5,s=-0.001,sHalf=0.00005,scale=1):
 
        rHom = self.tradeOffRidge(s,sHalf,scale)
        rHet = self.rec #self.tradeOffRidge(s,sHalf,scale)
        muHom = (self.mu+self.rec)-rHom
        muHet = (self.mu+self.rec)-rHet

        betaVec = [self.beta, self.beta, self.beta]
        rVec = [self.rec,rHet,rHom]
        muVec = [self.mu,muHet,muHom]        
        
        # initialize infected
        self.I0 = 1
        self.I1 = 0
        self.I2 = 0

        # need to consider how to account for diploid state here -- 3 classes for 0,1,2 copies of derived allele
        self.S2 = int(math.floor(p**2*self.N))
        self.S1 = int(math.ceil(2*p*(1-p)*self.N))-self.I1-self.I0 
        self.S0 = max(0,self.N-(self.S1+self.I0+self.I1+self.S2+self.I2))

        # differential eqns
        def SIR(t, z, mu, r, beta):
            I0, I1, I2, S0, S1, S2, R0, R1, R2, N = z
            return [-r[0]*I0-mu[0]*I0+S0*beta[0]*(I0+I1+I2),-r[1]*I1-mu[1]*I1+S1*beta[1]*(I0+I1+I2), -r[2]*I2-mu[2]*I2+S2*beta[2]*(I0+I1+I2),
                     -S0*beta[0]*(I0+I1+I2), -S1*beta[1]*(I0+I1+I2), -S2*beta[2]*(I0+I1+I2),
                      r[0]*I0, r[1]*I1, r[2]*I2,
                      -mu[0]*I0-mu[1]*I1-mu[2]*I2]
     
        radius = 1
        tmax = 5
        while radius > 1e-5:
            tmax *= 2
            sol = solve_ivp(SIR, [0, tmax], [self.I0, self.I1, self.I2, self.S0, self.S1, self.S2, self.R0, self.R1, self.R2, self.N], 
                            args=(muVec, rVec, betaVec), dense_output=True)
            t0 = np.linspace(0, tmax, 1000)
            radius = sum(np.absolute(SIR(tmax,np.transpose(sol.sol(t0))[-1],muVec, rVec, betaVec)))

        names = {}
        names[0] = "I0"
        names[1] = "I1"
        names[2] = "I2"
        names[3] = "S0"
        names[4] = "S1"
        names[5] = "S2"
        names[6] = "R0"
        names[7] = "R1"
        names[8] = "R2"
        names[9] = "N"
   
        for i in range(len(sol.sol(t0)[:-1])):
            j = 0
            for t in range(len(sol.sol(t0)[i])):
                print(sol.sol(t0)[i][t],t0[j], names[i])
                j += 1        
        return        

    def RStarFinal(self,beta,mu,r):
        lamb = beta[0]/(r[0]+mu[0]) 
        return np.real(self.N + lambertw(-self.N*lamb*math.exp(-self.N*lamb))/lamb)

    def SFinal(self,S,beta,mu,r,i):
        lamb = beta[0]/(r[0]+mu[0]) 
        return self.N*S[i]*math.exp(-lamb*self.RStarFinal(beta,mu,r))

    def RFinal(self,S,beta,mu,r,i):
        return (self.N*S[i]-self.SFinal(S, beta, mu, r, i))*r[i]/(r[i]+mu[i]) 

    def Nfinal(self,S,beta,mu,r):
        Nf = 0
        for i in range(0,3):
            Nf += self.SFinal(S,beta,mu,r,i)+self.RFinal(S,beta,mu,r,i)
        return Nf

    def dipGenoFreqFinal(self,S,beta,mu,r,i):
        return (self.SFinal(S,beta,mu,r,i) + self.RFinal(S,beta,mu,r,i))/self.Nfinal(S,beta,mu,r)

    def deltaFixProb(self,s,xi,xf):
        """
        returns the change in fixation probability for an allele that goes from xi to xf
        """
        if abs(s*self.N) > 50:
            return 0.
        return (math.exp(-4*self.N*s*xi)-math.exp(-4*self.N*s*xf))/(1-math.exp(-4*self.N*s))
 
    def tradeOffRidge(self,s,sHalf,scale):
        """ 
        returns the new rec for Homozygous ind given the proportion of (rec+mu) due to rec
        """
        propRec = 1.-(1/(1+math.exp(scale*(math.log(abs(s))-math.log(sHalf)))))
        return (propRec*(self.mu/(self.mu+self.rec))+(self.rec/(self.mu+self.rec)))*(self.mu+self.rec)

    def totSites(self, s, N):
        standingSFS = self.DiscSFSSelNegNotNorm(s,N)
        tot = 0
        for thing in standingSFS:
            tot += thing
        return (thing*self.theta*self.L)
        

    def fixOverAllFreq(self,s,sHalf,scale):
        """
        sum over all frequencies to get the increase in fixation rate per generation with an epidemic
        1. get frequency spectrum
        2. get beta for each genotype
        3. given beta get r,mu for each genotype
        4. for each frequency, compute the change in frequency given params
        5. for each frequency, compute change in prob fix given s, mu
        6. sum over change in prob fix * SFS for each freq
        7. return
        """
        # making sure s i neg
        s = -1.*abs(s)

        # get SFS and rec for tradeoff
        normSFS = self.DiscSFSSelNeg(s)  
        SFS = self.DiscSFSSelNegNotNorm(s,self.N)  
        rHom = self.tradeOffRidge(s,sHalf,scale)
        #rHet = self.rec - (self.rec - rHom)/2
        rHet = self.rec #self.tradeOffRidge(s,sHalf,scale)
        muHom = (self.mu+self.rec)-rHom      
        muHet = (self.mu+self.rec)-rHet      

        # get other betas and mu values
        lamb = self.beta/(self.rec+self.mu)

        beta = [self.beta,self.beta,self.beta]
        r = [rHom,rHet,self.rec]
        mu = [muHom,muHet,self.mu]

        #get frequency change
        fNew = np.zeros(2*self.N-1)
        pChange = np.zeros(2*self.N-1)
        pChangeWeightSFS = np.zeros(2*self.N-1)

        # prob of sampling selection coeff s
        sProbs = []
        sProbLoc = {}
        i = 0
        val = 0
        tot = 0
        for sVal in self.mutRatesGam:
            val = self.mutRatesGam[sVal]*self.relNumMuts[sVal]
            tot += val
            sProbs.append(val)
            sProbLoc[sVal] = i
            i += 1

        for j in range(len(sProbs)):
            sProbs[j]/=tot
       
        # get freqs to loop over
        small_freqs = np.arange(1./(2.*self.N),0.01,1./(2.*self.N))
        big_freqs = np.arange(0.01,1,0.01)
        freqs = np.concatenate([small_freqs,big_freqs])
        j = 0
        # main calc loop
        for i in freqs:
            # first compute the epi model
            S0 = i**2
            S1 = 2*i*(1-i)
            S2 = (1-i)**2
            S = [S0,S1,S2]
            Sf = self.dipGenoFreqFinal(S,beta,mu,r,0)
            Xf = Sf**0.5
            fNew[j] = Xf

            #next compute the change in the fix prob weighted by the number of alleles at this freq
            pChange[j] = self.deltaFixProb(s,i,Xf)
            pChangeWeightSFS[j] = pChange[j]*normSFS[j]
            
            j += 1        

        totAdd = np.sum(pChangeWeightSFS)

        #relRate = (totAdd+totOrig)/(totOrig)

        return [sProbs[sProbLoc[abs(round(s,15))]]*totAdd,self.mutRatesGam[abs(round(s,15))]*self.fixProb(s,1./(self.N*2.))*self.mutRate*2*self.N*self.L,sProbs[sProbLoc[abs(round(s,15))]]]

    def fixProb(self,s,f):
        """
        fixation probability of a selected allele
        """
        return (1-mpmath.exp(-4*self.N*s*f))/(1-mpmath.exp(-4*self.N*s))
 
    def fixProbBack(self):
  
        tot = 0
        for val in self.mutRatesGam:
            tot += self.mutRatesGam[val]*2*self.N*self.mutRate*self.L*self.fixProb(-val,1./(2*self.N))
        return tot

    def stochSim(self, pEpi, nGen, sMin, sMax, fac, sHalf, scale):

        self.mutRatesGamma(sMin,sMax,fac) 
        self.relNumMutsGamma(sMin,sMax,fac) 

        # make array that samples over the prob of getting a particular s
        sProbs = []
        sProbLoc = {}
        i = 0
        val = 0
        tot = 0
        for s in self.mutRatesGam:
            val = self.mutRatesGam[s]*self.relNumMuts[s]
            tot += val
            sProbs.append(val)
            sProbLoc[i] = s
            i += 1      
        
        for j in range(len(sProbs)):
            sProbs[j]/=tot

        # get background fix per gen
        rate = self.fixProbBack()

        # check how many epidemics happen this gen
        def getNumEpi(p):
            return np.random.poisson(p)        

        # get freq and s for effect allele
        def getSAndX():
            index = np.random.choice(len(sProbs), p=sProbs)
            s = sProbLoc[index]
            SFS = self.DiscSFSSelNeg(-s)
            index = np.random.choice(len(SFS), p=SFS)
            x = (index+1)/(2*self.N)           
          
            return -1*abs(s),x   
        
        # main loop
        totBack = 0
        totAdd = 0
        totMaxAdd = 0
        for gen in range(0,nGen):
            numEpi = getNumEpi(pEpi)
            while numEpi > 0:
                s,x = getSAndX()           
             
                # get params
                rHom = self.tradeOffRidge(s,sHalf,scale)
                #rHet = self.rec - (self.rec - rHom)/2
                rHet = self.rec #self.tradeOffRidge(s,sHalf,scale)
                muHom = (self.mu+self.rec)-rHom
                muHet = (self.mu+self.rec)-rHet

                # get other betas and mu values
                lamb = self.beta/(self.rec+self.mu)

                beta = [self.beta,self.beta,self.beta]
                r = [rHom,rHet,self.rec]
                mu = [muHom,muHet,self.mu]

                S0 = x**2
                S1 = 2*x*(1-x)
                S2 = (1-x)**2
                S = [S0,S1,S2]
                Sf = self.dipGenoFreqFinal(S,beta,mu,r,0)
                Xf = Sf**0.5
                pfixDiff = self.fixProb(s,Xf)-self.fixProb(s,x) 
                if np.random.random() < pfixDiff:
                    totAdd += 1
                totMaxAdd += 1-self.fixProb(s,x)
                numEpi -= 1
            totBack += rate

        totBackRand = np.random.poisson( totBack) 

        return [totAdd/(totAdd+totBackRand),totMaxAdd/(totMaxAdd+totBackRand),totAdd, totBackRand]

    def simTraj(self,s,x,maxF):
        x0 = x
        traj = [x]
        while x > 0 and x < maxF:
            next_gen = np.random.binomial(2*self.N,p=(1+s)*x)
            x = next_gen/(2.*self.N)
            if x == 0.:
                traj = []
                x = x0 
            traj.append(x)
        return traj

    def simTrajFinal(self,s,x):
        x0 = x
        traj = [x]
        while x > 0 and x < 1:
            next_gen = np.random.binomial(2*self.N,p=(1+s)*x)
            x = next_gen/(2.*self.N)
            traj.append(x)
        return traj

    def genSingleEvoEpi(self,s,f,sHalf,scale):
       
        traj = self.simTraj(s,f,0.5)
        
        x = traj[-1]
        S0 = x**2
        S1 = 2*x*(1-x)
        S2 = (1-x)**2
        S = [S0,S1,S2]
        
        rHom = self.tradeOffRidge(s,sHalf,scale)
        #rHet = self.rec - (self.rec - rHom)/2
        rHet = self.rec #self.tradeOffRidge(s,sHalf,scale)
        muHom = (self.mu+self.rec)-rHom
        muHet = (self.mu+self.rec)-rHet

        # get other betas and mu values
        lamb = self.beta/(self.rec+self.mu)

        beta = [self.beta,self.beta,self.beta]
        r = [rHom,rHet,self.rec]
        mu = [muHom,muHet,self.mu]
        
        Sf = self.dipGenoFreqFinal(S,beta,mu,r,0)
        Xf = Sf**0.5

        print("#"+str(len(traj)))

        traj.extend(self.simTrajFinal(s,Xf))

        return traj
