from numpy import *

class ising2d():
    def __init__(self, N, kT, flipratio = 0.1):
        self.kT = kT
        self.N = int(N)
        self.N2 = int(N**2)
        self.Z = ones(self.N**2, dtype = int)
        self.Z[random.choice(self.N2, self.N2/2, replace = False)] = -1

        self.Z_ind = arange(self.N2)
        self.n = int(flipratio*self.N2)
        self.flipratio = flipratio
        
        self.env = self.energy_env_torus() #_torus(arange(self.N2))
        
        self.n_logdist = 1000000
        self.logdist = self.kT*log(random.uniform(0,1,self.n_logdist))/2.0
        
        
        self.t = zeros(10, dtype = float)
        
    def reset(self, kT):
        self.kT = kT
        self.Z = ones(self.N**2, dtype = int)
        self.Z[random.choice(self.N2, self.N2/2, replace = False)] = -1
        self.Z_ind = arange(self.N2)
        self.n = int(self.flipratio*self.N2)
        self.env = self.energy_env_torus()
        
    def report(self):
        print "Z:"
        print self.Z.reshape((self.N,self.N))
        #print "E2:"
        #print self.energy_env2().reshape((self.N,self.N))
        print "E:"
        print self.energy_env(arange(self.N2)).reshape((self.N,self.N))
        
        print "dE:"
        dE = self.energy_env(arange(self.N2)).reshape((self.N,self.N))*self.Z.reshape((self.N,self.N))
        print self.energy_env(arange(self.N2)).reshape((self.N,self.N))*self.Z.reshape((self.N,self.N))
        
        print "prob:"
        print exp(-2*dE/self.kT)
        
        #print self.energy_env2().reshape((self.N,self.N))==self.energy_env(arange(self.N2)).reshape((self.N,self.N))
    def energy_env2(self):
        #env = self.Z[self.N-1:-self.N-1]+self.Z[self.N+1:-self.N+1]  #+ self.Z[2:]
        return roll(self.Z,1)+roll(self.Z,-1)+roll(self.Z,self.N) + roll(self.Z,-self.N)
        
    def energy_env(self, z_ind):
        #return energy environment for positions in z_ind
        return self.Z[z_ind-1] + self.Z[(z_ind+1)%self.N2] +self.Z[(z_ind + self.N)%self.N2] + self.Z[z_ind - self.N]  #+ self.Z[(z_ind - self.N)%self.N2]

    def impose_periodicity(self, Z):
        Z[0,:] = Z[-2,:]
        Z[:,0] = Z[:,-2]
        Z[-1,:] = Z[1,:]
        Z[:,-1] = Z[:,1]
    def energy_per_site(self,Z):
        #calculate energy per site
        return (Z[2:,1:-1]+Z[:-2,1:-1] +Z[1:-1,2:]+Z[1:-1,:-2]) #[1:-1,1:-1]
    
    
    
    def vec_env_helix(self, z_ind):
        return self.Z[z_ind-1] + self.Z[(z_ind+1)%self.N2] + self.Z[(z_ind+self.N)%self.N2] + self.Z[z_ind-self.N]
    
    def vec_env_torus(self,z_ind):
        #calculate energy environment in self.Z[z_ind]
        return self.Z[z_ind-1] + self.Z[(z_ind+1)%self.N2] + self.Z[(z_ind+self.N)%self.N2] + self.Z[z_ind-self.N]
    
    #def vec_energy_env_torus(self, z_ind):
    #    z = self.Z.reshape((self.N,self.N))
    #    return (roll(z, 1,0) + roll(z,-1,0) + roll(z,1,1) + roll(z,-1,1)).ravel()    
    
    
    def energy_env_torus(self):
        z = self.Z.reshape((self.N,self.N))
        return (roll(z, 1,0) + roll(z,-1,0) + roll(z,1,1) + roll(z,-1,1)).ravel()
        

    def energy_env_helix(self):
        #z = self.Z.reshape((self.N,self.N))
        return roll(self.Z, 1) + roll(self.Z,-1) + roll(self.Z,self.N) + roll(self.Z,-self.N)


    def energy_env_torus2(self,z_ind):
        z = self.Z.reshape((self.N, self.N))
        Z = zeros((self.N+2,self.N+2))
        Z[1:-1,1:-1] = z
        self.impose_periodicity(Z)
        eps = self.energy_per_site(Z)
        
        return eps.ravel() #[z_ind]
    def log_uniform(self, n):
        return self.logdist[random.randint(0,self.n_logdist,n)]
    
    
    def advance3(self):
        sel = random.choice(self.N2, self.n, replace = False) #bottleneck
        
        
        dz = self.Z[sel]
        denv = self.env[sel]
        
        dE = dz*denv
        
        flips = -dE>self.log_uniform(self.n)
        dz[flips] *= -1            
        self.Z[sel] = dz
        
        self.update_env(2*dz[flips], sel[flips]) #bottleneck
    
        
    def update_env(self, signs, z_ind):
        #update energy environment with a change in Z[z_ind]
        
        self.env[z_ind-1] += signs
        self.env[(z_ind+1)%self.N2] += signs
        self.env[z_ind-self.N] += signs
        self.env[(z_ind+self.N)%self.N2] += signs
        
        #self.env[z_ind-1] = self.
    
    
    
    def update_env_torus(self, signs, z_ind):
        #update energy environment with a change in Z[z_ind]
        #CHANGE THIS ONE
        #(z_ind-1)%self.N + (z_ind-1)//self.N
        
        
        
        self.env[z_ind-1] += signs
        self.env[(z_ind+1)%self.N2] += signs
        self.env[z_ind-self.N] += signs
        self.env[(z_ind+self.N)%self.N2] += signs
        
        #self.env[z_ind-1] = self.   
        
    def advance2(self):
        #flps = exp(-2.0*self.Z*self.env/float(self.kT))>random.uniform(0,1,self.N2)
        flps = -2.0*self.Z*self.env/float(self.kT)>self.log_uniform(self.N2) #cheating here
        
        flps[random.choice(self.N2, self.N2-self.n, replace = False)] = False
        
        self.Z[flps] *= -1
        
        self.env = self.energy_env_torus() #self.energy_env_torus(self.Z_ind)
        
        
    def advance(self):
        #choose a subset of n spins

        sel = random.choice(self.N2, self.n, replace = False)
        dz = self.Z[sel]
        denv = self.env[sel]
        
        dE = dz*denv
        
        flips = exp(-2.0*dE/self.kT)>random.uniform(0,1,self.n)

        dz[flips] *= -1            
        self.Z[sel] = dz
        
        denv[flips] = self.energy_env(sel[flips])
        
        
        self.env[sel] = denv
    def get_stats(self, Neq, Nsamples):
        #sample frequency
        
        fs = self.N2/self.n
        
        #print fs
        
        m = zeros(Nsamples/fs, dtype = float)
        e = zeros(Nsamples/fs, dtype = float)
        
        
        #print self.n, self.N2/self.n, Nsamples, Nsamples/fs
        
        #equilibrate
        for i in range(Neq):
            self.advance3()
        
        
        
        for i in range(Nsamples/fs):
            for j in range(fs):
                self.advance3()
            m[i]= self.mag()
            e[i]= self.energy()
        #print m
        
        q = mean((m/self.N**2)**2)/mean(abs(m/self.N**2))**2
        return mean(m)/self.N2, (var(m)/self.kT)/self.N2, mean(e)/self.N2, (var(e)/self.kT**2)/self.N2, q
    def sweep(self,kt,neq = 20000, ns = 50000):
        m = zeros(len(kt))
        ms = zeros(len(kt))
        e = zeros(len(kt))
        cv = zeros(len(kt))
        q = zeros(len(kt))
        
        
        for t in range(len(kt)):
            self.reset(kt[t])
            self.logdist = self.kT*log(random.uniform(0,1,self.n_logdist))/2.0
            #self.logdist/= self.kT
            
            m[t],ms[t],e[t],cv[t], q[t] = self.get_stats(neq,ns)
            #self.logdist *= self.kT
            #M[t] = m
        return m,ms,e,cv,q
    
    def mag(self):
        return sum(self.Z)
        
    def energy(self):
        return -sum(self.Z*self.env)/2.0
        #return -sum(self.Z*(roll(self.Z,1)+roll(self.Z,self.N)) )
    