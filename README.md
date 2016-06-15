# isingmod2d

A module for Metropolis-Hastings computations on the Ising model using helix-type boundaries.

To compile and install, simply run 
>python setup.py build_ext --inplace

Typical usage:

>from isingmod import *
>
>lab = ising2d(20,1.0) #a 20x20 lattice with kt= 1.0
>lab.kt = 3.0 #change temperature
>lab.advance3() #advance solution one step using helix-type boundaries
>
>m,s,e,c,q = lab.get_stats(neq = 30000, ns = 50000) #equalize for 30000 iterations, sample for 50000


m = magnetization
s = susceptibility
e = energy
c = heat capacity
q = binder ratio



