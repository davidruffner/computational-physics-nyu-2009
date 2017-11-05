#!/usr/bin/env python

import pylab

a=2.0
b=1.0
c=-3.0

def myfunc(x):
    return a*x**2 + b*x + c

def myFunc(x):
    return a*x**3/3 + b*x**2/2 + c*x

def L1_norm(f,g,dx):    
    return sum(abs(f-g))*dx
    

x0 = -11.0
x1 = 10.0

N_samples = 4


samples, dx_sample = pylab.linspace(x0,x1,N_samples,retstep=True)

values = [myfunc(x) for x in samples]

integral = pylab.zeros(N_samples)

true_int = [myFunc(x) for x in samples]

integral[0] = myFunc(x0)


#  for N_samples in pylab.logspace(1,4,10):


#    N_samples=int(N_samples)


for i in range(1,N_samples):

    dx = samples[i]-samples[i-1]

    x = (samples[i] + samples[i-1])/2
    
    integral[i] = integral[i-1] + myfunc(x) *dx

print L1_norm(integral,true_int,dx_sample)

pylab.plot(samples,values,label=r"$f(\lambda)$")

pylab.plot(samples, integral,'-o',label="F(x)")

pylab.plot(samples,true_int,'-s',label="true_int")

pylab.legend()

pylab.show()
