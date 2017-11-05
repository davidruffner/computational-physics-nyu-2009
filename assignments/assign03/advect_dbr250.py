#!/usr/bin/env python 

"""
Author: David Ruffner, NYU
Date: 9/26/09


This program integrates the advection equation using the Lax-Friedrichs Method.
"""

import pylab
import math

#The advection equation can be written as Ut + a*Ux = 0, where a is the rate of 
# advection, Ut is the derivative of U wrt time and Ux is the derivative of 
# U wrt x. It can be approximated by the Lax Friedrichs method which is given by
# Uj_n+1 = .5*(Uj-1_n + Uj+1_n) - a*k*(Uj+1_n-Uj-1_n)/h
# where a is the advection rate, k is the time step, h is the x step,j is the x index, and n is the time index.

# It works by discretizing the space, and then iteritavely finding the next 
#values of U.

#Intial functions we will consider
#----------------------------------------------------------------------------
def gaussian(x,a):
    return pylab.exp(-a*x**2)

def sin(x,w):
    return pylab.sin(w*x)

def square(x,b):
    if x< -(b/10.) or x>= (b/10):
        return 0
    else:
        return 1

def advection(f,name):
    L1_values=list()
    dx_values = list()
    n_dx = 6# exponent of the dx step
    for i in range(n_dx):
        
        #Discretizing Space
        #---------------------------------------------------------------------------
        A = 1.  #advection speed

        a = -1. # min x
        b = 1.  # max x
        nx = 100*2**i # number of x steps
        h = (b-a)/nx # the step size in x
        x, dx   =  pylab.linspace(a, b, nx+1, retstep = True)
        dx_values.append(dx)

        t_f = 2/A
        t_i = 0.
        k = abs(.5*h/A)

        n_graphs = 11# rough number (minus 1)of graphs of U that will be displayed

        #Intial Values for U
    #---------------------------------------------------------------------------


        U = [f(xi,pylab.pi) for xi in x]
        U_exact = [f(xi,pylab.pi) for xi in x]
        #U_un = pylab.array([gaussian(xi,pylab.pi) for xi in x])
        if i == n_dx-1:
            pylab.figure()
            pylab.subplot('211')
            pylab.xlabel("x")
            pylab.ylabel("U")
            pylab.title("U vs x at various times for initial condition: %s" % name)

        #Integration
        #---------------------------------------------------------------------------
        t = t_i
        iter = 0

        while t<= t_f:
            Uim1 = pylab.roll(U,1)
            Uip1 = pylab.roll(U,-1)
            U = .5*(Uim1 + Uip1) - .5*k*(Uip1-Uim1)/dx
    
            #Uim1_un = pylab.roll(U_un,1)
            #Uip1_un = pylab.roll(U_un,-1)
            #U_un = U_un - .5*k*(Uip1_un-Uim1_un)/dx
    
            t+=k
            iter+=1
            if iter % (2*nx/n_graphs) == 0 and i==n_dx-1:
                pylab.plot(x,U)
                #pylab.plot(x,U_un)

        L1_norm = dx*sum(abs(U - U_exact))
        L1_values.append(L1_norm)
    
    pylab.subplot('212')
    pylab.loglog(dx_values,L1_values,'-o')
    pylab.xlabel("dx")
    pylab.ylabel("L1_norm")
    

    slope = (pylab.log(L1_values[n_dx-1])-pylab.log(L1_values[0]))/(pylab.log(dx_values[n_dx-1])-pylab.log(dx_values[0]))
    print name
    print "The slope of convergence is : %f" % slope
  


functions = [sin, gaussian, square]
names = ["sin","gaussian","square"]


for j in range(0,3):
    advection(functions[j],names[j])
pylab.show()
