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


def L1_norm(f,g,dx):
    return dx*sum(abs(pylab.array(f)-pylab.array(g)))


def step_lax_fried(U,lam):
    Uim1 = pylab.roll(U,1)
    Uip1 = pylab.roll(U,-1)
    return .5*(Uim1 + Uip1) - .5*(lam)*(Uip1-Uim1)

def step_lax_wendroff(U,lam):
    Uip1 = pylab.roll(U, -1)
    Uim1 = pylab.roll(U, +1)
    return U - 0.5*A*lam * (Uip1 - Uim1) + 0.5*(A*lam)**2 * (Uip1 - 2*U + Uim1)


def advect(initial,Func_name,width,nx):#(initial waveform,narrowness of waveform,number of steps)
    #Will calculate one period of advection using lax Friedrichs and Lax Wendroff m
    #methods, and will plot the resulting U's vs the U_analytic
 
    #parameters for run
    #---------------------------------------------------------------------------
    a = -1. # min x
    b = 1.  # max x
    
    h = (b-a)/nx # the step size in x
    x, dx   =  pylab.linspace(a, b, nx+1, retstep = True)

    t_f = 2/A
    t_i = 0.
    k = abs(.5*h/A)



    steps = [step_lax_fried, step_lax_wendroff]
    names = ["Lax_friedrich","Lax_Wendroff"]
    #actual integration
    #----------------------------------------------------------------------------
    for step,name in zip(steps,names):
        t = t_i
        U = pylab.array([initial(xi,width) for xi in x])
        while t<= t_f:
            U = step(U,k/dx)
            t+= k
        pylab.title("U vs x, with w=%s" % width)
        pylab.plot(x,U,'-o', label = name)
        U_anal = [initial(xi,width) for xi in x]
    
    pylab.plot(x,U_anal, label = "analytical")
    pylab.xlabel("x (with nx=%s)" % nx)
    pylab.ylabel("y")
    pylab.legend()

def advectConv(initial, Func_name, width, order):
    

    #parameters for run
    #---------------------------------------------------------------------------
    a = -1. # min x
    b = 1.  # max x
    steps = [step_lax_fried, step_lax_wendroff]
    names = ["Lax_friedrich","Lax_Wendroff"]

    L1    = { }
    Nx    = { }

    for step,name in zip(steps,names):
        L1[name] = [ ]
        Nx[name] = [ ]
        print "Running method:", name

        for N_points in pylab.logspace(1,order,20):
            h = (b-a)/N_points # the step size in x
            x, dx   =  pylab.linspace(a, b, N_points+1, retstep = True)

            t_f = 2/A
            t_i = 0.
            k = abs(.5*h/A)   
    #actual integration
    #----------------------------------------------------------------------------
            t = t_i
            U = pylab.array([initial(xi,width) for xi in x])
            while t<= t_f:
                U = step(U,k/dx)
                t+= k
            U_anal = [initial(xi,width) for xi in x]

            L1[name].append(L1_norm(U, [initial(xi,width) for xi in x], dx))
            Nx[name].append(int(N_points))
    for name in L1:
        pylab.loglog(Nx[name], L1[name], '-o', label= name)
    pylab.title("Convergence rate for %s" % Func_name)
    pylab.xlabel("Number of points")
    pylab.ylabel("L1_norm")
    pylab.legend()
    
    


#------------------------------------------------------------------
A = 1.  #advection speed
#pylab.figure()
#advect(gaussian,5,500)
#pylab.figure()
#advect(gaussian,10,500)



functions = [sin, gaussian, square]
names = ["sin","gaussian","square"]


for j in range(0,3):
    pylab.subplot(2,3,j+1)
    advect(functions[j],names[j],pylab.pi,510)
    pylab.subplot(2,3,(j+3+1))
    advectConv(functions[j],names[j],pylab.pi,4)
pylab.show()
