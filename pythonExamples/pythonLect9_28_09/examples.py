#!/usr/din/env

import pylab

def demo1():
    for i in range(5):
        x = pylab.linspace(-1,1,100)
        y = [pylab.exp(-5*xi**2/i**2) for xi in x]

    pylab.plot(x,y)

pylab.show()

def demo2():
    i = 1
    x = pylab.linspace(-1,1,100)
    y = [pylab.exp(-5*xi**2/i**2) for xi in x]
    pylab.subplot('221')
    pylab.plot(x,y)

    i = 2
    x = pylab.linspace(-1,1,100)
    y = [pylab.exp(-5*xi**2/i**2) for xi in x]
    pylab.subplot('222')
    pylab.plot(x,y)

    i = 3
    x = pylab.linspace(-1,1,100)
    y = [pylab.exp(-5*xi**2/i**2) for xi in x]
    pylab.subplot('223')
    pylab.plot(x,y)

    i = 4
    x = pylab.linspace(-1,1,100)
    y = [pylab.exp(-5*xi**2/i**2) for xi in x]
    pylab.subplot('224')
    pylab.plot(x,y)
    
    pylab.show()

demo2()

#numpy examples
def demo4():

    Nx = 100
    Ny = 100
    mu = .4

    z = pylab.zeros([Nx,Ny])

    x = pylab.linspace(-1,1,Nx)
    y = pylab.linspace(-1,1,Ny)

    for i in range(100):
        for j in range(100):
            z[i,j]  = pylab.exp((x[i]**2 + y[i]**2)/mu) #works for arrays, not for lists
    pylab.imshow(z)
    pylab.show()
demo4()



