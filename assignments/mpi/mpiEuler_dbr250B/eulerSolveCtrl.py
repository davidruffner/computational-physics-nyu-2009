import pylab, optparse, os
import pickle

#Method that runs the advection1D eulerSolve program in C It makes a plot of
# the estimated solutions 
#----------------------------------------------------------------------------


#Set up initial Conditions
#------------------------------------------------------

def InitialPwave(N_points, dPx, outfile):
    x0 = -3.0
    x1 = 5.0
    dx = (x1 - x0)/N_points
    x  = pylab.arange(x0, x1, dx)

    values = list()
    #Small pressure difference
    #----------------------------------
    
    for i in range(N_points):
        if(i == N_points/2):
            Px = 2.0 + dPx
        else:
            Px = 2.0
    
        rhox = 1.0
        vxx = 0.0
        valuesAtx = (x[i],rhox,Px,vxx)
        values.append(valuesAtx)
        

    filename = outfile
    pfile = open(filename, 'w')

    for valueAtx in values:
        for comp in valueAtx:
            comp = str(comp)
        pfile.write("%s %s %s %s\n" % (valueAtx[0],valueAtx[1],valueAtx[2],
                                   valueAtx[3]) )
        
    pfile.close()

#  run c program
#-----------------------------
def runEuler():
    #  run c program
    #-----------------------------
    os.system("./eulerSolve")



def plotEulerAllVals(initial1st, initial2nd, final1st, final1stprim, final2nd,
                     final2ndprim):
    
    # Plot Results
    #-------------------------------
    #Initial Conditions   
    infile_name = initial1st
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x1i       = [float(v[0]) for v in strvals]
    rho1i     = [float(v[1]) for v in strvals]
    E1i       = [float(v[2]) for v in strvals]
    px1i      = [float(v[3]) for v in strvals]

        
    #Final Conditions
    infile_name = final1st
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x1f       = [float(v[0]) for v in strvals]
    rho1f     = [float(v[1]) for v in strvals]
    E1f       = [float(v[2]) for v in strvals]
    px1f      = [float(v[3]) for v in strvals]
        

    infile_name = final1stprim
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x1f       = [float(v[0]) for v in strvals]
    pre1f     = [float(v[1]) for v in strvals]
    e1f       = [float(v[2]) for v in strvals]
    vx1f      = [float(v[3]) for v in strvals]

    


    #Initial Conditions 2nd Order
        
    infile_name = initial2nd
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x2i       = [float(v[0]) for v in strvals]
    rho2i     = [float(v[1]) for v in strvals]
    E2i       = [float(v[2]) for v in strvals]
    px2i      = [float(v[3]) for v in strvals]
  
   
   
    #Final Conditions 2nd Order

    
    infile_name = final2nd
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x2f       = [float(v[0]) for v in strvals]
    rho2f     = [float(v[1]) for v in strvals]
    E2f       = [float(v[2]) for v in strvals]
    px2f      = [float(v[3]) for v in strvals]

    
    infile_name = final2ndprim
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x2f       = [float(v[0]) for v in strvals]
    pre2f     = [float(v[1]) for v in strvals]
    e2f       = [float(v[2]) for v in strvals]
    vx2f      = [float(v[3]) for v in strvals]

    #plotting
    pylab.figure()

    pylab.subplot(2,3,1)
    pylab.title("Initial Rho")
    pylab.plot(x2i,rho2i)
    pylab.plot(x1i,rho1i)
    pylab.subplot(2,3,2)
    pylab.title("Initial E")
    pylab.plot(x2i,E2i)
    pylab.plot(x1i,E1i)
    pylab.subplot(2,3,3)
    pylab.title("Initial px")
    pylab.plot(x2i,px2i)
    pylab.plot(x1i,px1i)

    pylab.subplot(2,3,4)
    pylab.title("Final Rho")
    pylab.plot(x2f,rho2f, label = '2nd Order')
    pylab.plot(x1f,rho1f,  label = '1st Order')
    pylab.legend()
    #pylab.subplot(2,3,5)
    #pylab.title("Final E")
    #pylab.plot(x,E)
    #pylab.subplot(2,3,6)
    #pylab.title("Final px")
    #pylab.plot(x,px)

    
    pylab.subplot(2,3,5)
    pylab.title("Final Pressure")
    pylab.plot(x2f,pre2f)
    pylab.plot(x1f,pre1f)
    pylab.subplot(2,3,6)
    pylab.title("Final Velocity")
    pylab.plot(x2f,vx2f)
    pylab.plot(x1f,vx1f)

def plotEulerRho(name):

    
    infile_name = name
    
    infile = open(infile_name, 'r')
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    rho     = [float(v[1]) for v in strvals]
    E       = [float(v[2]) for v in strvals]
    px      = [float(v[3]) for v in strvals]

    
    
    #pylab.subplot(2,3,4)
    
    pylab.plot(x,rho)
    pylab.xlabel("x")
    pylab.ylabel("rho")
    
def changeParamEntry(variable,value):
    #Initialize parmeter file
    # variable is the name of the variable in quotes
    # and values is the value of the variable
    #--------------------------------------------------------
    filename = "param.cfg"
    pfile = open(filename, 'r')

    
    #Time step
    
    
    lines = pfile.readlines()
    strvals = [l.strip().split() for l in lines]

    inputNameList = list();
    inputValuList = list();
    
    for v in strvals:
        inputNameList.append(v[0]);
        if(v[0] == variable):
            inputValuList.append(value)
        else:
            inputValuList.append(v[2]);

    
    
    pfile.close()

    
    pfile = open(filename, 'w')

    for infile_name,infile_value in zip(inputNameList , inputValuList):
        
        infile_value = str(infile_value)
        pfile.write("%s = %s \n" % (infile_name, infile_value) )
        
    pfile.close()

def runEulerNtimes(N, name):
    changeParamEntry("InitialType",2)
    changeParamEntry("runTime",.1)
    for i in range(N):
        runEuler()
        pylab.figure()
        plotEulerRho(name)
    
def runEuler6times(dt, name):
    changeParamEntry("InitialType",2)
    changeParamEntry("runTime",dt)
    pylab.figure()
    pylab.suptitle("Rho vs t")
    for i in range(6):
        runEuler()
        pylab.subplot(2,3,i+1)
        plotEulerRho(name)
N_points = 50;
dEx = .01;
InitialPwave(N_points+2, dEx, "run/initialIn.cfg")
InitialPwave(N_points+4, dEx, "run/initialIn2.cfg")

changeParamEntry("Nx", N_points)
changeParamEntry("runTime",.1)
changeParamEntry("InitialType",3)  
runEuler()
plotEulerAllVals("run/initialIn.cfg","run/initialIn2.cfg",
                 "run/testFinalCons.cfg","run/testFinalPrim.cfg",
                 "run/testFinalCons2nd.cfg","run/testFinalPrim2nd.cfg")
pylab.show()
#changeParamEntry("Nx", 50)
#changeParamEntry("runTime",.1)
#changeParamEntry("InitialType",1)  
#runEuler()
#plotEulerAllVals()
#name = "run/testFinalCons2nd.cfg"

#dt = .2
#runEuler6times(dt, name)

  
        
#pylab.show()

#put in the capability to output the time elapsed each time the
#c code runs. put right into the param.cfg file.
