import pylab, optparse, os
import pickle

#Method that runs the advection1D eulerSolve program in C It makes a plot of
# the estimated solutions 
#----------------------------------------------------------------------------


#Set up initial Conditions
#------------------------------------------------------

def InitialRhoWave(N_points, drhox, outfile):
    x0 = -3.0
    x1 = 5.0
    dx = (x1 - x0)/N_points
    x  = pylab.arange(x0, x1, dx)

    values = list()
    #Small pressure difference
    #----------------------------------
    
    for i in range(N_points):
        if(i == N_points/2):
            Rhox = 4.0/3.0 + drhox
        else:
            Rhox = 4.0/3.0
    
        Ex = 4.0/3.0
        pxx = 0.0
        valuesAtx = (x[i],Rhox,Ex,pxx)
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
    #  run c program with mpi
    #-----------------------------
    #with 2 processos
    os.system("mpiexec -np 8 ./eulerSolve")
    #os.system("./eulerSolve")


def plotEulerAllVals(initial2nd, final2nd,
                     final2ndprim):
    
    # Plot Results
    #-------------------------------
     


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
    pylab.title("Initial density")
    pylab.plot(x2i,rho2i)
   
    pylab.subplot(2,3,2)
    pylab.title("Initial energy density")
    pylab.plot(x2i,E2i)
    
    pylab.subplot(2,3,3)
    pylab.title("Initial px")
    pylab.plot(x2i,px2i)
    

    pylab.subplot(2,3,4)
    pylab.title("Final Rho")
    pylab.plot(x2f,rho2f, label = '2nd Order')
    
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
    
    pylab.subplot(2,3,6)
    pylab.title("Final Velocity")
    pylab.plot(x2f,vx2f)
    

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

def plotEulerPre(name):

    
    infile_name = name
    
    infile = open(infile_name, 'r')
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    pre     = [float(v[1]) for v in strvals]
    e       = [float(v[2]) for v in strvals]
    vx      = [float(v[3]) for v in strvals]

    
    
    #pylab.subplot(2,3,4)
    
    pylab.plot(x,pre)
    pylab.xlabel("x")
    pylab.ylabel("pressure")
    
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
    
def runEulerRho6times(dt, name):
    changeParamEntry("InitialType",2)
    changeParamEntry("runTime",dt)
    pylab.figure()
    pylab.suptitle("Rho vs t")
    for i in range(6):
        runEuler()
        pylab.subplot(2,3,i+1)
        plotEulerRho(name)
        
def runEuler6times(dt, name):
    changeParamEntry("InitialType",2)
    changeParamEntry("runTime",dt)
    pylab.figure()
    pylab.suptitle("Pressure vs t")
    for i in range(6):
        runEuler()
        pylab.subplot(2,3,i+1)
        plotEulerPre(name)

        
N_points = 480;

changeParamEntry("Nx", N_points)
changeParamEntry("runTime",.3)
changeParamEntry("InitialType",1)  
runEuler()
plotEulerAllVals("run/initialIn2nd.cfg",
                 "run/testFinalCons2nd.cfg","run/testFinalPrim2nd.cfg")
pylab.show()

#dt = .2
#runEulerRho6times(dt, "run/testFinalCons2nd.cfg")
#pylab.show()

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
