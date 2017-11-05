import pylab, optparse, os


#Method that runs the advection1D eulerSolve program in C It makes a plot of
# the estimated solutions
#----------------------------------------------------------------------------

def runEuler():

    #  run c program
    #-----------------------------
    os.system("./eulers")

    # Plot Results
    #-------------------------------
    #Initial Conditions
    infile_name = "testInitial.cfg"
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    rho     = [float(v[1]) for v in strvals]
    E       = [float(v[2]) for v in strvals]
    px      = [float(v[3]) for v in strvals]


    pylab.subplot(2,3,1)
    pylab.title("Initial Rho")
    pylab.plot(x,rho)
    pylab.subplot(2,3,2)
    pylab.title("Initial E")
    pylab.plot(x,E)
    pylab.subplot(2,3,3)
    pylab.title("Initial px")
    pylab.plot(x,px)

    #Final Conditions
    infile_name = "testFinalCons.cfg"
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    rho     = [float(v[1]) for v in strvals]
    E       = [float(v[2]) for v in strvals]
    px      = [float(v[3]) for v in strvals]

    #U = list(rho,E,px);

    pylab.subplot(2,3,4)
    pylab.title("Final Rho")
    pylab.plot(x,rho)
    #pylab.subplot(2,3,5)
    #pylab.title("Final E")
    #pylab.plot(x,E)
    #pylab.subplot(2,3,6)
    #pylab.title("Final px")
    #pylab.plot(x,px)

    infile_name = "testFinalPrim.cfg"
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    pre     = [float(v[1]) for v in strvals]
    e       = [float(v[2]) for v in strvals]
    vx      = [float(v[3]) for v in strvals]

    pylab.subplot(2,3,5)
    pylab.title("Final Pressure")
    pylab.plot(x,pre)
    pylab.subplot(2,3,6)
    pylab.title("Final Velocity")
    pylab.plot(x,vx)

runEuler()


'''def advect(initial, initialName, method, methodName, Nx):
    #Initialize parmeter file
    #--------------------------------------------------------
    filename = "param.cfg"
    pfile = open(filename, 'w')
    pfile.write("Nx = %d \n" % Nx)
    pfile.close()

    #  run c program
    #-----------------------------
    os.system("./advect2_dbr250")


    #Plot results
    #-------------------------------------------------------
    pylab.title("Advection of a %s wave using the %s method \n The L1_norm is calculated to be: %f" % (initialName ,methodName, norm))

    plotTitles = ["Initial Condition","Analytical Solution", "Estimated Solution"]
    plotFiles = ["advect_ini.dat", "advect_ana.dat", "advect_fin.dat"]
    for infile_title,infile_name in zip(plotTitles,plotFiles):
        infile = open(infile_name)
        lines = infile.readlines()

        strvals = [l.strip().split() for l in lines]
        x       = [float(v[0]) for v in strvals]
        y       = [float(v[1]) for v in strvals]


        pylab.plot(x,y, label = infile_title)

        pylab.xlabel("x (with Nx = %d )" % Nx )
        pylab.ylabel("u")

        pylab.legend()
    #----------------------------------------------------------
'''


pylab.show()
