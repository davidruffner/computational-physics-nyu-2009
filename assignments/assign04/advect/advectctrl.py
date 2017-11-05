import pylab, optparse, os

#opt parse arguements, type advect_ini.dat, advect_fin.dat, and norm.dat
#---------------------------------------------------------------------------------
# Not in use now
parser = optparse.OptionParser()
opts, args = parser.parse_args()



#Method that runs the advection program in C for a certain method, initial con
#-ditions and Nx. It makes a plot of the analytical and estimated solutions
# and it also 
#----------------------------------------------------------------------------
def advect(initial, initialName, method, methodName, Nx):
    #Initialize parmeter file
    #--------------------------------------------------------
    filename = "param.cfg"
    pfile = open(filename, 'w')
    pfile.write("Nx = %d \n" % Nx)
    pfile.write("x1 = 1.0 \n")
    pfile.write("x0 = -1.0 \n")
    pfile.write("CFL = 0.5 \n")
    pfile.write("a = 1.0 \n")
    pfile.write("t_max = 2.0 \n ")
    pfile.write("init = %d \n " % initial)
    pfile.write("meth = %d \n" % method)
    pfile.close()

    #  run c program
    #-----------------------------
    os.system("./advect_dbr250")


    #recover the calculation of the norm
    #-----------------------------------
    normFile = open("norm.dat", 'r')
    lines = normFile.readlines()
    norm = float(lines[0])
    Nx = int(lines[1])
    

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
        


#Method that runs the advection program in C for a certain method, initial con
#-ditions and values of Nx. It makes a plot of the convergence rate and prints
# the slope 
#----------------------------------------------------------------------------
def advectConv(initial, initialName, method, methodName, Nxs):
    L1 = list()
    Nx = list()

    for Nxvalue in Nxs:
        
        #Initialize parmeter file
        #--------------------------------------------------------
        filename = "param.cfg"
        pfile = open(filename, 'w')
        pfile.write("Nx = %d \n" % Nxvalue)
        pfile.write("x1 = 1.0 \n")
        pfile.write("x0 = -1.0 \n")
        pfile.write("CFL = 0.5 \n")
        pfile.write("a = 1.0 \n")
        pfile.write("t_max = 2.0 \n ")
        pfile.write("init = %d \n " % initial)
        pfile.write("meth = %d \n" % method)
        pfile.close()

        #  run c program
        #-----------------------------
        os.system("./advect_dbr250")


        #recover the calculation of the norm
        #-----------------------------------
        normFile = open("norm.dat", 'r')
        lines = normFile.readlines()
        norm = float(lines[0])
        
        
        #append result onto list of norms and nx values
        #---------------------------------
        L1.append(norm)
        Nx.append(Nxvalue)


    #Calculate slope of convergence
    #---------------------------------------------------------------
    slope = (pylab.log(L1[-1])-pylab.log(L1[1]))/(pylab.log(Nx[-1])-pylab.log(Nx[1]))
    #Plot Results
    #---------------------------------------------------------------
    pylab.loglog(Nx,L1,'-o', label = "%s with slope= %f" % (methodName, slope))
    pylab.xlabel("Nx")
    pylab.ylabel("L1_norm")
    pylab.legend()
    
        
        








#Run advection for Square wave and Gaussian with both methods
#-------------------------------------------------------------------
pylab.figure()
initialConditions = ["Gaussian","Square"]
methodNames = ["Lax-Friedrichs", "Lax-Wendroff"]

Nx = 40

count = 1
for j in range(0,2):
    for i in range(0,2):
        pylab.subplot(2,2,count)
        advect(i+1,initialConditions[i],j+1,methodNames[j], Nx)
        count = count +1
        

    
#Run advection and plot convergence for both with both methods
#-------------------------------------------------------------------
Nxs = [10,30,100,300,1000,3000]
initialConditions = ["Gaussian","Square"]
methodNames = ["Lax-Friedrichs", "Lax-Wendroff"]

pylab.figure()

for j in range(0,2):
    pylab.subplot(2,1,j+1)
    pylab.title("Convergence rates for %s wave" % initialConditions[j])
    for i in range(0,2):
        advectConv(j+1,initialConditions[j],i+1,methodNames[i], Nxs)
        
pylab.show()
