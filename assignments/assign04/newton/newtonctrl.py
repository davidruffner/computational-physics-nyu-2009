import pylab, optparse, os

#Method that runs the newton program in C for a certain method, initial con
#-ditions and NP. It makes a plot of the analytical and estimated solutions
# and it also results
#----------------------------------------------------------------------------

def newton():
    #Initialize parmeter file
    #--------------------------------------------------------
    filename = "param.cfg"
    pfile = open(filename, 'w')
    pfile.write("NP = 50 \n")
    pfile.write("P0 = 1.0 \n")
    pfile.write("P1 = 10.0 \n")
    pfile.write("a = 1.0 \n")
    pfile.write("b = 0.1 \n")
    pfile.write("T0 = 1.0 \n ")
    pfile.write("tol = 0.000001 ")
    pfile.close()

    #  run c program
    #-----------------------------
    os.system("./newton_dbr250")

    #Plot results
    #-------------------------------------------------------
    pylab.title("T vs P, for given values of P")

    infile_name = "TvsP.dat"
    infile = open(infile_name)
    lines = infile.readlines()

    strvals = [l.strip().split() for l in lines]
    x       = [float(v[0]) for v in strvals]
    y       = [float(v[1]) for v in strvals]

    
    pylab.plot(x,y,'-o')

    pylab.xlabel("P")
    pylab.ylabel("T")
    
    #----------------------------------------------------------

#Run progam
newton()
pylab.show()
