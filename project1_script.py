# This python script reads from the command line the filename (its root) and
# the largest exponent of 10. This defines the number of mesh points.
# This script calls then an executable from a c++ or fortran code that solves
# a set of linear equations with a tridiagonal matrix defining the second derivative.
# It makes in turn the various plots as pdf files and finally sets up the basis for
# a report and its pertinent latex file

import sys, os
from  matplotlib import pyplot as plt
import numpy as np
import tabulate

print("Upon executing this script you will be prompted for a matrix size, n.")
print(" This size will determine the resolution of the stepsize used to solve our problem.")
print(" A secondary prompt will be launched to ask if you would solution plots or error plots")


# Command line arguments using sys.argv[]
try:
    filename = np.array(["Guassian", "LUdecomp"])
    matrix_size = []
except:
    print ("Usage of this script", sys.argv[0], "infile", sys.argv[1], "Exponent", sys.argv[2]); sys.exit(1)

# executes the c++ exe file, appending the matrix size i onto the command prompt
# loops through all members of matrix_size and executes main.exe for each i
'''
for i in matrix_size:
    cmdline = "main.exe " + str(i)
    cmd = cmdline
    failure = os.system(cmd)
    if failure:
        print('running main.exe failed'); sys.exit(1)
'''
# for a user interactive experience comment out the above and enable the below statement
# appends user input matrix sizes to the array matrix_size
print("Please enter matrix sizes seperated by a space")
gate = 1;
while gate == 1: # asks for a list of matrix sizes
    size = input("size list: ").split(' ');
    if type(size) == list:
        matrix_size=[int(num) for num in size];
        gate = 0;
    else:
        print("please enter an array")
        gate = 1;
for i in matrix_size: # Now run code, here c++ code  which has been compiled and linked
    cmdline  = "main.exe " + str(i)
    cmd = cmdline
    failure = os.system(cmd)
    if failure:
        print ('running project1 failed'); sys.exit(1)

# Start making figures looping over all exponents

# user prompt for figure creation
# solution generates a plot of the numerical and exact solutions for each matrix in matrix_size
# error produces a plot that looks at a log-log plot of the relative error vs matrix size

gate = 1;
while gate == 1:
        # declaring which figures are going to be made
    try:
        print("Let me know which plots you would like to make")
        figure_type = input("solution or error or flops: ")
        if figure_type == "solution":
            fig_type = 2;
        elif figure_type == "error":
            fig_type = 3;
        elif figure_type == "flops":
            fig_type = 4;
    except:
        print("you must enter, 'solution' or 'error', not whatever you put in just now")
    if fig_type==2:
        # iterates through the matrix sizes, creating figures for each method and matrix size
        i = j = 0;
        for i in np.arange(0,np.size(matrix_size)): # loops through matrix sizes
            for j in np.arange(0,np.size(filename)): # loops through methods
        #   define files to open data and make plots to
                if j == 0:
                    fout = filename[j]+str(int(matrix_size[i]+1)) # for the guassian that uses n + 1 points
                else:
                    fout = filename[j]+str(int(matrix_size[i])) # for the LU decomp that uses n points
                figfile = fout+"Solution.pdf"
                data = np.loadtxt(fout)
                x = data[:,0]
                solution = data[:,1]
                exact = data[:,2]
                plt.axis([0,1,0, 1.0])
                numericalplot = plt.plot(x, solution, 'r:.', linewidth = 2.0, label = 'Numerical')
                exactplot = plt.plot(x, exact, 'm:v', linewidth = 2.0, label = 'Exact')
                plt.xlabel(r'$x$')
                plt.ylabel(r'Numerical and Exact Solutions')
                plt.title('Solution Using Matrix size: ' + str(int(matrix_size[i])) + ' x ' + str(int(matrix_size[i])))
                plt.savefig(figfile)
        #       Then clean up
                plt.clf()
                j += 1;
        gate = 0; 
    if fig_type==3:
        error_list = np.zeros((np.size(matrix_size),np.size(filename)+1));
        for i in np.arange(0,np.size(matrix_size)):
            for j in np.arange(0,np.size(filename)):
                if j == 0:
                    fout = filename[j]+str(int(matrix_size[i]+1))
                else:
                    fout = filename[j]+str(int(matrix_size[i]))
                data = np.loadtxt(fout)
                avg_error = np.average(data[:,3]);
                error_list[i,j] = avg_error;
            error_list[i,2] = matrix_size[i]
        figfile = "errorlist.pdf"
        print(error_list)
        fig = plt.figure('error')
        errorplotGuass = plt.semilogx(matrix_size,error_list[:,0], 'r:.', linewidth=2.0, label ='Guass')
        errorplotLU = plt.semilogx(matrix_size,error_list[:,1], 'm:v', linewidth = 2.0, label = 'LUD')
        plt.xlabel("Logarithmic Matrix Size")
        plt.ylabel("Logarithmic Error")
        plt.title("log-log plot of Error and Matrix Size")
        plt.savefig(figfile)
        plt.clf()
        gate = 0;
    if fig_type == 4:
        flop_list = np.zeros((np.size(matrix_size),np.size(filename)+3));
        for i in np.arange(0,np.size(matrix_size)):
            for j in np.arange(0,np.size(filename)):
                if j == 0:
                    fout = filename[j]+str(int(matrix_size[i]+1))
                else:
                    fout = filename[j]+str(int(matrix_size[i]))
                data = np.loadtxt(fout)
                flops = np.average(data[:,4])
                time = np.average(data[:,5])
                flop_list[i,j] = flops
                flop_list[i,j+2] = time
            flop_list[i,4] = matrix_size[i]
        figfile = "flopplot.pdf"
        fig = plt.figure('flops')
        plt.subplot('211')
        flopsGuass = plt.plot(matrix_size,flop_list[:,0], 'r:.', linewidth=2.0, label ='Guass')
        plt.title("Total FLOPS to matrix size comparison")
        plt.ylabel('Guassian Elim')
        plt.subplot('212')
        flopsLU = plt.plot(matrix_size,flop_list[:,1]/10000, 'm:v', linewidth = 2.0, label = 'LUD')
        plt.xlabel("Matrix Size")
        plt.ylabel("LU Decomp [x10^3]")
        plt.savefig(figfile)
        plt.clf()
        polyflop = np.polyfit(matrix_size,flop_list[:,1],3)
        polynomial = ''
        for n in np.arange(0,4):
            temp = str(np.round(polyflop[n], 3))
            polynomial += str(temp)
            if n < 3:
                polynomial += 'x^'+str(3-n) + ' + '
            else:
                pass
        print('For the LU decomp method we find the polynomial of flops to size n as:')
        print(polynomial) 
        print(np.roll(flop_list.tolist().append(['Guassian FLOPS', 'LUDecomp FLOPS', 'Guassian Time', 'LUDecomp time', 'Matrix Size']),-1))
        #writing the flop list to a file
        temp = open('flop_list', 'w')
        header = "Matrix Size      Guassian: FLOPS     TIME      LU Decmop: FLOPS      TIME"
        temp.write(header + '\n')
        for jj in np.arange(0,np.size(matrix_size)):
            string = str(matrix_size[jj])
            for ii in [0,2,1,3]:
                string += ' ' + str(flop_list[jj,ii])
            temp.write(string)
            temp.write('\n')
        temp.close()
       
   # user prompt to enable the production of another set of figures
    # used if you want both solution/exact and error figures
    print("would you like another set of figures?")
    check = input("y/n: ")
    if check == "y":
        gate = 1;
    else:
        gate = 0;
    
# Now prepare latex file, r in front avoids backslashes being treated
# as control chars in strings. What follows are plain  latex commands
preamb = r"""\documentclass[10pt,showpacs,preprintnumbers,footinbib,amsmath,amssymb,aps,prl,twocolumn,groupedaddress,superscriptaddress,showkeys]{revtex4-1}
\usepackage{graphicx}
\usepackage{dcolumn}
\usepackage{bm}
\usepackage[colorlinks=true,urlcolor=blue,citecolor=blue]{hyperref}
\usepackage{color}
\begin{document}
\title{Project 1}
\author{A.~N.~Author}
\affiliation{Department of Something, University of Somewhere, Outer Space}
\begin{abstract}
We present our Ferrari algorithm for solving linear equations. Our best algorithm runs as $4n$ FLOPS with $n$ the dimensionality of the matrix.
\end{abstract}
\maketitle

"""

figure = r"""\begin{figure}[hbtp]
\includegraphics[scale=0.4]{test1.pdf}
\caption{Exact and numerial solutions for $n=10$ mesh points.} 
\label{fig:n10points}
\end{figure}

"""


introduction = r"""\section{Introduction}

"""

theory = r"""\section{Theory, algorithms and methods}

"""

results = r"""\section{Results and discussions}

"""

conclusions = r"""\section{Conclusions}

"""

references = r"""\begin{thebibliography}{99}
\bibitem{miller2006} G.~A.~Miller, A.~K.~Opper, and E.~J.~Stephenson, Annu.~Rev.~Nucl.~Sci.~{\bf 56}, 253 (2006).
\end{thebibliography}

"""

# Dump to file:
filename = 'ReportProject1'
f = open(filename + '.tex', "w")
f.write(preamb)
f.write(introduction)
f.write(theory)
f.write(results)
f.write(figure)
f.write(conclusions)
f.write(references)
f.write("""\end{document}""")
f.close()









