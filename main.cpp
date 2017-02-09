//
//  main.cpp
//  project uno
//
//  Created by Colin Gordon on 1/26/17.
//  Copyright © 2017 Colin Gordon. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "time.h"

// initialize namespace
using namespace std;
// setup an output file class
ofstream outfile;

void morten_wrfile(string filename, double *u, double *x, double* exact, int n); //write function defined below

void writeoutput (double * u, double * x, int n); //output function defined below

void urkle(double * u, double * q,int n); //error calculation function defined below

void gauss(int n); //backward-forward substitution, ie: gaussian elimination

//function definitions
//double exact(double x) { return 1.0 - (1 - exp(-10))*x - exp(-10 * x); }

int main(int argc, const char * argv[]){
    
    gauss(11);
	cout << "begin gause 101" << endl;
    gauss(101);
	cout << "begin gause 1001" << endl;
	gauss(1001);
    /*
    for(int i=0;i<101;i++){
    int n = (int) (10000.0*pow(10,(double) i/100.0));
        gauss(n);
    }
	*/
    return 0;
}

void gauss(int n) {

   double h = 1.0/(double)n; //step size
   double* f = new double[n]; //function
   double* b = new double[n]; //diagonal
   double* a = new double[n-1]; //upper
   double* c = new double[n-1]; //lower
   double* u = new double[n+1]; // numerical solution
   double* exact = new double[n];  // exact solution
   double* x = new double[n+1]; //x axis
    u[0]=0.0;
    u[n]=0.0;
    x[0]=0.0;
    x[n]=1.0;
	
	// establish initial data
	for (int i = 0; i <= n; i++) {
		x[i] = i*h; //creates our xdata
		f[i] = h*h*100.0*exp(-10.0*x[i]); //creates our funtion data
		exact[i] = 1.0 - (1 - exp(-10))*x[i] - exp(-10 * x[i]); //creates our exact solutions data
		b[i] = 2.0; //the diagonal terms
	}
	// Start timing
	clock_t start, finish;
	start = clock();

	//begin for-loop to do guassian solution
    
    for (int j=1; j< n-1; j++) //tridiagonal terms
    {
        a[j]= -1.0; //upper
        
        c[j] = -1.0; //lower
    }
    
    for (int i=2; i<n; i++) //foward substitution
    {
        // these calculations will set the lower terms to zero while altering the rest of the matrix accordingly
        b[i]=b[i]-a[i-1]*c[i-1]/b[i-1];
        
        f[i]=f[i]-f[i-1]*c[i-1]/b[i-1];
    }
    for (int q=2; q<n; q++) //backward substitution
    {
		// backward substitution will put A into reduced row eschelon form
        f[n-q]=f[n-q]-a[n-q]*f[n-q+1]/b[n-q+1];
    }
    
    
    for (int l=1; l<n; l++)
    {
        u[l]= f[l]/b[l];
        
        
        //std::cout << "result at step "<< l << " " << u[l] << "  "<< "error "<< (q[l]-u[l])/q[l] << std::endl;
    }
	// finish timing and print time
	finish = clock();
	double timeused = (double)(finish - start) / ((double)CLOCKS_PER_SEC);
	cout << setiosflags(ios::showpoint | ios::uppercase);
	cout << setprecision(10) << setw(20) << "Time used  for  computation with " << n << " elements =" << timeused << endl;

	//begin data output
	/*
    writeoutput(u,x,n); //output function made by colin
	cout << "writeoutput" << endl;
    urkle(u,exact,n); //error output function made by colin
	cout << "urkle" << endl;
	*/
	//data output by spencer
	string filename = "Guassian"; //declares a filename fitting the solution method
	string elements = to_string(n); //converts the number of elements, n, to a string
	morten_wrfile(filename, u, x, exact, n);
	cout << "completed morten writefile " << n << " elements" << endl;
	//empty the arrays to free memory
	delete[] f;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] exact;
	delete[] u;
	delete[] x;
	cout << "complted deleting stuff" << endl;

}

// writing function using morten's code
void morten_wrfile(string filename, double *u, double *x, double* exact,int n) {
	double* solution = u; //brings in our solution from u
	string elements = to_string(n); //converts the number of elements to a string
	string fileout = filename + "_n=_" + elements + ".txt"; //builds a filename as "filename_n=_elements.txt"

	outfile.open(fileout); // opens the file of interest
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	outfile << "       x:             approx:          exact:       relative error" << endl;
	for (int i = 1; i < n; i++) {
		double RelativeError = fabs((exact[i] - solution[i]) / exact[i]);
		outfile << setw(15) << setprecision(8) << x[i];
		outfile << setw(15) << setprecision(8) << solution[i];
		outfile << setw(15) << setprecision(8) << exact[i];
		outfile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
	}
	outfile.close();
}

void writeoutput (double * u, double * x, int n){
    
    std::string filename;
   
    if (n==11){
        filename += "file10.csv";
    }
    else if (n==101){
            filename += "file100.csv";
        }
    else if (n==1001){
            filename += "file1000.csv";
        }
    else{
        return;
    }
    
    std::ofstream myfile;
    myfile.open(filename);
    for (int i=0; i<n+1; i++){
        myfile << x[i] << "," << u[i] << std::endl;
    }
    myfile.close();
    }
    
    
    
    
void urkle(double * u, double * q,int n){ //error function
    
    double e = std::abs((q[1]-u[1])/q[1]);
    
    for (int i=2; i<n; i++){
        double g=std::abs((q[i]-u[i])/q[i]);
        if(g>e){
            e=g;
        }
    }
    //std::cout << "result at step "<< n << " " << " max error is " u[l] << std::endl
    std::string filename = "error.csv";
    
    std::ofstream myfile;
    myfile.open(filename,std::ofstream::out | std::ofstream::app);
    myfile << -log10((double) n) << "," << log10(e) << std::endl;
    
    
    myfile.close();
}


