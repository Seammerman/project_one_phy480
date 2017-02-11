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
#include <vector>
#include "MatrixHead.h"

// initialize namespace
using namespace std;
// setup an output file class
ofstream outfile;

void morten_wrfile(string filename, vector<double>& u, vector<double>& x, vector<double>& exact, int n); //write function defined below
void LUdecomp(int n);
void gauss(int n); //backward-forward substitution, ie: gaussian elimination
void LUdecomp(int n); //solution by LU composition

//function definitions
//double exact(double x) { return 1.0 - (1 - exp(-10))*x - exp(-10 * x); }
//funtion 100.0*exp(-10.0*x.at(i))

void gauss(int n) {

   double h = 1.0/(double)n; //step size
   /*
   double *f = new double[n]; //function
   double *b = new double[n]; //diagonal
   double *a = new double[n-1]; //upper
   double *c = new double[n-1]; //lower
   double *u = new double[n+1]; // numerical solution
   double *exact = new double[n];  // exact solution
   double *x = new double[n+1]; //x axis
    u[0]=0.0;
    u[n]=0.0;
    x[0]=0.0;
    x[n]=1.0;
	*/

   vector<double> f(n);
   vector<double> b(n);
   vector<double> a(n-1);
   vector<double> c(n-1);
   vector<double> u(n+1);
   vector<double> exact(n);
   vector<double> x(n+1);
  
   u.at(0) = 0.0;
   u.at(n) = 0.0;
   x.at(0) = 0.0;
   x.at(n) = 1.0;
   
	// establish initial data
	for (int i = 0; i < n; i++) {
		x.at(i) = i*h; //creates our xdata
		f.at(i) = h*h*100.0*exp(-10.0*x.at(i)); //creates our funtion data
		exact.at(i) = 1.0 - (1 - exp(-10))*x.at(i) - exp(-10 * x.at(i)); //creates our exact solutions data
		b.at(i) = 2.0; //the diagonal terms
	}
	// Start timing
	clock_t start, finish;
	start = clock();

	//begin for-loop to do guassian solution
    
    for (int j=1; j< n-1; j++) //tridiagonal terms
    {
        a.at(j)= -1.0; //upper
        
        c.at(j) = -1.0; //lower
    }
    for (int i=2; i<n; i++) //foward substitution
    {
        // these calculations will set the lower terms to zero while altering the rest of the matrix accordingly
        b.at(i)=b.at(i)-a.at(i-1)*c.at(i-1)/b.at(i-1);
        
        f.at(i)=f.at(i)-f.at(i-1)*c.at(i-1)/b.at(i-1);
    }
    for (int q=2; q<n; q++) //backward substitution
    {
		// backward substitution will put A into reduced row eschelon form
        f.at(n-q)=f.at(n-q)-a.at(n-q)*f.at(n-q+1)/b.at(n-q+1);
    }
    for (int l=1; l<n; l++)
    {
        u.at(l)= f.at(l)/b.at(l);
       
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
	string filename = "Guassian";										//declares a filename fitting the solution method
	string elements = to_string(n);										//converts the number of elements, n, to a string
	morten_wrfile(filename, u, x, exact, n);
	cout << "completed morten writefile " << n << " elements" << endl;
}



//print matricies

void printmatrix(Matrix& amatrix, int n){
    
    for(int j=0;j<n;j++){
        for(int k=0;k<n;k++){
            cout<<amatrix(j,k) << "  ";
        
        }
        cout<<endl;
    }
    cout << endl;
    
    
}










//definition of the LUdecomp function

void LUdecomp(int n) {
	
	Matrix mat_A(n, n);
	Matrix mat_L(n, n);
	Matrix mat_U(n, n);
    Matrix mat_LU(n,n);
	vector<double> x(n);
	vector<double> solution(n);
	vector<double> f(n);
	vector<double> analyticSol(n);

	double h = 1.0 / (double)n;
	double hh = h*h;
	// filling the matrices with zeros
	for (int i = 0; i < n; i++) {
		x.at(i) = i*h;
		analyticSol.at(i) = 1.0 - (1 - exp(-10))*x.at(i) - exp(-10 * x.at(i));
		f.at(i) = hh*100.0*exp(-10.0*x.at(i));

		for (int j = 0; j < n; j++) {
			mat_A(i, j) = mat_L(i, j) = mat_U(i, j)= mat_LU(i,j) = 0.0;
		}
	}
	//filling the tridiagonal elements
	for (int i = 0; i < n-1; i++) {
		mat_A(i, i) = 2.0; // center diagonal
		mat_A(i, i + 1) = -1.0; //right diagonal
	}
	//left diagonal
	for (int i = 1; i < n; i++) {
		mat_A(i, i - 1) = -1.0;
	}
	//fill in the last diagonal
	mat_A(n-1, n-1) = 2.0;
	
    for (int i=0; i<n; i++) {
        mat_U(0,i)=mat_A(0,i);
        mat_L(i,0)=mat_A(i,0)/mat_U(0,0);
    }
    for (int i=1;i<n;i++){
        mat_L(i,i)=1.0;
    }
    
    for(int j=1;j<n;j++){
        for(int i=1;i<n;i++){
         
            if(i>j) {
                double s=0.0;
                for(int p=0;p<j;p++){
                    s+=mat_L(i,p)*mat_U(p,j);
                }
                mat_L(i,j)=(mat_A(i,j)-s)/mat_U(j,j);
            }
            
            else if (i==j){
                double g=0.0;
                for(int p=0;p<i;p++){
                     g+=mat_L(i,p)*mat_U(p,j);
                }
                mat_U(i,j)=mat_A(i,j)-g;
            
            }
            else{
                double h=0.0;
                for(int q=0; q<i;q++){
                     h+=mat_L(i,q)*mat_U(q,j);
                }
                mat_U(i,j)=mat_A(i,j)-h;
            
            }
            
        }
    }
    
    
    for (int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            double s=0.0;
            for(int k=0;k<n;k++){
                s+=mat_L(i,k)*mat_U(k,j);
            }
            mat_LU(i,j)=s;
         }
    }
    
    printmatrix(mat_U,n);
    printmatrix(mat_L,n);
    printmatrix(mat_A,n);
    printmatrix(mat_LU,n);
    
    //cout << "Leaving LUdecomp"<< endl;
    
}

// writing function using morten's code
void morten_wrfile(string filename,vector<double>& u, vector<double>& x, vector<double>& exact,int n) {
	vector<double> solution = u; //brings in our solution from u
	string elements = to_string(n); //converts the number of elements to a string
	string fileout = filename + "_n=_" + elements + ".txt"; //builds a filename as "filename_n=_elements.txt"

	outfile.open(fileout); // opens the file of interest
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	outfile << "       x:             approx:          exact:       relative error" << endl; //table heading: labels columns
	for (int i = 1; i < n; i++) {
		double RelativeError = fabs((exact.at(i) - solution.at(i)) / exact.at(i));
		outfile << setw(15) << setprecision(8) << x.at(i);
		outfile << setw(15) << setprecision(8) << solution.at(i);
		outfile << setw(15) << setprecision(8) << exact.at(i);
		outfile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
	}
	outfile.close();
}

int main(int argc, const char * argv[]) {
	//gauss(11);
	//gauss(101);
	//gauss(1001);

	LUdecomp(10);
	return 0;
}
