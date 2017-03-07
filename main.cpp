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
//double funtion(double x) { return 100.0*exp(-10.0*x.at(i)); }
void print(Matrix& mat,int n){
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << mat(i, j) << " ";
		}
		cout << endl;
	}
}
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

	double h = 1.0 / (((double)n)+1);
	double hh = h*h;
	// initial matrix and vector values
	for (int i = 0; i < n; i++) {
		x.at(i) = (i+1)*h;
		analyticSol.at(i) = 1.0 - (1 - exp(-10))*x.at(i) - exp(-10 * x.at(i));
		f.at(i) = hh*100.0*exp(-10.0*x.at(i));
		
		for (int j = 0; j < n; j++) {
			mat_A(i, j) = 0.0;
				mat_L(i, j) = mat_U(i, j)= mat_LU(i,j) = 0.0;
		}
	}
	

	// center and right diagonals
	for (int i = 0; i < (n-1); i++) {
		mat_A(i,i) = 2.0; // center diagonal
		mat_A(i, i + 1) = -1.0; //right diagonal
	}

	for (int i = 1; i < n; i++) {
		mat_A(i, i - 1) = -1.0;
	}

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
                double gg=0.0;
                for(int q=0; q<i;q++){
                     gg+=mat_L(i,q)*mat_U(q,j);
                }
                mat_U(i,j)=mat_A(i,j)-gg;
            
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

	// Start timing
	clock_t start, finish;
	start = clock();

	// begin solution process
	vector<double> y(n); //defines an intermediate vector
	y.at(0) = f.at(0);
	cout << y.at(0) << " compared to: " << f.at(0) << endl;
	
	for (int i = 1; i < n; i++) { //solves for the intermediate vector values
        double temp1=0.0;
        for (int j = 0 ; j < i; j++) {
            temp1+=mat_L(i,j)*y.at(j);
        }
        
       y.at(i) = f.at(i) - temp1;
		
	}
	//new algorithm
    solution.at(n-1)=y.at(n-1)/mat_U(n-1,n-1);
    double temp2=0.0;
    for(int i=2;i<n+1;i++){
        temp2=0.0;
        for(int j=1;j<i;j++){
            temp2+=mat_U(n-i,n-j)*solution.at(n-j);
        }
        solution.at(n-i)=(y.at(n-i)-temp2)/mat_U(n-i,n-i);
    }
    
    
    
    
    /*
    for (int i = 0; i < n; i++) {
		y.at(i) = y.at(i) / mat_U(i, i);
	}
	for (int i = 2; i < n ; i++) {
		for (int j = 1; j < i; j++) {
			mat_U(n-i, n-j) = mat_U(n-i, n-j) / mat_U(n-i, n-i);
		}
	}
	// backward reduced row eschelon form
	for (int q = 2; q<n; q++){
		double temp = 0.0;
		// backward substitution will put A into reduced row eschelon form
		for (int j = 1; j < q; j++) {
			temp += mat_U(n - q, n - j)*y.at(n - q + 1) / mat_U(n - q + 1, n - q + 1);
			
		}
		y.at(n - q) = y.at(n - q) - temp;
	}
	for (int l = 1; l<n; l++){
		solution.at(l) = y.at(l);
	}
    
     */
	cout << "Solution       y value       U(n,n)       exact value" << endl;
	cout << solution.at(n - 1) <<" "<< y.at(n-1) <<" "<< mat_U(n-1,n-1) <<" "<< f.at(n - 1) << endl;
	/*
	cout << "matrix: U" << endl;
	print(mat_U, n);
	cout << "matrix: L" << endl;
	print(mat_L, n);
	cout << "matrix: LU" << endl;
	print(mat_LU, n);
	cout << "matrix: A" << endl;
	print(mat_A, n);
	*/
	// finish timing and print time
	finish = clock();
	double timeused = (double)(finish - start) / ((double)CLOCKS_PER_SEC);
	cout << setiosflags(ios::showpoint | ios::uppercase);
	cout << setprecision(10) << setw(20) << "Time used  for  computation with " << n << " elements =" << timeused << endl;
	
	//saving calculations and outputing to file
	cout << "begin writing to file - " <<  to_string(n) << " - elements" << endl;
	string filename = "LUdecomp";
	morten_wrfile(filename, solution, x, analyticSol, n);
	cout << "finished writing" << endl;
	
}

// writing function using morten's code
void morten_wrfile(string filename,vector<double>& u, vector<double>& x, vector<double>& exact,int n) {
	vector<double> solution = u; //brings in our solution from u
	string elements = to_string(n); //converts the number of elements to a string
	string fileout = filename + elements; //builds a filename as "filename_n=_elements.txt"

	outfile.open(fileout); // opens the file of interest
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	//outfile << "       x:             approx:          exact:       relative error" << endl; //table heading: labels columns
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
	
	gauss(11);
	gauss(101);
	gauss(1001);
	
	LUdecomp(10);
	LUdecomp(100);
	LUdecomp(1000);
	
	return 0;
}
