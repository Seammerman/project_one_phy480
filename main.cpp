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


using namespace std;

void writeoutput (double * u, double * x, int n); //output function defined below

void urkle(double * u, double * q,int n); //error calculation function defined below

void gauss(int n); //backward-forward substitution, ie: gaussian elimination

int main(int argc, const char * argv[]){
    
    gauss(11);
    gauss(101);
    gauss(1001);
    
    
    
    for(int i=0;i<101;i++){
    int n = (int) (10000.0*pow(10,(double) i/100.0));
        gauss(n);
        
        
        
    }

    return 0;
}

void gauss(int n) {

   double h = 1.0/(double)n; //step size
   double* f = new double[n]; //function
   double* b = new double[n]; //diagonal
   double* a = new double[n-1]; //upper
   double* c = new double[n-1]; //lower
   double* u = new double[n+1]; // numerical solution
   double* q = new double[n];  // exact solution
   double* x = new double[n+1]; //x axis
    const double e1=1.0-exp(-10.0);
    u[0]=0.0;
    u[n]=0.0;
    x[0]=0.0;
    x[n]=1.0;
    
    
    for (int i=1; i<n;i++)
    
    {
        f[i]= h*h*100.0*exp(-10.0*(double)i*h); //the inhomogenous function times step size squared.
        
        q[i]= -exp(-10.0*(double)i*h)+1.0-e1*(double)i*h; //exact solution
        
        x[i]= (double)i*h; //x axis
        
        b[i]= 2.0; //the diagonal terms
    }
    
    
    
    
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
    writeoutput(u,x,n);
    urkle(u,q,n);
    
	//empty the arrays to free memory
	delete[] f;
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] q;
	delete[] u;
	delete[] x;

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


