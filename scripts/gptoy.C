#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrix.h"
#include "TDecompChol.h"
#include "TRandom3.h"

#include <iostream>
using namespace std;

double kern(double x, double y, int type){

   if( type==0 ) return 1.0*x*y;
   else if( type==1 ) return exp(-0.5*pow( (x-y)/0.15, 2));
   else return -1;
}

void gptoy(){

   int n = 100;
   TRandom3 rand(0);

   // sampling points
   double x [n];
   for(int i=0; i < n; i++) x[i] = 1.0*i/n;

   // construct covariance matrix
   TMatrixD C(n,n);
   for(int i=0; i < n; i++)
      for(int j=0; j < n; j++)
         C[i][j] = kern(x[i],x[j],1);

   // add small number to diagonal (positive definite)
   for(int i=0; i < n; i++)
      C[i][i] += 0.00001;

   // use Cholesky decomposition to obtain vector distributed as multivariate gaussian
   TDecompChol D(C);
   TMatrixD L(n,n);
   L = D.GetU();

   // random rumbers
   TVectorD u(n);
   for(int i=0; i < n; i++) u[i] = rand.Gaus(0,1);

   // result -- random vector
   TVectorD y(n);
   y = L*u;

   // plot
   /*
   TGraph *grand = new TGraph();
   for(int i=0; i < n; i++){
      grand->SetPoint(i, x[i], y[i]);
   }

   grand->SetLineWidth(2);
   grand->Draw("ACP");
   */

   //
   // prediction w/ noise-free observables
   //

   int m = 5;
   double t [m];
   t[0] = 0.05;
   t[1] = 0.25;
   t[2] = 0.35;
   t[3] = 0.55;
   t[4] = 0.75;

   TVectorD ft(m);
   ft[0] = 1.1;
   ft[1] = 0.1;
   ft[2] = 0.5;
   ft[3] = -1.0;
   ft[4] = -0.5;

   TMatrixD K11(m,m);
   TMatrixD K12(m,n);
   TMatrixD K21(n,m);
   TMatrixD K22(n,n);

   for(int i=0; i < n; i++){
      for(int j=0; j < n; j++){

         if( i < m and j < m ){
            K11[i][j] = kern(t[i],t[j],1);
            if( i==j ) K11[i][j] += 0.00001;
         }
         if( i < m and j < n ){
            K12[i][j] = kern(t[i],x[j],1);
            if( i==j ) K12[i][j] += 0.00001;
         }
         if( i < n and j < m ){
            K21[i][j] = kern(x[i],t[j],1);
            if( i==j ) K21[i][j] += 0.00001;
         }
         if( i < n and j < n ){
            K22[i][j] = kern(x[i],x[j],1);
            if( i==j ) K22[i][j] += 0.00001;
         }

      }
   }

   // use Cholesky decomposition to obtain vector distributed as multivariate gaussian
   
   // compute covariance matrix
   TMatrixD fC(n,n);
   TMatrixD K11inv(m,m);
   K11inv = K11;
   K11inv.Invert();
   fC = K22 - K21*K11inv*K12;
   TDecompChol fD(fC);
   TMatrixD fL(n,n);
   fL = fD.GetU();

   //
   // result
   //
   
   TVectorD v(n);
   TVectorD fy(n);
   TVectorD fy_err(n);
   // mean value and error band
   for(int i=0; i < n; i++) v[i] = 1.0;
   fy = K21*K11inv*ft;
   fy_err = fL*v;

   TGraphErrors *gtest = new TGraphErrors();
   for(int i=0; i < n; i++){
      gtest->SetPoint(i, x[i], fy[i]);
      gtest->SetPointError(i, 0.0, fy_err[i]);
   }

   TCanvas * canvas = new TCanvas("canvas","canvas",800,800);
   canvas->cd();

   TGraph *gtrain = new TGraph();
   for(int i=0; i < m; i++) gtrain->SetPoint(i, t[i], ft[i]);

   gtest->SetFillStyle(3002);
   gtest->Draw("AC3");
   gtrain->SetMarkerStyle(20);
   gtrain->SetMarkerColor(2);
   gtrain->Draw("P");

   canvas->Draw();

   return;
}


