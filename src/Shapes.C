#include "Shapes.h"

#include "TH1.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;


//
// constructor and destructor
//

Shapes::Shapes(TH1D *&hist){
   hmbl_bkg = hist;
}

Shapes::~Shapes(){
}


//
// member definitions
//

double Shapes::Fmbl_tot(double *px, double *pp){

   double x = px[0];
   double mt = pp[0];
   double k = pp[1];
   double norm = pp[2];
   double integral = pp[3];

   return norm*(k*Fmbl_sig(x, mt)/integral + (1-k)*Fmbl_bkg(x));
}

double Shapes::Fmbl_sig(double x, double mt){

   double par [9];
   par[0] = 0.04314 - 0.0001872*mt;
   par[1] = -89.26 + 1.139*mt;
   par[2] = 20.86 + 0.004963*mt;
   par[3] = -0.005268 + 7.442E-05*mt;
   par[4] = -90.66 + 0.9241*mt;
   par[5] = -26.65 + 0.2591*mt;
   par[6] = 0.02192 + -6.172E-05*mt;
   par[7] = -23.11 + 0.3386*mt;
   par[8] = -15.22 + 0.1252*mt;

   double gaus1 = par[0]*exp( -0.5*pow((x-par[1])/par[2],2) );
   double gaus2 = par[3]*exp( -0.5*pow((x-par[4])/par[5],2) );
   double landau = par[6]*TMath::Landau(x,par[7],par[8],false);

   return gaus1+gaus2+landau;
}

double Shapes::Fmbl_bkg(double x){

   return hmbl_bkg->Interpolate(x) / hmbl_bkg->Integral("width");
}
