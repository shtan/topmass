#include "Shapes.h"

#include "TH1.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


//
// constructor and destructor
//

Shapes::Shapes( TH1D *&hmbl_bkg_temp ){

   //hists_ = temphists_;
   //hmbl_bkg = (TH1D*)hists_["mbl_fit"]["data_bkgcontrol"]->Clone("hmbl_bkg");
   hmbl_bkg = (TH1D*)hmbl_bkg_temp->Clone("hmbl_bkg");
   norm_mbl_bkg = hmbl_bkg->Integral("width");

   // GP options
   lmbl = 25.0;
   lmass = 5.0;
   int ntrain = 100;
   double rtrain = 1000;
   for(int i=0; i < ntrain; i++) ptrain.push_back( (i+0.5)*rtrain/ntrain );

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

   //return norm*(k*Fmbl_sig_param(x, mt)/integral + (1-k)*Fmbl_bkg(x));
   return norm*(k*Fmbl_sig_gp(x, mt)/integral + (1-k)*Fmbl_bkg(x));
}

double Shapes::Fmbl_sig_param(double x, double mt){

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

double Shapes::Fmbl_sig_gp(double x, double mt){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;

   double fgp = 0;
   for(unsigned int i=0; i < ptrain.size(); i++){
      for(int j=0; j < nmasses; j++){
         fgp += aGP[i+j*ptrain.size()]*GPkern( x, ptrain[i], lmbl, mt, masspnts[j], lmass );
      }
   }

   return fgp;
}

double Shapes::Fmbl_bkg(double x){

   return hmbl_bkg->Interpolate(x) / norm_mbl_bkg;
}

double Shapes::GPkern(double x1, double x2, double lx, double m1, double m2, double lm){

   double kernel = exp(-0.5*( pow( (x1-x2)/lx, 2) + pow( (m1-m2)/lm, 2) ));
   return kernel;
}

void Shapes::TrainGP( map< string, map<string, TH1D*> > & hists_ ){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;
   int ntrain = ptrain.size();

   // histograms
   vector<TH1D*> hgp;
   for(int i=0; i < nmasses; i++){

      stringstream ssmass;
      ssmass << floor(masspnts[i]);
      string smass = ssmass.str();

      hgp.push_back( (TH1D*)hists_["mbl_fit"]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp"+smass).c_str()) );
      hgp[i]->Scale( 1.0/hgp[i]->Integral("width") );

   }

   // compute covariance matrix
   TMatrixD K(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      for(int j=0; j < ntrain*nmasses; j++){
         int im = i % ntrain;
         int jm = j % ntrain;
         int imass = i / ntrain;
         int jmass = j / ntrain;
         K[i][j] = GPkern( ptrain[im], ptrain[jm], lmbl, masspnts[imass], masspnts[jmass], lmass );
     }
   }
   // compute noise matrix
   TMatrixD N(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      double binerr = hgp[imass]->GetBinError( hgp[imass]->FindBin(ptrain[im]) );
      for(int j=0; j < ntrain*nmasses; j++){
         if( i==j ){
            N[i][j] = pow( max(binerr,0.001), 2 );
         }else{
            N[i][j] = 0;
         }
     }
   }

   // inverse of sum
   TMatrixD A(ntrain*nmasses,ntrain*nmasses);
   A = K + N;
   A.Invert();

   // vector of training points
   TVectorD y(ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      y[i] = hgp[imass]->GetBinContent( hgp[imass]->FindBin(ptrain[im]) );
   }

   // alpha vector
   aGP.ResizeTo( ntrain*nmasses );
   aGP = A*y;
   
}
