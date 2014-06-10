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

Shapes::Shapes( double gplength_mbl, double gplength_mt,
     double lbound_mbl, double rbound_mbl, double norm ){

   // GP options
   lmbl = gplength_mbl;
   lmass = gplength_mt;
   lbmbl = lbound_mbl;
   rbmbl = rbound_mbl;
   // TODO
   gnorm = norm;
   int ntrain = 100;
   double rtrain = 300;
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
   double integralsig = pp[3];
   double integralbkg = pp[4];

   double val = norm*(k*Fmbl_gp(x, mt, "sig")/integralsig + (1-k)*Fmbl_gp(x, mt, "bkg")/integralbkg);
   if( val <= 0 or (x > lbmbl and x < rbmbl) ) return 1E-10;
   else return val;
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

double Shapes::Fmbl_gp(double x, double mt, string sb){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;

   double fgp = 0;
   for(unsigned int i=0; i < ptrain.size(); i++){
      for(int j=0; j < nmasses; j++){
         double agp = 0;
         if( sb.compare("sig") == 0 ) agp = aGPsig[i+j*ptrain.size()];
         else if( sb.compare("bkg") == 0 ) agp = aGPbkg[i+j*ptrain.size()];
         else{
            cout << "ERROR in GP shape." << endl;
            return -1;
         }
         fgp += agp*GPkern( x, ptrain[i], lmbl, mt, masspnts[j], lmass );
      }
   }

   return fgp;
}

// TODO
double Shapes::Fmbl_gp_var(double x, double mt, string sb){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;
   int ntrain = ptrain.size();

   // vector of covariances
   TVectorD k(ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      k[i] = GPkern( x, ptrain[im], lmbl, mt, masspnts[imass], lmass );
   }
   TVectorD kT = k;

   double c1=0, c2=0;
   c1 = GPkern( x, x, lmbl, mt, mt, lmass );
   if( sb.compare("sig") == 0 ) k *= Asig;
   else if( sb.compare("bkg") == 0 ) k *= Abkg;
   else{
      cout << "ERROR in GP shape." << endl;
      return -1;
   }
   c2 = kT*k;

   return (c1-c2);
}

double Shapes::GPkern(double x1, double x2, double lx, double m1, double m2, double lm ){

   double kernel = exp(-0.5*(pow( (x1-x2)/lx, 2)+pow( (m1-m2)/lm, 2)));
   return kernel;
}

void Shapes::TrainGP( map< string, map<string, TH1D*> > & hists_ ){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;
   int ntrain = ptrain.size();

   // histograms
   vector<TH1D*> hgp_sig;
   vector<TH1D*> hgp_bkg;
   for(int i=0; i < nmasses; i++){

      stringstream ssmass;
      ssmass << floor(masspnts[i]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)hists_["mbl"]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sig"+smass).c_str()) );
      hgp_sig[i]->Scale( 1.0/hgp_sig[i]->Integral("width") );

      // background shape
      hgp_bkg.push_back( (TH1D*)hists_["mbl"]["ttbar"+smass+"_mistag"]
            ->Clone( ("hgp_bkg"+smass).c_str()) );
      hgp_bkg[i]->Add( hists_["mbl"]["ttbar"+smass+"_hadronic"] );
      hgp_bkg[i]->Add( hists_["mbl"]["ttbar"+smass+"_taus"] );
      hgp_bkg[i]->Add( hists_["mbl"]["other"] );
      hgp_bkg[i]->Scale( 1.0/hgp_bkg[i]->Integral("width") );

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
   TMatrixD Nsig(ntrain*nmasses,ntrain*nmasses);
   TMatrixD Nbkg(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      double binerr_sig = hgp_sig[imass]->GetBinError( hgp_sig[imass]->FindBin(ptrain[im]) );
      double binerr_bkg = hgp_bkg[imass]->GetBinError( hgp_bkg[imass]->FindBin(ptrain[im]) );
      // TODO
      binerr_sig *= sqrt(gnorm);
      binerr_bkg *= sqrt(gnorm);
      cout << binerr_sig << " " << binerr_bkg << endl;
      for(int j=0; j < ntrain*nmasses; j++){
         if( i==j ){
            Nsig[i][j] = pow( max(binerr_sig,0.001), 2 );
            Nbkg[i][j] = pow( max(binerr_bkg,0.001), 2 );
         }else{
            Nsig[i][j] = 0;
            Nbkg[i][j] = 0;
         }
     }
   }

   // inverse of sum
   Asig.ResizeTo( ntrain*nmasses, ntrain*nmasses );
   Abkg.ResizeTo( ntrain*nmasses, ntrain*nmasses );
   Asig = K + Nsig;
   Abkg = K + Nbkg;
   Asig.Invert();
   Abkg.Invert();

   // TODO
   TMatrixD Asigtemp = Asig;
   TMatrixD Abkgtemp = Abkg;

   // vector of training points
   TVectorD ysig(ntrain*nmasses);
   TVectorD ybkg(ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      ysig[i] = hgp_sig[imass]->GetBinContent( hgp_sig[imass]->FindBin(ptrain[im]) );
      ybkg[i] = hgp_bkg[imass]->GetBinContent( hgp_bkg[imass]->FindBin(ptrain[im]) );
   }

   // alpha vector
   aGPsig.ResizeTo( ntrain*nmasses );
   aGPsig = Asigtemp*ysig;

   aGPbkg.ResizeTo( ntrain*nmasses );
   aGPbkg = Abkgtemp*ybkg;

}
