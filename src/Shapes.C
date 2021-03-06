#include "Shapes.h"

#include "TH1.h"
#include "TMath.h"
#include "TDecompLU.h"
#include "TDecompChol.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;


//
// constructor and destructor
//

Shapes::Shapes( string var, double gplength_x, double gplength_mt, double norm1, double norm2 ){

   name = var;
   // GP options
   lx = gplength_x;
   lmass = gplength_mt;
   gnorm1 = norm1;
   gnorm2 = norm2;
   int ntrain = 100;
   double rtrain = 300;
   for(int i=0; i < ntrain; i++) ptrain.push_back( (i+0.5)*rtrain/ntrain );

   // right and left bounds -- set to zero unless needed
   lbx = 0.0;
   rbx = 0.0;

   // flag for cross validation two-stage fit
   do_gpvar = false;

}

Shapes::~Shapes(){
}


//
// member definitions
//

double Shapes::Ftot(double *px, double *pp){

   double x = px[0];
   double mt = pp[0];
   double k = pp[1];
   double norm = pp[2];
   double integralsig = pp[3];
   double integralbkg = pp[4];

   double val = norm*(k*Fmbl_gp(x, mt, "sig")/integralsig + (1-k)*Fmbl_gp(x, mt, "bkg")/integralbkg);
   if( val <= 0 or (x > lbx and x < rbx) ) return 1E-10;
   else return val;
}

double Shapes::Fsig_param(double x, double mt){

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
         fgp += agp*GPkern( x, ptrain[i], lx, mt, masspnts[j], lmass );
      }
   }

   return fgp;
}

double Shapes::Fmbl_gp_var(double x, double mt, string sb){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;
   int ntrain = ptrain.size();

   // vector of covariances
   TVectorD k(ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      k[i] = GPkern( x, ptrain[im], lx, mt, masspnts[imass], lmass );
   }
   TVectorD kT = k;

   double c1=0, c2=0;
   c1 = GPkern( x, x, lx, mt, mt, lmass );
   if( sb.compare("sig") == 0 ) k *= Ainv_sig;
   else if( sb.compare("bkg") == 0 ) k *= Ainv_bkg;
   else{
      cout << "ERROR in GP shape." << endl;
      return -1;
   }
   c2 = kT*k;

   return (c1-c2);
}

double Shapes::GPkern(double x1, double x2, double lsx, double m1, double m2, double lsm ){

   double kernel = 1E-06*gnorm2*gnorm1*exp(-0.5*(pow( (x1-x2)/lsx, 2)+pow( (m1-m2)/lsm, 2)));
   return kernel;
}

void Shapes::TrainGP( map< string, map<string, TH1D*> > & hists_,
     double &m2llsig, double &m2llbkg ){

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
      hgp_sig.push_back( (TH1D*)hists_[name]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sig"+smass).c_str()) );
      hgp_sig[i]->Scale( 1.0/hgp_sig[i]->Integral("width") );

      // background shape
      hgp_bkg.push_back( (TH1D*)hists_[name]["ttbar"+smass+"_mistag"]
            ->Clone( ("hgp_bkg"+smass).c_str()) );
      hgp_bkg[i]->Add( hists_[name]["ttbar"+smass+"_hadronic"] );
      hgp_bkg[i]->Add( hists_[name]["ttbar"+smass+"_taus"] );
      hgp_bkg[i]->Add( hists_[name]["other"] );
      hgp_bkg[i]->Scale( 1.0/hgp_bkg[i]->Integral("width") );

      for(int n=0; n < hgp_sig[i]->GetNbinsX(); n++){
         if( hgp_sig[i]->GetBinError(n) < 5E-06 ) hgp_sig[i]->SetBinError(n, 5E-06);
         if( hgp_bkg[i]->GetBinError(n) < 5E-06 ) hgp_bkg[i]->SetBinError(n, 5E-06);
      }
      

   }

   // compute covariance matrix
   TMatrixD K(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      for(int j=0; j < ntrain*nmasses; j++){
         int im = i % ntrain;
         int jm = j % ntrain;
         int imass = i / ntrain;
         int jmass = j / ntrain;
         K[i][j] = GPkern( ptrain[im], ptrain[jm], lx, masspnts[imass], masspnts[jmass], lmass );
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
      binerr_sig *= sqrt(gnorm2);
      binerr_bkg *= sqrt(gnorm2);
      for(int j=0; j < ntrain*nmasses; j++){
         if( i==j ){
            Nsig[i][j] = binerr_sig*binerr_sig;//pow( max(binerr_sig,0.001), 2 );
            Nbkg[i][j] = binerr_bkg*binerr_bkg;//pow( max(binerr_bkg,0.001), 2 );
         }else{
            Nsig[i][j] = 0;
            Nbkg[i][j] = 0;
         }
     }
   }

   // inverse of sum
   TMatrixD Asig = K + Nsig;
   TMatrixD Abkg = K + Nbkg;
   TDecompChol Cholsig(Asig);
   TDecompChol Cholbkg(Abkg);
   TMatrixD AsigU = Cholsig.GetU();
   TMatrixD AbkgU = Cholbkg.GetU();
   bool status = 0;
   TMatrixDSym Asinv_sig = Cholsig.Invert(status);
   TMatrixDSym Asinv_bkg = Cholbkg.Invert(status);
   Ainv_sig.Clear();
   Ainv_bkg.Clear();
   Ainv_sig.ResizeTo( ntrain*nmasses, ntrain*nmasses );
   Ainv_bkg.ResizeTo( ntrain*nmasses, ntrain*nmasses );
   Ainv_sig = (TMatrixD)Asinv_sig;
   Ainv_bkg = (TMatrixD)Asinv_bkg;


   TMatrixD Ainv_sigtemp = Ainv_sig;
   TMatrixD Ainv_bkgtemp = Ainv_bkg;

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
   aGPsig.Clear();
   aGPbkg.Clear();

   aGPsig.ResizeTo( ntrain*nmasses );
   aGPsig = Ainv_sigtemp*ysig;

   aGPbkg.ResizeTo( ntrain*nmasses );
   aGPbkg = Ainv_bkgtemp*ybkg;

   // compute marginal likelihood
   //TDecompLU lusig(Asig);
   //TDecompLU lubkg(Abkg);
   //TMatrixD AsigLU = lusig.GetLU();
   //TMatrixD AbkgLU = lubkg.GetLU();

   double ldetsig = 0.0;
   double ldetbkg = 0.0;
   for(int i=0; i < Cholsig.GetNrows(); i++){
      ldetsig += 2*log(AsigU[i][i]);
      ldetbkg += 2*log(AbkgU[i][i]);
   }

   double term1sig = -0.5*ysig*aGPsig;
   double term2sig = -0.5*ldetsig;
   double term3sig = -0.5*ntrain*log(2*TMath::Pi());

   double term1bkg = -0.5*ybkg*aGPbkg;
   double term2bkg = -0.5*ldetbkg;
   double term3bkg = -0.5*ntrain*log(2*TMath::Pi());

   m2llsig = -2.0*(term1sig+term2sig+term3sig);
   m2llbkg = -2.0*(term1bkg+term2bkg+term3bkg);

   return;
}

void Shapes::LearnGPparams( map< string, map<string, TH1D*> > & hists_ ){

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetPrintLevel(3);

   // set training hist
   hists_train_ = &hists_;

   fFunc = new ROOT::Math::Functor ( this, &Shapes::GPm2llX, 4 );
   gMinuit->SetTolerance(1000.0);
   gMinuit->SetFunction( *fFunc );

   // stage 1 -- use bin errros as variance
   cout << "**************** Cross Validation Round 1 *****************" << endl;
   do_gpvar = false;
   gMinuit->SetLowerLimitedVariable(0, "gpnorm1", 5.0, 0.1, 0.0);
   gMinuit->SetLowerLimitedVariable(1, "gpnorm2", 10.0, 0.1, 0.0);
   gMinuit->SetLowerLimitedVariable(2, "lx", 15, 1, 0.0);
   gMinuit->SetLowerLimitedVariable(3, "lmass", 30, 1, 0.0);

   gMinuit->Minimize();

   cout << "**************** Cross Validation Round 2 *****************" << endl;
   do_gpvar = false;
   gMinuit->SetLimitedVariable(0, "gpnorm1", xstmp[0], 1, 0.5*xstmp[0], 2*xstmp[0]);
   gMinuit->SetLowerLimitedVariable(1, "gpnorm2", xstmp[1], 1, 0.0);
   gMinuit->SetLimitedVariable(2, "lx", xstmp[2], 1, 0.8*xstmp[2], 1.2*xstmp[2]);
   gMinuit->SetLimitedVariable(3, "lmass", xstmp[3], 1, 0.8*xstmp[3], 1.2*xstmp[3]);

   gMinuit->Minimize();

   const double *xs = gMinuit->X();
   gnorm1 = xs[0];
   gnorm2 = xs[1];
   lx = xs[2];
   lmass = xs[3];

   return;
}

double Shapes::GPm2ll( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt: " << x[1] << ", " << x[2] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];
   double m2llsig, m2llbkg;
   TrainGP( *hists_train_, m2llsig, m2llbkg );

   return m2llsig;
}

double Shapes::GPm2llX( const double *x ){
   cout << "gnorm: " << x[0] << ", " << x[1] << endl;
   cout << "lx, lmt: " << x[2] << ", " << x[3] << endl;

   gnorm1 = x[0];
   gnorm2 = x[1];
   lx = x[2];
   lmass = x[3];

   int nmasses = 8;
   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};

   // histograms
   vector<TH1D*> hgp_sig;
   for(int j=0; j < nmasses; j++){

      stringstream ssmass;
      ssmass << floor(masspnts[j]);
      string smass = ssmass.str();

      // signal shape
      hgp_sig.push_back( (TH1D*)(*hists_train_)["mbl"]["ttbar"+smass+"_signal"]
            ->Clone( ("hgp_sigx"+smass).c_str()) );
      hgp_sig[j]->Scale( 1.0/hgp_sig[j]->Integral("width") );

   }

   // set up training vector for cross validation
   int ntrain = 100;
   double rtrain = 300;

   int nval = 3;
   double m2ll_tot = 0;
   for(int c=0; c < nval; c++){ // n-fold cross validation

      ptrain.clear();
      vector<double> ptrainX;
      for(int i=0; i < ntrain; i++){ // exclude every nth point
         if( i%nval != c ){
            ptrain.push_back( (i+0.5)*rtrain/ntrain );
         }else{
            ptrainX.push_back( (i+0.5)*rtrain/ntrain );
         }
      }
      double m2llsig, m2llbkg;
      TrainGP( *hists_train_, m2llsig, m2llbkg );

      // evaluate shape at excluded points
      double m2ll = 0;
      for(unsigned int i=0; i < ptrainX.size(); i++){
         for(int j=0; j < nmasses; j++){

            double mean = Fmbl_gp(ptrainX[i],masspnts[j],"sig");
            double yi = hgp_sig[j]->GetBinContent( hgp_sig[j]->FindBin(ptrainX[i]) );
            double var = -1;
            if( do_gpvar ) var = Fmbl_gp_var(ptrainX[i],masspnts[j],"sig");
            else var = pow(hgp_sig[j]->GetBinError( hgp_sig[j]->FindBin(ptrainX[i]) ), 2);

            if( var <= 0 ){
               //cout << "NEGATIVE VARIANCE IN GP: " << var << "  ----> setting to minimum" << endl;
               var = 1E-10;
            }

            m2ll += log(var) + pow(mean-yi,2)/var + log(2*TMath::Pi());

         }
      }

      m2ll_tot += m2ll;
   }

   return m2ll_tot;

}
