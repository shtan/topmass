#include "Shapes.h"

#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
   lmass = 3.0;
   gnorm = 1.0;
   gnoise = 1.0E-05;
   int ntrain = 100;
   //double rtrain = 1000;
   double rtrain = 300;
   ptrain.push_back(0);
   //for(int i=0; i < ntrain; i++) ptrain.push_back( (i+0.5)*rtrain/ntrain );
   double x = 1.5;
   while( x < 180 ){
      ptrain.push_back(x);
      x+=3;
   }
   ptrain.push_back(187.5);
   ptrain.push_back(202.5);
   ptrain.push_back(217.5);
   ptrain.push_back(232.5);
   ptrain.push_back(247.5);
   ptrain.push_back(262.5);
   ptrain.push_back(277.5);
   ptrain.push_back(292.5);
   //ptrain.push_back(195);
   //ptrain.push_back(225);
   //ptrain.push_back(255);
   //ptrain.push_back(285);

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

   //return norm*(k*Fmbl_sig_param(x, mt)/integral + (1-k)*Fmbl_bkg(x));
   //return norm*(k*Fmbl_sig_gp(x, mt)/integralsig + (1-k)*Fmbl_bkg(x));
   double val = norm*(k*Fmbl_sig_gp(x, mt)/integralsig + (1-k)*Fmbl_bkg_gp(x, mt)/integralbkg);
   if( val <= 0 ) return 1E-10;
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

double Shapes::Fmbl_sig_gp(double x, double mt){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;

   double fgp = 0;
   for(unsigned int i=0; i < ptrain.size(); i++){
      for(int j=0; j < nmasses; j++){
         fgp += aGPsig[i+j*ptrain.size()]
            *GPkern( x, ptrain[i], lmbl, mt, masspnts[j], lmass, gnorm, gnoise );
      }
   }

   return fgp;
}

double Shapes::Fmbl_bkg_gp(double x, double mt){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;

   double fgp = 0;
   for(unsigned int i=0; i < ptrain.size(); i++){
      for(int j=0; j < nmasses; j++){
         fgp += aGPbkg[i+j*ptrain.size()]
            *GPkern( x, ptrain[i], lmbl, mt, masspnts[j], lmass, gnorm, gnoise );
      }
   }

   return fgp;
}

double Shapes::Fmbl_bkg(double x){

   return hmbl_bkg->Interpolate(x) / norm_mbl_bkg;
}

double Shapes::GPkern(double x1, double x2, double lx, double m1, double m2, double lm,
      double norm, double noise){

   double delta = (x1==x2 and m1==m2) ? 1.0 : 0.0;
   double kernel = norm*exp(-0.5*(pow( (x1-x2)/lx, 2)+pow( (m1-m2)/lm, 2))) + noise*delta;
   return kernel;
}

void Shapes::TrainGP( map< string, map<string, TH1D*> > & hists_ ){

   double masspnts [] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5};
   int nmasses = 8;

   // histograms
   vector<TH1D*> hgp_sig;
   vector<TH1D*> hgp_bkg;
   for(int i=0; i < nmasses; i++){

      stringstream ssmass;
      ssmass << floor(masspnts[i]);
      string smass = ssmass.str();

      // TEMP modify sig shape
      /*
      int bin = hists_["mbl"]["ttbar"+smass+"_signal"]->FindBin(170.0);
      double val = hists_["mbl"]["ttbar"+smass+"_signal"]->GetBinContent(bin);
      for(int x=bin; x <= hists_["mbl"]["ttbar"+smass+"_signal"]->GetNbinsX(); x++)
         hists_["mbl"]["ttbar"+smass+"_signal"]->SetBinContent(x,val*exp(-0.15*(x-bin)));
         */
      //hists_["mbl"]["ttbar"+smass+"_signal"]->Rebin(10);
      // TEMP modify sig shape
      
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

   // define TGraphs for training (variable bin size)
   //for(int b=0; b <= hgp_sig[0]->GetNbinsX(); b++){
   //   cout << b << ": [" << hgp_sig[0]->GetBinLowEdge(b)
   //      << ", " << hgp_sig[0]->GetBinLowEdge(b)+hgp_sig[0]->GetBinWidth(b) << "]" << endl;
   //}
   
   vector<TGraphErrors*> ggp_sig;
   vector<TGraphErrors*> ggp_bkg;
   for(int sb=0; sb < 2; sb++){
      for(int i=0; i < nmasses; i++){

         // get histogram
         TH1D* htemp;
         if( sb==0 ) htemp = (TH1D*)hgp_sig[i]->Clone("htemp");
         if( sb==1 ) htemp = (TH1D*)hgp_bkg[i]->Clone("htemp");

         TGraphErrors *gtemp = new TGraphErrors();
         gtemp->SetPoint(0, 0.0, 0.0);
         gtemp->SetPointError(0, 0.0, htemp->GetBinError(1));
         int pnt = 1;
         for(int b=1; b <= htemp->GetNbinsX(); b++){
            if( htemp->GetBinCenter(b) < 180.0 ){
               gtemp->SetPoint(pnt, htemp->GetBinCenter(b), htemp->GetBinContent(b));
               gtemp->SetPointError(pnt, 0.0, htemp->GetBinError(b));
               pnt++;
               //cout << b << ": [" << htemp->GetBinLowEdge(b)
               //   << ", " << htemp->GetBinLowEdge(b)+htemp->GetBinWidth(b) << "]" << endl;
            }
         }
         // now take care of tail bins
         int firstbin = htemp->FindBin(180);
         int numbins = htemp->GetNbinsX() - firstbin + 1;
         int grpsize = numbins / 8;
         for(int j=0; j < 8; j++){
            double cent = 0;
            double val = 0;
            double errsq = 0;
            int begin = firstbin + j*grpsize;
            int count = 0;
            for(int b = begin; b < begin + grpsize; b++){
               cent += htemp->GetBinCenter(b);
               val += htemp->GetBinContent(b);
               errsq += pow(htemp->GetBinError(b),2);
               count++;
               cout << b << "(" << j << "): [" << htemp->GetBinLowEdge(b)
                  << ", " << htemp->GetBinLowEdge(b)+htemp->GetBinWidth(b) << "]" << endl;
            }
            cent /= count;
            val /= count;
            errsq /= count*count;
            cout << "***** (cent, val) = " << cent << ", " << val << endl;
            gtemp->SetPoint(pnt, cent, val);
            gtemp->SetPointError(pnt, 0.0, sqrt(errsq));
            pnt++;
         }

         if( sb==0 ) ggp_sig.push_back( gtemp );
         if( sb==1 ) ggp_bkg.push_back( gtemp );

         delete htemp;
      }
   }

   // set training point locations
   double xtemp, ytemp;
   for(int i=0; i < ggp_sig[0]->GetN(); i++){
      ggp_sig[0]->GetPoint(i,xtemp,ytemp);
      //ptrain.push_back( xtemp );
      if( ptrain[i] != xtemp ) cout << "ERROR IN PTRAIN PNTS" << endl;
   }
   int ntrain = ptrain.size();
   

   // compute covariance matrix
   TMatrixD K(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      for(int j=0; j < ntrain*nmasses; j++){
         int im = i % ntrain;
         int jm = j % ntrain;
         int imass = i / ntrain;
         int jmass = j / ntrain;
         K[i][j] = GPkern( ptrain[im], ptrain[jm], lmbl, masspnts[imass], masspnts[jmass], lmass,
              gnorm, gnoise );
     }
   }
   // compute noise matrix
   
   TMatrixD Nsig(ntrain*nmasses,ntrain*nmasses);
   TMatrixD Nbkg(ntrain*nmasses,ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      double binerr_sig = ggp_sig[imass]->GetErrorY( im );
      double binerr_bkg = ggp_bkg[imass]->GetErrorY( im );
      for(int j=0; j < ntrain*nmasses; j++){
         if( i==j ){
            Nsig[i][j] = binerr_sig*binerr_sig;//pow( max(binerr_sig,0.001), 2 );
            Nbkg[i][j] = binerr_bkg*binerr_bkg;//pow( max(binerr_bkg,0.001), 2 );
         }else{
            Nsig[i][j] = 0.0;
            Nbkg[i][j] = 0.0;
         }
     }
   }

   // inverse of sum
   TMatrixD Asig(ntrain*nmasses,ntrain*nmasses);
   TMatrixD Abkg(ntrain*nmasses,ntrain*nmasses);
   Asig = K + Nsig;
   Abkg = K + Nbkg;
   Asig.Invert();
   Abkg.Invert();

   // vector of training points
   TVectorD ysig(ntrain*nmasses);
   TVectorD ybkg(ntrain*nmasses);
   for(int i=0; i < ntrain*nmasses; i++){
      int im = i % ntrain;
      int imass = i / ntrain;
      ggp_sig[imass]->GetPoint(im, xtemp, ytemp);
      ysig[i] = ytemp;
      //cout << x << ", " << y << endl;
      ggp_bkg[imass]->GetPoint(im, xtemp, ytemp);
      ybkg[i] = ytemp;
   }

   // alpha vector
   aGPsig.ResizeTo( ntrain*nmasses );
   aGPsig = Asig*ysig;

   aGPbkg.ResizeTo( ntrain*nmasses );
   aGPbkg = Abkg*ybkg;

   //for(int i=0; i < nmasses; i++){
   //   delete ggp_sig[i];
   //   delete ggp_bkg[i];
   //}
   
}
