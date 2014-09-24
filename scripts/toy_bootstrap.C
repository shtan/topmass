#include "TH1.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLegend.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

using namespace std;

// random number engine
ROOT::Math::Random<ROOT::Math::GSLRngMT> GSLr;

// global event vectors
vector<double> vx1, vx2;
vector<int> fitevents;

int numevents = 1E5;
int numPE = 1E3;
int numBPE = 100;

//
// set variable properties
//
// sigma
double sx1 = 10;
double sx2 = 10;
// mean
double mx1 = 100;
double mx2 = 100;
// rho
double rho12 = 0.0;

//
// fit results
//
double X = 0.0;
double Xerr = 0.0;
double Xup = 0.0;
double Xdn = 0.0;

// measurement results
double X1Dmean = 0.0;
double X1Drms = 0.0;
double X2Dmean = 0.0;
double X2Drms = 0.0;
double XI1Dmean = 0.0;
double XI1Drms = 0.0;
double XI2Dmean = 0.0;
double XI2Drms = 0.0;
double XB1Dmean = 0.0;
double XB1Drms = 0.0;
double XB2Dmean = 0.0;
double XB2Drms = 0.0;

double Gaus1D( double x, double m, double s );
double Gaus2D( double x1, double m1, double s1, double x2, double m2, double s2, double rho12 );

double m2ll_1d( const double *xx );
double m2ll_2d( const double *xx );

void Fit1D();
void Fit2D();

void run_measurement();

void toy_bootstrap(){

   // random number seed
   GSLr.SetSeed(1);

   // run measurements, fill graphs
   TGraph *g1D = new TGraph();
   TGraphErrors *g2D = new TGraphErrors();
   TGraph *g1DIndPE = new TGraph();
   TGraph *g2DIndPE = new TGraph();
   TGraph *g1DBPE = new TGraph();
   TGraph *g2DBPE = new TGraph();

   TGraph *g1Dm = new TGraph();
   TGraphErrors *g2Dm = new TGraphErrors();
   TGraph *g1DIndPEm = new TGraph();
   TGraph *g2DIndPEm = new TGraph();
   TGraph *g1DBPEm = new TGraph();
   TGraph *g2DBPEm = new TGraph();

   for(int i=0; i < 10; i++){

      rho12 = 1.0*i/10;
      cout << "#################### rho = " << rho12 << " ####################" << endl;

      run_measurement();

      g2D->SetPoint(i, rho12+0.05, X2Drms);
      g2D->SetPointError(i, 0.05, 0.0);

      g1D->SetPoint(i, rho12+0.05, X1Drms);
      g1DIndPE->SetPoint(i, rho12+0.02, XI1Drms);
      g2DIndPE->SetPoint(i, rho12+0.04, XI2Drms);
      g1DBPE->SetPoint(i, rho12+0.06, XB1Drms);
      g2DBPE->SetPoint(i, rho12+0.08, XB2Drms);

      g2Dm->SetPoint(i, rho12+0.05, X2Dmean);
      g2Dm->SetPointError(i, 0.05, 0.0);

      g1Dm->SetPoint(i, rho12+0.05, X1Dmean);
      g1DIndPEm->SetPoint(i, rho12+0.02, XI1Dmean);
      g2DIndPEm->SetPoint(i, rho12+0.04, XI2Dmean);
      g1DBPEm->SetPoint(i, rho12+0.06, XB1Dmean);
      g2DBPEm->SetPoint(i, rho12+0.08, XB2Dmean);

   }

   TCanvas *canvas = new TCanvas("canvas","canvas",800,800);
   canvas->cd();

   g2D->SetMarkerStyle(20);
   g1D->SetLineColor(2);
   g1D->SetMarkerColor(2);
   g1D->SetMarkerStyle(5);
   g1DIndPE->SetLineColor(3);
   g1DIndPE->SetMarkerColor(3);
   g1DIndPE->SetMarkerStyle(20);
   g2DIndPE->SetLineColor(4);
   g2DIndPE->SetMarkerColor(4);
   g2DIndPE->SetMarkerStyle(20);
   g1DBPE->SetLineColor(7);
   g1DBPE->SetMarkerColor(7);
   g1DBPE->SetMarkerStyle(20);
   g2DBPE->SetLineColor(6);
   g2DBPE->SetMarkerColor(6);
   g2DBPE->SetMarkerStyle(20);

   g2D->Draw("APE");
   g1D->Draw("P");
   g1DIndPE->Draw("P");
   g2DIndPE->Draw("P");
   g1DBPE->Draw("P");
   g2DBPE->Draw("P");

   TLegend * legend = new TLegend(0.717,0.650,0.874,0.870);
   legend->AddEntry( g2D, "2D fit MINUIT error", "p" );
   legend->AddEntry( g1D, "1D fit MINUIT error", "p" );
   legend->AddEntry( g2DIndPE, "2D fit, Ind. PE RMS", "p" );
   legend->AddEntry( g1DIndPE, "1D fit, Ind. PE RMS", "p" );
   legend->AddEntry( g2DBPE, "2D fit, Boot PE RMS", "p" );
   legend->AddEntry( g1DBPE, "1D fit, Boot PE RMS", "p" );
   legend->Draw("same");

   canvas->Draw();

   // plot means
   TCanvas *canvasm = new TCanvas("canvasm","canvasm",800,800);
   canvasm->cd();

   g2Dm->SetMarkerStyle(20);
   g1Dm->SetLineColor(2);
   g1Dm->SetMarkerColor(2);
   g1Dm->SetMarkerStyle(5);
   g1DIndPEm->SetLineColor(3);
   g1DIndPEm->SetMarkerColor(3);
   g1DIndPEm->SetMarkerStyle(20);
   g2DIndPEm->SetLineColor(4);
   g2DIndPEm->SetMarkerColor(4);
   g2DIndPEm->SetMarkerStyle(20);
   g1DBPEm->SetLineColor(7);
   g1DBPEm->SetMarkerColor(7);
   g1DBPEm->SetMarkerStyle(20);
   g2DBPEm->SetLineColor(6);
   g2DBPEm->SetMarkerColor(6);
   g2DBPEm->SetMarkerStyle(20);

   g2Dm->Draw("APE");
   g1Dm->Draw("P");
   g1DIndPEm->Draw("P");
   g2DIndPEm->Draw("P");
   g1DBPEm->Draw("P");
   g2DBPEm->Draw("P");

   TLegend * legendm = new TLegend(0.717,0.650,0.874,0.870);
   legendm->AddEntry( g2D, "2D fit MINUIT central value", "p" );
   legendm->AddEntry( g1D, "1D fit MINUIT central value", "p" );
   legendm->AddEntry( g2DIndPE, "2D fit, Ind. PE mean", "p" );
   legendm->AddEntry( g1DIndPE, "1D fit, Ind. PE mean", "p" );
   legendm->AddEntry( g2DBPE, "2D fit, Boot PE mean", "p" );
   legendm->AddEntry( g1DBPE, "1D fit, Boot PE mean", "p" );
   legendm->Draw("same");

   canvasm->Draw();

   return;
}

void run_measurement(){
   // event loop
   vx1.clear(); vx2.clear();
   for(int n=0; n < numevents; n++){
      double dx1=0, dx2=0;
      GSLr.Gaussian2D( sx1, sx2, rho12, dx1, dx2 );
      vx1.push_back( mx1 + dx1 );  vx2.push_back( mx2 + dx2 );
   }

   // entire dataset
   fitevents.clear();
   for(int i=0; i < numevents; i++) fitevents.push_back(i);
   Fit1D();
   cout << "========== 1D MINIMIZATION ==========" << endl;
   cout << "MINIMUM = " << X << " +- " << Xerr << endl;
   cout << "MINOS ERROR: " << Xdn << " + " << Xup << endl;
   X1Dmean = X;
   X1Drms = Xerr;

   Fit2D();
   cout << "========== 2D MINIMIZATION ==========" << endl;
   cout << "MINIMUM = " << X << " +- " << Xerr << endl;
   cout << "MINOS ERROR: " << Xdn << " + " << Xup << endl;
   X2Dmean = X;
   X2Drms = Xerr;

   // independent pseudoexperiments
   vector<double> X1D, X2D;
   for(int p=0; p < numPE; p++){
      fitevents.clear();
      for(int i=(1.0*p/numPE)*numevents; i < numevents*(1.0*(p+1)/numPE); i++){
         fitevents.push_back(i);
      }
      Fit1D(); X1D.push_back( X );
      Fit2D(); X2D.push_back( X );
   }

   XI1Dmean = TMath::Mean(X1D.begin(), X1D.end());
   XI1Drms = TMath::RMS(X1D.begin(), X1D.end())/sqrt(numPE);

   XI2Dmean = TMath::Mean(X2D.begin(), X2D.end());
   XI2Drms = TMath::RMS(X2D.begin(), X2D.end())/sqrt(numPE);

   cout << "========== 1D IND PEs ==========" << endl;
   cout << "MINIMUM = " << XI1Dmean << " +- " << XI1Drms << endl;
   cout << "========== 2D IND PEs ==========" << endl;
   cout << "MINIMUM = " << XI2Dmean << " +- " << XI2Drms << endl;

   // bootstrap pseudoexperiments
   vector<double> XB1D, XB2D;
   for(int p=0; p < numBPE; p++){
      fitevents.clear();
      for(int i=0; i < numevents; i++){
         fitevents.push_back(floor(GSLr.Rndm()*numevents));
      }
      Fit1D(); XB1D.push_back( X );
      Fit2D(); XB2D.push_back( X );
   }

   XB1Dmean = TMath::Mean(XB1D.begin(), XB1D.end());
   XB1Drms = TMath::RMS(XB1D.begin(), XB1D.end());

   XB2Dmean = TMath::Mean(XB2D.begin(), XB2D.end());
   XB2Drms = TMath::RMS(XB2D.begin(), XB2D.end());

   cout << "========== 1D BOOT PEs ==========" << endl;
   cout << "MINIMUM = " << XB1Dmean << " +- " << XB1Drms << endl;
   cout << "========== 2D BOOT PEs ==========" << endl;
   cout << "MINIMUM = " << XB2Dmean << " +- " << XB2Drms << endl;

   return;
}

void Fit1D(){

   // setup minimizer
   ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
   ROOT::Math::Functor f( &m2ll_1d, 1 );
   min.SetFunction( f );
   min.SetVariable(0, "topMass", 75.0, 1.0);

   // minimize
   min.Minimize();
   min.GetMinosError(0,Xdn,Xup);

   // get results
   const double *xs = min.X();
   const double *es = min.Errors();
   X = xs[0];
   Xerr = es[0];

   return;
}

void Fit2D(){
  
   // setup minimizer
   ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
   ROOT::Math::Functor f( &m2ll_2d, 1 );
   min.SetFunction( f );
   min.SetVariable(0, "topMass", 75.0, 1.0);

   // minimize
   min.Minimize();
   min.GetMinosError(0,Xdn,Xup);

   // get results
   const double *xs = min.X();
   const double *es = min.Errors();
   X = xs[0];
   Xerr = es[0];

   return;
}

double m2ll_1d( const double *xx ){

   double m2ll = 0;
   for(vector<int>::iterator ev = fitevents.begin(); ev != fitevents.end(); ev++){
      m2ll -= 2.0*log(Gaus1D( vx1[*ev], xx[0], sx1));
   }
   for(vector<int>::iterator ev = fitevents.begin(); ev != fitevents.end(); ev++){
      m2ll -= 2.0*log(Gaus1D( vx2[*ev], xx[0], sx2));
   }

   return m2ll;

}

double m2ll_2d( const double *xx ){

   double m2ll = 0;
   //for(unsigned int i=0; i < vx1.size(); i++){
   for(vector<int>::iterator ev = fitevents.begin(); ev != fitevents.end(); ev++){
      m2ll -= 2.0*log(Gaus2D( vx1[*ev], xx[0], sx1, vx2[*ev], xx[0], sx2, rho12));
   }

   return m2ll;

}

double Gaus1D( double x, double m, double s ){
   double norm = sqrt(2*TMath::Pi())*s;
   double exp = TMath::Exp( -(x-m)*(x-m)/(2*s*s) );
   return exp/norm;
}

double Gaus2D( double x1, double m1, double s1, double x2, double m2, double s2, double rho ){
   double norm = 2*TMath::Pi()*s1*s2*sqrt(1-rho*rho);
   double exp0 = -1.0/(2*(1-rho*rho));
   double exp1 = (x1-m1)*(x1-m1)/(s1*s1);
   double exp2 = (x2-m2)*(x2-m2)/(s2*s2);
   double exp3 = (-2*rho*(x1-m1)*(x2-m2))/(s1*s2);
   return TMath::Exp(exp0*(exp1+exp2+exp3))/norm;
}

