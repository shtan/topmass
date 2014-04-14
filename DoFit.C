#include "TopMass.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <vector>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]){

   // declarations
   Fitter fitter;
   map<string, Dataset> datasets;
   vector<Event> eventvec;

   // output fit results
   int fitstatus=-1;
   double mt=0, mt_err=0;
   double kmbl=0, kmbl_err=0;
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);


   // option flags
   char c;
   bool do_fit = true;
   bool do_diagnostics = true;

   while( (c = getopt(argc, argv, "fdh")) != -1 ) {
      switch(c)
      {
         case 'f' :
            do_fit = false;
            break;

         case 'd' :
            do_diagnostics = false;
            break;

         case 'h' :
            cout << "Usage: ./DoFit <flags>\n";
            cout << "Flags: \n";
            cout << "\t-f\t          Turn off fit.\n";
            cout << "\t-d\t          Turn off diagnostics.\n";
            cout << "\t-h\t          Display this menu.\n";
            return -1;
            break;

         default : 
            continue;
      }
   }

   fitter.LoadDatasets( datasets );

   cout << "\nLoading datasets" << endl;
   for(map<string, Dataset>::iterator it = datasets.begin(); it != datasets.end(); it++){

      string name = it->first;
      Dataset *dat = &(it->second);

      fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
            "RealData", eventvec );
      if( name.compare("data") == 0 ){ // bkg control sample
         fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
               "buBkg", eventvec );
      }

   }

   fitter.GetVariables( eventvec );
   cout << endl;

   fitter.DeclareHists();
   fitter.FillHists( eventvec );
   if( do_diagnostics ){
      fitter.PrintHists();
   }
   cout << endl;

   fitter.PlotTemplates();

   //
   // ################ GP regression ###################
   //
   TFile *filegp = new TFile("testgp.root", "RECREATE");
   filegp->cd();

   int ntrain = fitter.hists_["mbl"]["ttbar172_signal"]->GetNbinsX();
   TH2D *hmbl = new TH2D("hmbl","hmbl",ntrain,0,300,20,160,183);

   // get training points
   vector<double> px; 
   vector< vector<double> > pz;
   double masspnts [] = {161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5};
   for(int i=0; i < ntrain; i++){
      vector<double> pmass;
      for(int j=0; j < 8; j++){

         stringstream dstr;
         dstr << floor(masspnts[j]);
         string dname = "ttbar"+dstr.str();

         TH1D *htemp = (TH1D*)fitter.hists_["mbl"][dname+"_signal"]->Clone("htemp");
         htemp->Scale( 1.0/htemp->Integral("width") );

         double bincent = htemp->GetBinCenter(i+1);
         double bincont = htemp->GetBinContent(i+1);
         if( j==0 ) px.push_back( bincent );
         pmass.push_back( bincont );

         hmbl->SetBinContent( hmbl->FindBin(bincent,masspnts[j]), bincont );

      }
      pz.push_back( pmass );
   }

   TCanvas *ctest = new TCanvas("ctest","ctest",800,800);
   ctest->cd();
   hmbl->Draw("LEGO2Z");
   ctest->Write("cmbl2d");

   //
   // test 1D case
   //
   double lgp = 25.0;
   ntrain = 100;
   double rtrain = 300;

   vector<double> ptrain;
   for(int i=0; i < ntrain; i++) ptrain.push_back( (i+0.5)*rtrain/ntrain );

   TH1D *htemp = (TH1D*)fitter.hists_["mbl"]["ttbar172_signal"]->Clone("htemp");
   //htemp->Scale( 1.0/htemp->Integral("width") );
   cout << "norm = " << htemp->Integral("width") << endl;

   // compute covariance matrix
   TMatrixD K(ntrain,ntrain);
   for(int i=0; i < ntrain; i++){
      for(int j=0; j < ntrain; j++){
         K[i][j] = exp(-0.5*pow( (ptrain[i]-ptrain[j])/lgp, 2));
     }
   }
   // compute noise matrix
   TMatrixD N(ntrain,ntrain);
   for(int i=0; i < ntrain; i++){
      double binerr = htemp->GetBinError( htemp->FindBin(ptrain[i]) );
      cout << htemp->GetBinCenter( htemp->FindBin(ptrain[i]) ) << ": " << binerr << endl;
      for(int j=0; j < ntrain; j++){
         if( i==j ){
            N[i][j] = pow( max(binerr,2E-05), 2 );
         }else{
            N[i][j] = 0;
         }
     }
   }
   // inverse of sum
   TMatrixD A(ntrain,ntrain);
   A = K + N;
   A.Invert();

   // vector of training points
   TVectorD y(ntrain);
   for(int i=0; i < ntrain; i++) y[i] = htemp->GetBinContent( htemp->FindBin(ptrain[i]) );
   // alpha vector
   TVectorD a(ntrain);
   a = A*y;

   // cout training points
   for(int i=0; i < ntrain; i++){
      cout << "(x,y,a) = " << ptrain[i] << ", " << y[i] << ", " << a[i] << endl;
   }

   // graph of training pnts
   TGraph *gtrain = new TGraph();
   for(int i=0; i < ntrain; i++) gtrain->SetPoint(i, ptrain[i], y[i]);

   // evaluate GP at 1k points
   TGraph *gp = new TGraph();
   for(int i=0; i < 1000; i++){
      double xgp = i*rtrain/1000;
      double fgp = 0;
      for(int k=0; k < ntrain; k++) fgp += a[k]*exp(-0.5*pow( (xgp-ptrain[k])/lgp, 2));
      gp->SetPoint(i, xgp, fgp);
   }

   TCanvas *hmblgp = new TCanvas("hmblgp","hmblgp",800,800);
   hmblgp->cd();

   htemp->SetMarkerColor(2);
   htemp->Draw();
   gtrain->SetMarkerStyle(3);
   gtrain->SetMarkerColor(4);
   gtrain->Draw("P");
   gp->Draw("C");
   hmblgp->Write("hmblgp");

   filegp->Close();

   return 0;
   if( do_fit ){
      fitter.RunMinimizer( eventvec );
      fitter.PlotResults();
   }



   //
   // write fit results
   //
   if( do_fit ){
      fitstatus = fitter.gMinuit->Status();
      const double *par = fitter.gMinuit->X();
      const double *par_err = fitter.gMinuit->Errors();
      mt = par[0];
      kmbl = par[1];
      mt_err = par_err[0];
      kmbl_err = par_err[1];

      tree->Fill();

      // set up output file path
      std::string pathstr;
      char* path = std::getenv("WORKING_DIR");
      if (path==NULL) {
         pathstr = "./results/";
      }else {
         pathstr = path;
      }

      TFile *file = new TFile((pathstr+"fitresults.root").c_str(), "RECREATE");
      file->cd();
      tree->Write();
      file->Write();
      file->Close();
   }

   return 0;
}
