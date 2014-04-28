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
<<<<<<< HEAD
#include <sstream>
=======
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <getopt.h>
>>>>>>> 06531be

using namespace std;

void print_usage(){

   cout << "\nUsage: ./DoFit <flags>\n";
   cout << "Flags: \n";
   cout << setiosflags(ios::left);
   cout << setw(25) << "\t-f --fit" << "Turn on fit.\n";
   cout << setw(25) << "\t-d --diagnostics" << "Turn on diagnostics.\n";
   cout << setw(25) << "\t-a --data" << "Run the fit on data (use full mc for training).\n";
   cout << setw(25) << "\t-m --masspnt  <value>" << "If running on mc, use masspoint indicated.\n";
   cout << setw(25) << "\t-h --help" << "Display this menu.\n";
   cout << endl;
   return;
}

int main(int argc, char* argv[]){

   // declarations
   Fitter fitter;
   map<string, Dataset> datasets;
   vector<Event> eventvec_datamc;
   vector<Event> eventvec_train;
   vector<Event> eventvec_test;
   map< string, map<string, TH1D*> > hists_all_;
   map< string, map<string, TH1D*> > hists_train_;
   map< string, map<string, TH1D*> > hists_test_;

   // output fit results
   int fitstatus=-1;
   double mt=0, mt_err=0;
   double kmbl=0, kmbl_err=0;
   double mcmass=0;
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);
   tree->Branch("mcmass", &mcmass);


   // option flags
   int c;
   int do_fit = 0;
   int do_diagnostics = 0;
   int use_data = 0;
   float masspnt = 0;

   struct option longopts[] = {
      { "fit",          no_argument,         &do_fit,          'f' },
      { "diagnostics",  no_argument,         &do_diagnostics,  'd' },
      { "data",         no_argument,         &use_data,        'a' },
      { "masspnt",      required_argument,   0,                'm' },
      { "profile",      no_argument,         0,                'p' },
      { "help",         no_argument,         NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdahpm:", longopts, NULL)) != -1 ) {
      switch(c)
      {
         case 'f' :
            do_fit = 1;
            break;

         case 'd' :
            do_diagnostics = 1;
            break;

         case 'a' :
            use_data = 1;
            break;

         case 'm' :
            masspnt = atof(optarg);
            break;

         case 'p' :
            fitter.compute_profile = true;
            break;

         case 'h' :
            print_usage();
            return -1;
            break;

         case 0:     /* getopt_long() set a variable, just keep going */
            break;

         case ':':   /* missing option argument */
            fprintf(stderr, "%s: option `-%c' requires an argument\n",
                  argv[0], optopt);
            return -1;

         case '?':
         default:    /* invalid option */
            fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                  argv[0], optopt);
            print_usage();
            return -1;

      }
   }


   fitter.LoadDatasets( datasets );

   // for event counting
   map<string, int> datacount;

   cout << "\nLoading datasets" << endl;
   for(map<string, Dataset>::iterator it = datasets.begin(); it != datasets.end(); it++){

      string name = it->first;
      Dataset *dat = &(it->second);

      datacount[name] = 0;
      datacount[name+"_bkgcontrol"] = 0;

      TFile file( (dat->path+dat->file).c_str() );
      TTree *trees = (TTree*)file.Get("RealData");
      TTree *treeb = (TTree*)file.Get("buBkg");

      cout << setiosflags(ios::left);
      cout << "... " << setw(25) << name
         << ": " << trees->GetEntries() << " events" << endl;
      cout << "... " << setw(25) << name+"_bkgcontrol"
         << ": " << treeb->GetEntries() << " events" << endl;

      if( do_diagnostics or use_data ){
         fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
               "RealData", eventvec_datamc );
         // bkg control sample
         fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
               "buBkg", eventvec_datamc );
      }

      // events for training and testing
      if( name.compare("data") != 0 ){

         if( use_data ){ // train on full mc set
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_train );
            fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
                  "buBkg", eventvec_train );
         }else{
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_train, 1 );
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_test, 2 );
            // bkg control sample
            fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
                  "buBkg", eventvec_train, 1 );
            fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
                  "buBkg", eventvec_test, 2 );
         }

      }

   }

   // data/mc plots, kinematic distributions
   fitter.GetVariables( eventvec_datamc );
   fitter.GetVariables( eventvec_train );
   fitter.GetVariables( eventvec_test );

   fitter.DeclareHists( hists_train_, "train" );
   fitter.FillHists( hists_train_, eventvec_train );

   //fitter.DeclareHists( hists_test_, "test" );
   //fitter.FillHists( hists_test_, eventvec_test );

   if( do_diagnostics ){ 
      fitter.DeclareHists( hists_all_, "all" );
      fitter.FillHists( hists_all_, eventvec_datamc );
      fitter.PrintHists( hists_all_ );
   }

   fitter.PlotTemplates();
/*
   return 0;

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

   TH1D *htemp = (TH1D*)fitter.hists_["mbl"]["ttbar161_signal"]->Clone("htemp");
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

   for(int i=0; i < ntrain; i++) cout << i << " " << a[i] << endl;

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
   */
   if( do_fit ){

<<<<<<< HEAD
=======
      vector<Event> eventvec_fit;
      vector<Event> eventvec_fit_bkgcontrol;
      map< string, map<string, TH1D*> > hists_fit_;
      map< string, map<string, TH1D*> > hists_fit_bkgcontrol_;

      if( use_data ){
         for(vector<Event>::iterator ev = eventvec_datamc.begin(); ev < eventvec_datamc.end();ev++){
            if( ev->type.find("data") != string::npos ){
               if( ev->type.find("bkgcontrol") == string::npos ){
                  eventvec_fit.push_back(*ev);
               }else{
                  eventvec_fit_bkgcontrol.push_back(*ev);
               }
            }
         }

         // flag events to be fitted
         for( vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
            ev->fit_event = true;
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
            }
         }
         for( vector<Event>::iterator ev = eventvec_fit_bkgcontrol.begin();
               ev < eventvec_fit_bkgcontrol.end(); ev++){
            ev->fit_event = true;
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
            }
         }
         cout << "Fit event count: " << endl;
         cout << setiosflags(ios::left);
         for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
            cout << "... " << setw(25) << it->first << ": " << it->second << " events" << endl;
         }

         fitter.DeclareHists( hists_fit_, "fit" );
         fitter.FillHists( hists_fit_, eventvec_fit, true );
         fitter.DeclareHists( hists_fit_bkgcontrol_, "fit_bkgcontrol" );
         fitter.FillHists( hists_fit_bkgcontrol_, eventvec_fit_bkgcontrol, true );

         double wgt=0;
         for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
            wgt += ev->weight;
         }
         cout << "wgt = " << wgt << endl;

         // do GP training
         Shapes * fptr = new Shapes( hists_fit_bkgcontrol_["mbl"]["fitevts"] );
         fptr->TrainGP( hists_train_ );
         fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
         fitter.aGPsig = fptr->aGPsig;
         fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
         fitter.aGPbkg = fptr->aGPbkg;

         /*
         TFile *ftemp = new TFile("ftemp.root","RECREATE");
         ftemp->cd();
         TH1D *htemp = new TH1D("hmbl","hmbl",100,0,300);
         for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
            for(int m=0; m < ev->mbls.size(); m++){
               htemp->Fill(ev->mbls[m]);
            }
         }
         htemp->Write("hmbl");
         ftemp->Close();
         return 0;
         */

         fitter.PlotTemplates( hists_train_ );
         return 0;

         // events for fitting, hists for training
         fitter.RunMinimizer( eventvec_fit, hists_fit_bkgcontrol_["mbl"]["fitevts"] );
         fitter.PlotResults( hists_fit_ ); // plot fitted events

         // fill results tree
         mcmass = 0;
         fitstatus = fitter.gMinuit->Status();
         const double *par = fitter.gMinuit->X();
         const double *par_err = fitter.gMinuit->Errors();
         mt = par[0];
         kmbl = par[1];
         mt_err = par_err[0];
         kmbl_err = par_err[1];

         tree->Fill();

         eventvec_fit.clear();
         eventvec_fit_bkgcontrol.clear();
         fitter.DeleteHists( hists_fit_ );
         fitter.DeleteHists( hists_fit_bkgcontrol_ );

      }else{ // loop over mc masses

         double masspnts [] = {161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5};
         for(int i=0; i < 8; i++){

            double mass = masspnts[i];
            // masspoint from command line
            if( masspnt != 0 ){
               bool check = false;
               for(int j=0; j < 8; j++){
                  if( masspnts[j] == masspnt ) check = true;
               }
               if(check){
                  mass = masspnt;
                  i = 7;
               }else{
                  cout << "masspoint " << masspnt << " not found!" << endl;
                  return -1;
               }
            }

            cout << "############# Fitting Masspoint " << mass << " ###############" << endl;

            stringstream dstr;
            dstr << floor(mass);
            string dname = "ttbar"+dstr.str();

            // load events to be fitted
            for(vector<Event>::iterator ev = eventvec_test.begin(); ev < eventvec_test.end(); ev++){
               if( ev->type.find(dname) != string::npos or ev->type.find("other") != string::npos ){
                  if( ev->type.find("bkgcontrol") == string::npos ){
                        eventvec_fit.push_back(*ev);
                  }else{
                     eventvec_fit_bkgcontrol.push_back(*ev);
                  }
               }
            }

            // define mc event weights -- ttbar events must have an average weight of 1
            double weight_norm = 0;
            double nevts_ttbar = 0;
            for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
               if( ev->type.find(dname) != string::npos ){
                  weight_norm += ev->weight;
                  nevts_ttbar++;
               }
            }
            cout << "nevts = " << nevts_ttbar << " weight_norm = " << weight_norm << endl;
            // reweight
            double wgt=0;
            for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
               ev->weight *= 1.0*nevts_ttbar/weight_norm;
               wgt += ev->weight;
            }
            cout << "wgt = " << wgt << endl;

            // flag events to be fitted
            for( vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
               ev->fit_event = true;
               for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
                  if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
               }
            }
            for( vector<Event>::iterator ev = eventvec_fit_bkgcontrol.begin();
                  ev < eventvec_fit_bkgcontrol.end(); ev++){
               ev->fit_event = true;
               for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
                  if( ev->process.compare(it->first) == 0 ) datacount[it->first]+=1;
               }
            }

            cout << "Fit event count: " << endl;
            cout << setiosflags(ios::left);
            for(map<string, int>::iterator it = datacount.begin(); it != datacount.end(); it++){
               cout << "... " << setw(25) << it->first << ": " << it->second << " events" << endl;
            }

            fitter.DeclareHists( hists_fit_, "fit" );
            fitter.FillHists( hists_fit_, eventvec_fit, true );
            fitter.DeclareHists( hists_fit_bkgcontrol_, "fit_bkgcontrol" );
            fitter.FillHists( hists_fit_bkgcontrol_, eventvec_fit_bkgcontrol, true );

            // do GP training
            Shapes * fptr = new Shapes( hists_fit_bkgcontrol_["mbl"]["fitevts"] );
            fptr->TrainGP( hists_train_ );
            fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            fitter.aGPsig = fptr->aGPsig;
            fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            fitter.aGPbkg = fptr->aGPbkg;

            /*
            TFile *ftemp = new TFile("ftemp.root","RECREATE");
            ftemp->cd();
            TH1D *htemp = new TH1D("hmbl","hmbl",100,0,300);
            for(vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
               for(int m=0; m < ev->mbls.size(); m++){
                  htemp->Fill(ev->mbls[m],ev->weight);
               }
            }
            htemp->Write("hmbl");
            ftemp->Close();
            return 0;
            */
            

            // events for fitting, hists for training
            fitter.RunMinimizer( eventvec_fit, hists_fit_bkgcontrol_["mbl"]["fitevts"] );
            fitter.PlotResults( hists_fit_ ); // plot fitted events

            fitter.PlotTemplates( hists_train_ );

            // fill results tree
            mcmass = mass;
            fitstatus = fitter.gMinuit->Status();
            const double *par = fitter.gMinuit->X();
            const double *par_err = fitter.gMinuit->Errors();
            mt = par[0];
            kmbl = par[1];
            mt_err = par_err[0];
            kmbl_err = par_err[1];

            tree->Fill();

            eventvec_fit.clear();
            eventvec_fit_bkgcontrol.clear();
            fitter.DeleteHists( hists_fit_ );
            fitter.DeleteHists( hists_fit_bkgcontrol_ );

         }

      }

   }
>>>>>>> 06531be


   //
   // write fit results
   //
   if( do_fit ){
      // set up output file path
      std::string pathstr;
      char* path = std::getenv("WORKING_DIR");
      if (path==NULL) {
         pathstr = "./results";
      }else {
         pathstr = path;
      }

      TFile *file = new TFile((pathstr+"/fitresults.root").c_str(), "RECREATE");
      file->cd();
      tree->Write();
      file->Write();
      file->Close();
   }

   return 0;
}
