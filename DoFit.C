#include "TopMass.h"
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <getopt.h>

using namespace std;

void print_usage(){

   cout << "\nUsage: ./DoFit <flags>\n";
   cout << "Flags: \n";
   cout << "\t-f --fit\t            Turn on fit.\n";
   cout << "\t-d --diagnostics\t    Turn on diagnostics.\n";
   cout << "\t-a --data\t           Run the fit on data (use full mc for training).\n";
   cout << "\t-h --help\t            Display this menu.\n";
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
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);


   // option flags
   int c;
   int do_fit = 0;
   int do_diagnostics = 0;
   int use_data = 0;

   struct option longopts[] = {
      { "fit",          no_argument,   &do_fit,          'f' },
      { "diagnostics",  no_argument,   &do_diagnostics,  'd' },
      { "data",         no_argument,   &use_data,        'a' },
      { "help",         no_argument,   NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdah", longopts, NULL)) != -1 ) {
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

   cout << "\nLoading datasets" << endl;
   for(map<string, Dataset>::iterator it = datasets.begin(); it != datasets.end(); it++){

      string name = it->first;
      Dataset *dat = &(it->second);

      if( do_diagnostics or use_data ){
         fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
               "RealData", eventvec_datamc );
         if( name.compare("data") == 0 ){ // bkg control sample
            fitter.ReadNtuple( dat->path+dat->file, name+"_bkgcontrol", dat->mc_xsec/dat->mc_nevts,
                  "buBkg", eventvec_datamc );
         }
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

   if( do_fit ){

      vector<Event> eventvec_fit;
      vector<Event> eventvec_fit_bkgcontrol;
      map< string, map<string, TH1D*> > hists_fit_;
      map< string, map<string, TH1D*> > hists_fit_bkgcontrol_;

      if( use_data ){
         for(vector<Event>::iterator ev = eventvec_datamc.begin(); ev < eventvec_datamc.end(); ev++){
            if( ev->type.find("data") != string::npos ){
               if( ev->type.find("bkgcontrol") == string::npos ){
                  eventvec_fit.push_back(*ev);
               }else{
                  eventvec_fit_bkgcontrol.push_back(*ev);
               }
            }
         }
      }else{
         for( vector<Event>::iterator ev = eventvec_test.begin(); ev < eventvec_test.end(); ev++){
            if(ev->type.find("ttbar181") != string::npos or ev->type.find("other") != string::npos ){
               if( ev->type.find("bkgcontrol") == string::npos ){
                  if( ev->type.find("signal") != string::npos )
                     eventvec_fit.push_back(*ev);
               }else{
                  eventvec_fit_bkgcontrol.push_back(*ev);
               }
            }
         }
      }
      for( vector<Event>::iterator ev = eventvec_fit.begin(); ev < eventvec_fit.end(); ev++){
         ev->fit_event = true;
      }
      for( vector<Event>::iterator ev = eventvec_fit_bkgcontrol.begin();
            ev < eventvec_fit_bkgcontrol.end(); ev++){
         ev->fit_event = true;
      }

      fitter.DeclareHists( hists_fit_, "fit" );
      fitter.FillHists( hists_fit_, eventvec_fit, true );
      fitter.DeclareHists( hists_fit_bkgcontrol_, "fit_bkgcontrol" );
      fitter.FillHists( hists_fit_bkgcontrol_, eventvec_fit_bkgcontrol, true );

      // do GP training
      Shapes * fptr = new Shapes( hists_fit_bkgcontrol_["mbl_fit"]["fitevts"] );
      fptr->TrainGP( hists_train_ );
      fitter.aGP.ResizeTo( fptr->aGP.GetNoElements() );
      fitter.aGP = fptr->aGP;

      // events for fitting, hists for training
      fitter.RunMinimizer( eventvec_fit, hists_fit_bkgcontrol_["mbl_fit"]["fitevts"] );
      fitter.PlotResults( hists_fit_ ); // plot fitted events

      fitter.PlotTemplates( hists_train_ );
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
