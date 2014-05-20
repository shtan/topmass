#include "TopMass.h"
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <getopt.h>

using namespace std;

void print_usage(){

   cout << "\nUsage: ./DoFit <flags>\n";
   cout << "Flags: \n";
   cout << setiosflags(ios::left);
   cout << setw(25) << "\t-f --fit" << "Turn on fit.\n";
   cout << setw(25) << "\t-d --diagnostics" << "Turn on diagnostics.\n";
   cout << setw(25) << "\t-a --data" << "Run the fit on data (use full mc for training).\n";
   cout << setw(25) << "\t-p --profile" << "Run the likelihood profile.\n";
   cout << setw(25) << "\t-m --masspnt  <value>" << "If running on mc, use masspoint indicated.\n";
   cout << setw(25) << "\t-b --lmbl" << "Set the mbl lengthscale.\n";
   cout << setw(25) << "\t-t --lmbl" << "Set the mt lengthscale.\n";
   cout << setw(25) << "\t-l --lbnd" << "Left bound for exclusion.\n";
   cout << setw(25) << "\t-r --rbnd" << "Right bound for exclusion.\n";
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
   double fitchi2=0;
   double gplength_mbl=0;
   double gplength_mt=0;
   double lbound_mbl=0;
   double rbound_mbl=0;
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);
   tree->Branch("mcmass", &mcmass);
   tree->Branch("fitchi2", &fitchi2);
   tree->Branch("gplength_mbl", &gplength_mbl);
   tree->Branch("gplength_mt", &gplength_mt);
   tree->Branch("lbound_mbl", &lbound_mbl);
   tree->Branch("rbound_mbl", &rbound_mbl);


   // option flags
   int c;
   int do_fit = 0;
   int do_diagnostics = 0;
   int use_data = 0;
   float masspnt = 0;
   float lengthscale_mbl = 13;
   float lengthscale_mt = 32;

   struct option longopts[] = {
      { "fit",          no_argument,         &do_fit,          'f' },
      { "diagnostics",  no_argument,         &do_diagnostics,  'd' },
      { "data",         no_argument,         &use_data,        'a' },
      { "profile",      no_argument,         0,                'p' },
      { "masspnt",      required_argument,   0,                'm' },
      { "lmbl",         required_argument,   0,                'b' },
      { "lmt",          required_argument,   0,                't' },
      { "lbnd",         required_argument,   0,                'l' },
      { "rbnd",         required_argument,   0,                'r' },
      { "help",         no_argument,         NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdahpm:b:t:", longopts, NULL)) != -1 ) {
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

         case 'b' :
            lengthscale_mbl = atof(optarg);
            break;

         case 't' :
            lengthscale_mt = atof(optarg);
            break;

         case 'l' :
            fitter.lbnd = atof(optarg);
            break;

         case 'r' :
            fitter.rbnd = atof(optarg);
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

   fitter.gplength_mbl = lengthscale_mbl;
   fitter.gplength_mt = lengthscale_mt;

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

   fitter.DeclareHists( hists_test_, "test" );
   fitter.FillHists( hists_test_, eventvec_test );

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
         cout << "Training GP." << endl;
         Shapes * fptr = new Shapes( hists_fit_bkgcontrol_["mbl"]["fitevts"],
              fitter.gplength_mbl, fitter.gplength_mt, fitter.lbnd, fitter.rbnd );
         fptr->TrainGP( hists_train_ );
         fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
         fitter.aGPsig = fptr->aGPsig;
         fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
         fitter.aGPbkg = fptr->aGPbkg;

         fitter.PlotTemplates( hists_train_ );

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
         fitchi2 = fitter.fitchi2;

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
            Shapes * fptr = new Shapes( hists_fit_bkgcontrol_["mbl"]["fitevts"],
                 fitter.gplength_mbl, fitter.gplength_mt, fitter.lbnd, fitter.rbnd );
            fptr->TrainGP( hists_train_ );
            fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            fitter.aGPsig = fptr->aGPsig;
            fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            fitter.aGPbkg = fptr->aGPbkg;

            typedef map<string, TH1D*> tmap;
            typedef map<string, tmap> hmap;
            for( hmap::iterator h = hists_train_.begin(); h != hists_train_.end(); h++){
               for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
                  for(int n=1; n < t->second->GetNbinsX(); n++){
                     if( t->second->GetBinContent(n) == 0 )
                        t->second->SetBinError(n, 1.0/35000);
                  }
               }
            }
            fitter.PlotTemplates( hists_train_ );

            // events for fitting, hists for training
            fitter.RunMinimizer( eventvec_fit, hists_fit_bkgcontrol_["mbl"]["fitevts"] );
            fitter.PlotResults( hists_fit_ ); // plot fitted events

            cout << "Fit Chi2 = " << fitter.fitchi2 << endl;


            // fill results tree
            mcmass = mass;
            fitstatus = fitter.gMinuit->Status();
            const double *par = fitter.gMinuit->X();
            const double *par_err = fitter.gMinuit->Errors();
            mt = par[0];
            kmbl = par[1];
            mt_err = par_err[0];
            kmbl_err = par_err[1];
            fitchi2 = fitter.fitchi2;
            gplength_mbl = lengthscale_mbl;
            gplength_mt = lengthscale_mt;
            lbound_mbl = fitter.lbnd;
            rbound_mbl = fitter.rbnd;

            tree->Fill();

            eventvec_fit.clear();
            eventvec_fit_bkgcontrol.clear();
            fitter.DeleteHists( hists_fit_ );
            fitter.DeleteHists( hists_fit_bkgcontrol_ );

         }

      }

   }


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
