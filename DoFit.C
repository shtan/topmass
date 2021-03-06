#include "TopMass.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

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
   cout << setw(25) << "\t-n --run_number" << "Run number.\n";
   cout << setw(25) << "\t-f --fit" << "Turn on fit.\n";
   cout << setw(25) << "\t-d --diagnostics" << "Turn on diagnostics.\n";
   cout << setw(25) << "\t-e --templates" << "Turn on template plots.\n";
   cout << setw(25) << "\t-x --learnparams" << "Do a fit to learn the GP hyperparameters.\n";
   cout << setw(25) << "\t-a --data" << "Run the fit on data (use full mc for training).\n";
   cout << setw(25) << "\t-p --profile" << "Run the likelihood profile.\n";
   cout << setw(25) << "\t-m --masspnt  <value>" << "If running on mc, use masspoint indicated.\n";
   cout << setw(25) << "\t-o --bootstrap" << "Turn on bootstrapping.\n";
   cout << setw(25) << "\t-c --fracevts" << "Fit fraction of events.\n";
   cout << setw(25) << "\t-b --mbl" << "Activate Mbl distribution.\n";
   cout << setw(25) << "\t-t --mt2" << "Activate MT2 220 distribution.\n";
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
   int run_number=-1;
   int fitstatus=-1;
   double mt=0, mt_err=0;
   double kmbl=0, kmbl_err=0;
   double k220=0, k220_err=0;
   double mcmass=0;
   double fitchi2=0;
   double tsig_mbl_chi2 [8] = {0};
   double tbkg_mbl_chi2 [8] = {0};
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("runNumber", &run_number);
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);
   tree->Branch("k220", &k220);
   tree->Branch("k220_err", &k220_err);
   tree->Branch("mcmass", &mcmass);
   tree->Branch("fitchi2", &fitchi2);
   tree->Branch("tsig_mbl_chi2", tsig_mbl_chi2, "tsig_mbl_chi2[8]/D");
   tree->Branch("tbkg_mbl_chi2", tbkg_mbl_chi2, "tbkg_mbl_chi2[8]/D");


   // option flags
   int c;
   int do_fit = 0;
   int do_diagnostics = 0;
   int use_data = 0;
   float masspnt = 0;
   int do_bootstrap = 0;
   int do_templates = 0;
   int do_learnparams = 0;
   int do_mbl = 0;
   int do_mt2 = 0;
   double fracevts = -1;

   struct option longopts[] = {
      { "run_number",   required_argument,   0,                'n' },
      { "fit",          no_argument,         &do_fit,          'f' },
      { "diagnostics",  no_argument,         &do_diagnostics,  'd' },
      { "templates",    no_argument,         &do_templates,    'e' },
      { "learnparams",  no_argument,         &do_learnparams,  'x' },
      { "data",         no_argument,         &use_data,        'a' },
      { "profile",      no_argument,         0,                'p' },
      { "masspnt",      required_argument,   0,                'm' },
      { "bootstrap",    no_argument,         &do_bootstrap,    'o' },
      { "fracevts",     required_argument,   0,                'c' },
      // If the lmbl flag is not entered, lengthscale_mbl has default value -1.
      // This instructs the code to not use mbl in the fit.
      // The same goes for each other kinematic variable.
      { "mbl",          no_argument,         &do_mbl,          'b' },
      { "mt2",          no_argument,         &do_mt2,          't' },
      { "help",         no_argument,         NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdexahponbtm:c:", longopts, NULL)) != -1 ) {
      switch(c)
      {
         case 'n' :
            run_number = atoi(optarg);
            break;

         case 'f' :
            do_fit = 1;
            break;

         case 'd' :
            do_diagnostics = 1;
            break;

         case 'e' :
            do_templates = 1;
            break;

         case 'x' :
            do_learnparams = 1;
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

         case 'o' :
            do_bootstrap = 1;
            break;

         case 'c' :
            fracevts = atof(optarg);
            break;

         case 'b' :
            do_mbl = true;
            break;

         case 't' :
            do_mt2 = true;
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

   fitter.dists["mbl"].activate = do_mbl;
   fitter.dists["mt2_220_nomatchmbl"].activate = do_mt2;

   // Check that at least one kinematic variable's lengthscale has been entered.
   // Any additional distributions need to be added here
   if (!(fitter.dists["mbl"].activate) and !(fitter.dists["mt2_220_nomatchmbl"].activate) and do_fit == 1){
      std::cout << "At least one variable needed to do fit.  Input at least one lengthscale." << std::endl;
      print_usage();
      return -1;
   }

   fitter.LoadDatasets( datasets );

   // for event counting
   map<string, int> datacount;

   // random number seed for bootstrapping (turns on when nonzero)
   int randseed = 0;
   if( do_bootstrap ) randseed = run_number+1+10E6;

   cout << "\nLoading datasets" << endl;
   for(map<string, Dataset>::iterator it = datasets.begin(); it != datasets.end(); it++){

      string name = it->first;
      Dataset *dat = &(it->second);

      datacount[name] = 0;

      TFile file( (dat->path+dat->file).c_str() );
      TTree *trees = (TTree*)file.Get("RealData");

      cout << setiosflags(ios::left);
      cout << "... " << setw(25) << name
         << ": " << trees->GetEntries() << " events" << endl;

      if( do_diagnostics or use_data ){
         fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
               "RealData", eventvec_datamc, 0, 0, -1 );
      }

      // events for training and testing
      if( name.compare("data") != 0 ){

         if( use_data ){ // train on full mc set
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_train, 0, randseed, -1 );
         }else{
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_train, 1, randseed, -1 );
            fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts,
                  "RealData", eventvec_test, 2, randseed, fracevts );
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

      if( use_data ){ // added 220 distribution to the data fit too, but haven't tested it

         //
         // turn this feature off for now -- will need to clean up later.
         //

         /*
         cout << "REMINDER: EVENT WEIGHTS IN MC" << endl;
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

         // GP training for mbl
         if ( lengthscale_mbl != -1){
            Shapes * fptr = new Shapes( "mbl", fitter.gplength_mbl, fitter.gplength_mt_mbl,
                  fitter.lbnd, fitter.rbnd );
            fptr->TrainGP( hists_train_ );
            fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            fitter.aGPsig = fptr->aGPsig;
            fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            fitter.aGPbkg = fptr->aGPbkg;
         }

         // GP training for 220
         if ( lengthscale_220 != -1){
            Shapes * fptr220 = new Shapes( "mt2_220_nomatchmbl", fitter.gplength_220, fitter.gplength_mt_220,
                  fitter.lbnd, fitter.rbnd );
            fptr220->TrainGP( hists_train_ );
            fitter.aGPsig220.ResizeTo( fptr220->aGPsig.GetNoElements() );
            fitter.aGPsig220 = fptr220->aGPsig;
            fitter.aGPbkg220.ResizeTo( fptr220->aGPbkg.GetNoElements() );
            fitter.aGPbkg220 = fptr220->aGPbkg;
         }

         // events for fitting, hists for training
         fitter.RunMinimizer( eventvec_fit );
         fitter.PlotResults( hists_fit_ ); // plot fitted events

         // fill results tree
         mcmass = 0;
         fitstatus = fitter.gMinuit->Status();
         const double *par = fitter.gMinuit->X();
         const double *par_err = fitter.gMinuit->Errors();
         mt = par[0];
         kmbl = par[1];
         k220 = par[2];
         mt_err = par_err[0];
         kmbl_err = par_err[1];
         k220_err = par_err[2];
         fitchi2 = fitter.fitchi2;
         gplength_mbl = lengthscale_mbl;
         gplength_220 = lengthscale_220;
         gplength_mt_mbl = lengthscale_mt_mbl;
         gplength_mt_220 = lengthscale_mt_220;

         // TODO
         //tree->Fill();

         eventvec_fit.clear();
         eventvec_fit_bkgcontrol.clear();
         fitter.DeleteHists( hists_fit_ );
         fitter.DeleteHists( hists_fit_bkgcontrol_ );
         */

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

            fitter.ReweightMC( eventvec_fit, dname );

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
            for( map<string, Distribution>::iterator it = fitter.dists.begin(); it != fitter.dists.end(); it++ ){

               string name = it->first;
               Distribution *dist = &(it->second);

               double m2llsig, m2llbkg;

               if( dist->activate ){
                  Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2 );
                  fptr->TrainGP( hists_train_, m2llsig, m2llbkg );

                  dist->aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
                  dist->aGPsig = fptr->aGPsig;
                  dist->aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
                  dist->aGPbkg = fptr->aGPbkg;

                  delete fptr;

               }

            }

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

            // events for fitting, hists for training
            fitter.RunMinimizer( eventvec_fit );
            fitter.PlotResults( hists_fit_ ); // plot fitted events

            cout << "Fit Chi2 = " << fitter.fitchi2 << endl;


            // fill results tree
            // any additional variables need to be added here
            mcmass = mass;
            fitstatus = fitter.gMinuit->Status();
            const double *par = fitter.gMinuit->X();
            const double *par_err = fitter.gMinuit->Errors();
            mt = par[0];
            kmbl = par[1];
            mt_err = par_err[0];
            kmbl_err = par_err[1];
            k220 = par[2];
            k220_err = par_err[2];
            fitchi2 = fitter.fitchi2;

            tree->Fill();

            eventvec_fit.clear();
            eventvec_fit_bkgcontrol.clear();
            fitter.DeleteHists( hists_fit_ );
            fitter.DeleteHists( hists_fit_bkgcontrol_ );

         }

      }

   }

   if( do_templates ){

            // do GP training
            for( map<string, Distribution>::iterator it = fitter.dists.begin(); it != fitter.dists.end(); it++ ){

               string name = it->first;
               Distribution *dist = &(it->second);

               double m2llsig, m2llbkg;

               if( dist->activate ){
                  Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2 );
                  fptr->TrainGP( hists_train_, m2llsig, m2llbkg );

                  dist->aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
                  dist->aGPsig = fptr->aGPsig;
                  dist->aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
                  dist->aGPbkg = fptr->aGPbkg;

                  dist->Ainv_sig.ResizeTo( fptr->aGPsig.GetNoElements(), fptr->aGPsig.GetNoElements() );
                  dist->Ainv_sig = fptr->Ainv_sig;
                  dist->Ainv_bkg.ResizeTo( fptr->aGPbkg.GetNoElements(), fptr->aGPbkg.GetNoElements() );
                  dist->Ainv_bkg = fptr->Ainv_bkg;

                  delete fptr;
               }

            }

            fitter.PlotTemplates( hists_train_ );

            for(int j=0; j < 8; j++){
               tsig_mbl_chi2[j] = fitter.tsig_mbl_chi2[j];
               tbkg_mbl_chi2[j] = fitter.tbkg_mbl_chi2[j];
            }

   }

   if( do_learnparams ){
      string name = "mbl";
      Distribution *dist = &(fitter.dists[name]);
      Shapes * fptr = new Shapes( name, dist->glx, dist->glmt, dist->gnorm1, dist->gnorm2 );
      fptr->LearnGPparams( hists_train_ );
      delete fptr;
   }

   if( do_fit or do_templates ){
      tree->Fill();
   }


   //
   // write fit results
   //
   if( do_fit or do_templates ){
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

