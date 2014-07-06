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
   cout << setw(25) << "\t-a --data" << "Run the fit on data (use full mc for training).\n";
   cout << setw(25) << "\t-p --profile" << "Run the likelihood profile.\n";
   cout << setw(25) << "\t-m --masspnt  <value>" << "If running on mc, use masspoint indicated.\n";
   cout << setw(25) << "\t-b --lmbl" << "Set the mbl lengthscale.\n";
   cout << setw(25) << "\t-t --lmt" << "Set the mt lengthscale.\n";
   cout << setw(25) << "\t-g --gnorm" << "Set noise term normalization parameter.\n";
   cout << setw(25) << "\t-l --lbnd" << "Left bound for exclusion.\n";
   cout << setw(25) << "\t-r --rbnd" << "Right bound for exclusion.\n";
   cout << setw(25) << "\t-o --bootstrap" << "Turn on bootstrapping.\n";
   cout << setw(25) << "\t-c --fracevts" << "Fit fraction of events.\n";
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
   double mcmass=0;
   double fitchi2=0;
   double gplength_mbl=0;
   double gplength_mt=0;
   double gnorm=0;
   double lbound_mbl=0;
   double rbound_mbl=0;
   double tsig_mbl_chi2 [8] = {0};
   double tbkg_mbl_chi2 [8] = {0};
   
   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("runNumber", &run_number);
   tree->Branch("fitStatus", &fitstatus);
   tree->Branch("mt", &mt);
   tree->Branch("mt_err", &mt_err);
   tree->Branch("kmbl", &kmbl);
   tree->Branch("kbml_err", &kmbl_err);
   tree->Branch("mcmass", &mcmass);
   tree->Branch("fitchi2", &fitchi2);
   tree->Branch("gplength_mbl", &gplength_mbl);
   tree->Branch("gplength_mt", &gplength_mt);
   tree->Branch("gnorm", &gnorm);
   tree->Branch("lbound_mbl", &lbound_mbl);
   tree->Branch("rbound_mbl", &rbound_mbl);
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
   double fracevts = -1;

   struct option longopts[] = {
      { "run_number",   required_argument,   0,                'n' },
      { "fit",          no_argument,         &do_fit,          'f' },
      { "diagnostics",  no_argument,         &do_diagnostics,  'd' },
      { "templates",    no_argument,         &do_templates,    'e' },
      { "data",         no_argument,         &use_data,        'a' },
      { "profile",      no_argument,         0,                'p' },
      { "masspnt",      required_argument,   0,                'm' },
      { "lmbl",         required_argument,   0,                'b' },
      { "lmt",          required_argument,   0,                't' },
      { "gnorm",       required_argument,   0,                'g' },
      { "lbnd",         required_argument,   0,                'l' },
      { "rbnd",         required_argument,   0,                'r' },
      { "bootstrap",    no_argument,         &do_bootstrap,    'o' },
      { "fracevts",     required_argument,   0,                'c' },
      { "help",         no_argument,         NULL,             'h' },
      { 0, 0, 0, 0 }
   };

   while( (c = getopt_long(argc, argv, "fdeahpon:m:b:t:g:l:r:c:", longopts, NULL)) != -1 ) {
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

         case 'a' :
            use_data = 1;
            break;

         case 'm' :
            masspnt = atof(optarg);
            break;

         case 'b' :
            fitter.gplength_mbl = atof(optarg);
            break;

         case 't' :
            fitter.gplength_mt = atof(optarg);
            break;

         case 'l' :
            fitter.lbnd = atof(optarg);
            break;

         case 'r' :
            fitter.rbnd = atof(optarg);
            break;

         case 'g' :
            fitter.gnorm2 = atof(optarg);
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

      if( use_data ){
         // TODO
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
         cout << "Training GP... ";
         // TODO
         double m2llsig, m2llbkg;
         Shapes * fptr = new Shapes( fitter.gplength_mbl, fitter.gplength_mt,
               fitter.lbnd, fitter.rbnd, fitter.gnorm1, fitter.gnorm2 );
         fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
         cout << " done." << endl;
         fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
         fitter.aGPsig = fptr->aGPsig;
         fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
         fitter.aGPbkg = fptr->aGPbkg;

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
         mt_err = par_err[0];
         kmbl_err = par_err[1];
         fitchi2 = fitter.fitchi2;

         // TODO
         //tree->Fill();

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

            // TODO
            fitter.ReweightMC( eventvec_fit, dname );
            /* 
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
            */

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

            // TODO
            
            fitter.gplength_mbl = 13;
            //fitter.gplength_mt = 18;
            fitter.gplength_mt = 32;
            //fitter.gnorm1 = 1.5;
            //fitter.gnorm2 = 12;
            //fitter.gnorm1 = 30;
            //fitter.gnorm2 = 1;
            fitter.gnorm1 = 1E6;
            fitter.gnorm2 = 3E4;
            
            /*
            fitter.gplength_mbl = 3.98;
            fitter.gplength_mt = 12.49;
            fitter.gnorm1 = 15.8;
            fitter.gnorm2 = 11.8;
            */


            // do GP training
            // TODO
            double m2llsig, m2llbkg;
            Shapes * fptr = new Shapes( fitter.gplength_mbl, fitter.gplength_mt,
                  fitter.lbnd, fitter.rbnd, fitter.gnorm1, fitter.gnorm2 );
            fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
            fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
            fitter.aGPsig = fptr->aGPsig;
            fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
            fitter.aGPbkg = fptr->aGPbkg;

            // TODO
            //
            /*
            cout << "Testing integral..." << endl;
            //double masspnts [] = {161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5};
            for( int m=160; m < 182; m++ ){
               cout << "mass = " << m << ": " << endl;

               TF1 *fmbl_tot = new TF1("fmbl_tot", fptr, &Shapes::Fmbl_tot, 0, 300, 5);
               fmbl_tot->SetParameters( m, 1.0, 1.0, 1.0, 1.0 );
               double integralsig = fmbl_tot->Integral(0,300);
               fmbl_tot->SetParameters( m, 0.0, 1.0, 1.0, 1.0 );
               double integralbkg = fmbl_tot->Integral(0,300);

               // signal integral
               double pmbl [] = {m, 1.0, 1.0, integralsig, integralbkg};
               double integral = 0;
               for(double x=0; x < 300; x+=0.1){
                  integral += 0.1*fptr->Fmbl_tot( &x, pmbl );
               }
               cout << " -----> sig integral = " << integral << endl;

               // bkg integral
               pmbl[1] = 0.0;
               integral = 0;
               for(double x=0; x < 300; x+=0.1){
                  integral += 0.1*fptr->Fmbl_tot( &x, pmbl );
               }
               cout << " -----> bkg integral = " << integral << endl;

               // s+b integral
               pmbl[1] = 0.5;
               integral = 0;
               for(double x=0; x < 300; x+=0.1){
                  integral += 0.1*fptr->Fmbl_tot( &x, pmbl );
               }
               cout << " -----> s+b integral = " << integral << endl;

            }
            return 0;
            */

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
            mcmass = mass;
            fitstatus = fitter.gMinuit->Status();
            const double *par = fitter.gMinuit->X();
            const double *par_err = fitter.gMinuit->Errors();
            mt = par[0];
            kmbl = par[1];
            mt_err = par_err[0];
            kmbl_err = par_err[1];
            fitchi2 = fitter.fitchi2;
            gplength_mbl = fitter.gplength_mbl;
            gplength_mt = fitter.gplength_mt;
            gnorm = fitter.gnorm2;
            lbound_mbl = fitter.lbnd;
            rbound_mbl = fitter.rbnd;

            // TODO
            //tree->Fill();

            eventvec_fit.clear();
            eventvec_fit_bkgcontrol.clear();
            fitter.DeleteHists( hists_fit_ );
            fitter.DeleteHists( hists_fit_bkgcontrol_ );

         }

      }

   }

   if( do_templates ){
      
      
      /*
      Shapes * fptr2 = new Shapes( fitter.gplength_mbl, fitter.gplength_mt,
            fitter.lbnd, fitter.rbnd, fitter.gnorm1, fitter.gnorm2 );
      fptr2->LearnGPparams( hists_train_ );
      return 0;
*/

      // do GP training
      cout << "Training GP... ";
      double m2llsig, m2llbkg;
      // TODO
      fitter.gplength_mbl = 13;
      //fitter.gplength_mt = 18;
      fitter.gplength_mt = 32;
      //fitter.gnorm1 = 1.5;
      //fitter.gnorm2 = 12;
      fitter.gnorm1 = 1;
      fitter.gnorm2 = 1;
      Shapes * fptr = new Shapes( fitter.gplength_mbl, fitter.gplength_mt,
            fitter.lbnd, fitter.rbnd, fitter.gnorm1, fitter.gnorm2 );
      /*
      for( double g=1.9E-7; g <= 3.0E-7; g+=1E-8 ){
         fptr->gnorm1 = g;
         fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
         cout << "gnorm = " << g << " ---> " << "M2LL (sig, bkg): "
            << m2llsig << ", " << m2llbkg << endl;
      }
         fptr->gnorm1 = 2.8E-07;
         fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
         */
      
      /*
      for( double g=0.1; g <= 2; g+=0.1 ){
         fptr->gnorm1 = 2.8E-07;
         fptr->gnorm2 = g;
         fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
         cout << "gnorm = " << g << " ---> " << "M2LL (sig, bkg): "
            << m2llsig << ", " << m2llbkg << endl;
      }
      return 0;
      */
         fptr->TrainGP( hists_train_, m2llsig, m2llbkg );
      fitter.aGPsig.ResizeTo( fptr->aGPsig.GetNoElements() );
      fitter.aGPsig = fptr->aGPsig;
      fitter.aGPbkg.ResizeTo( fptr->aGPbkg.GetNoElements() );
      fitter.aGPbkg = fptr->aGPbkg;
      // TODO
      fitter.Ainv_sig.ResizeTo( fptr->aGPsig.GetNoElements(), fptr->aGPsig.GetNoElements() );
      fitter.Ainv_sig = fptr->Ainv_sig;
      fitter.Ainv_bkg.ResizeTo( fptr->aGPbkg.GetNoElements(), fptr->aGPbkg.GetNoElements() );
      fitter.Ainv_bkg = fptr->Ainv_bkg;

      fitter.PlotTemplates( hists_train_ );
      for(int j=0; j < 8; j++){
         tsig_mbl_chi2[j] = fitter.tsig_mbl_chi2[j];
         tbkg_mbl_chi2[j] = fitter.tbkg_mbl_chi2[j];
      }

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

