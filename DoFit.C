#include "TopMass.h"
#include "TTree.h"
#include "TFile.h"

#include <vector>

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

   while( (c = getopt(argc, argv, "fd")) != -1 ) {
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

   if( do_fit ){
      fitter.RunMinimizer( eventvec );
      fitter.PlotResults(fitter.gMinuit->X());
   }

   fitter.PlotTemplates();


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
