#include "TopMass.h"

#include <vector>

using namespace std;

int main(int argc, char* argv[]){

   // declarations
   Fitter fitter;
   map<string, Dataset> datasets;
   vector<Event> eventvec;


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
   //fitter.PrintHists();

   fitter.RunMinimizer( eventvec );

   //fitter.PlotTemplates();

   return 0;
}
