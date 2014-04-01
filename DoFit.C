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

      fitter.ReadNtuple( dat->path+dat->file, name, dat->mc_xsec/dat->mc_nevts, eventvec );

   }

   cout << endl;
   return 0;
}
