#ifndef TOPMASS_H
#define TOPMASS_H

#include "Shapes.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TVectorD.h"

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>

using namespace std;

struct Dataset {

   // dataset info
   string path;
   string file;

   // monte carlo
   double mc_xsec;
   int mc_nevts;

   // constructor
   Dataset( string p="", string f="" ): path(p), file(f) {
      mc_xsec = 1.0;
      mc_nevts = 1;
   }

};

struct Event {

   string process;

   // classify events
   string type;

   int nvertices;
   double weight;

   // kinematic variables
   double mt2_220;
   double mt2_210;
   double mt2_221;
   vector<double> mbls;

   // reconstructed objects
   TLorentzVector jet1, jet2, lep1, lep2, met;

   // for fit
   bool fit_event;

   Event(){
      process = "";
      type = "";

      nvertices = 0;
      weight = 0;

      jet1 = TLorentzVector();
      jet2 = TLorentzVector();
      lep1 = TLorentzVector();
      lep2 = TLorentzVector();
      met = TLorentzVector();

      mt2_220 = 0;
      mt2_210 = 0;
      mt2_221 = 0;

      fit_event = false;
   }

};

class Fitter{

   public:

      Fitter();
      ~Fitter();

      void ReadNtuple(string, string, double, string, vector<Event>&, int=0, int=0);
      void LoadDatasets(map<string, Dataset>&);
      void GetVariables(vector<Event>&);

      void RunMinimizer(vector<Event>&);
      void PlotResults(map< string, map<string, TH1D*> >&);

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      // diagnostics
      void DeclareHists( map< string, map<string, TH1D*> >&, string label );
      void DeleteHists( map< string, map<string, TH1D*> >& );
      void FillHists( map< string, map<string, TH1D*> >&, vector<Event>&, bool=false );
      void PrintHists( map< string, map<string, TH1D*> >& );
      void PlotTemplates( map< string, map<string, TH1D*> >& );
      void PlotTemplatesOld( map< string, map<string, TH1D*> >& );

      TVectorD aGPsig;
      TVectorD aGPbkg;
      TH1D* hmbl_bkg;
      TVectorD aGPsig220;
      TVectorD aGPbkg220;
      TH1D* h221_bkg;

      double rangembl;
      double range220;
      double rbnd, lbnd;
      bool compute_profile;
      double fitchi2;
      double gplength_mbl;
      double gplength_220;
      double gplength_mt;
      
   private:

      double Min2LL(const double*);
      double Min2LL_individual(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<Event>* eventvec_fit;

};

#endif

