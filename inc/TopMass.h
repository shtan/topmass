#ifndef TOPMASS_H
#define TOPMASS_H

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

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
   }

};

struct Hists {

   string name;
   string title;
   int nbins;
   double xlow;
   double xup;

   map<string, TH1D*> proc;
   TH1D *httbar_signal;
   TH1D *httbar_mistag;
   TH1D *httbar_taus;
   TH1D *httbar_hadronic;
   TH1D *hother;

   Hists(){

   }

};

class Fitter{
   public:

      Fitter();
      ~Fitter();

      void ReadNtuple(string, string, double, string, vector<Event>&);
      void LoadDatasets(map<string, Dataset>&);
      void GetVariables(vector<Event>&);

      void RunMinimizer(vector<Event>&);
      double Fmbl_sig(double*, double*);
      TH1D *Hmbl_bkg;

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      // diagnostics
      map< string, map<string, TH1D*> > hists_;
      void DeclareHists();
      void FillHists(vector<Event>&);
      void PrintHists();
      void PlotTemplates();

   private:

      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<Event>* eventvec_fit;

};

#endif

