#include "TH1D.h"
#include "TVectorD.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

      //Shapes( map< string, map<string,TH1D*> >& );
      Shapes( TH1D*& );
      ~Shapes();

      double Fmbl_tot(double*, double*);
      double Fmbl_sig_param(double, double);
      double Fmbl_sig_gp(double, double);
      double Fmbl_bkg_gp(double, double);
      double Fmbl_bkg(double);


      double lmbl, lmass, gnorm;
      vector<double> ptrain;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      void SetGPopts();
      void TrainGP( map< string, map<string, TH1D*> >& );
      double GPkern(double, double, double, double, double, double, double);

      //map< string, map<string, TH1D*> > hists_; // copy from TopMass class

      TH1D* hmbl_bkg;
      double norm_mbl_bkg;

};

#endif

