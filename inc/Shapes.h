#include "TH1D.h"
#include "TVectorD.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

      Shapes( double, double, double, double );
      ~Shapes();

      double Fmbl_tot(double*, double*);
      double Fmbl_sig_param(double, double);
      double Fmbl_sig_gp(double, double);
      double Fmbl_bkg_gp(double, double);

      double lmbl, lmass, gnorm;
      double lbmbl, rbmbl;
      vector<double> ptrain;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      void SetGPopts();
      void TrainGP( map< string, map<string, TH1D*> >& );
      double GPkern(double, double, double, double, double, double, double);

};

#endif

