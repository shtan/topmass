#include "TH1D.h"
#include "TVectorD.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

      // TODO
      Shapes( double, double, double, double, double );
      ~Shapes();

      double Fmbl_tot(double*, double*);
      double Fmbl_sig_param(double, double);
      double Fmbl_gp(double, double, string);
      double Fmbl_gp_var(double, double, string);

      double lmbl, lmass, gnorm;
      double lbmbl, rbmbl;
      vector<double> ptrain;
      TMatrixD Asig;
      TMatrixD Abkg;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      void SetGPopts();
      void TrainGP( map< string, map<string, TH1D*> >& );
      double GPkern(double, double, double, double, double, double);

};

#endif

