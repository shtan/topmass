#include "TH1D.h"

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

      Shapes(TH1D*&);
      ~Shapes();

      double Fmbl_tot(double*, double*);
      double Fmbl_sig(double, double);
      double Fmbl_bkg(double);

      TH1D* hmbl_bkg;
      double norm_mbl_bkg;

};

#endif

