#include "TH1D.h"
#include "TVectorD.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

      Shapes( string, double, double, double, double );
      ~Shapes();

      double Ftot(double*, double*);
      double Fsig_param(double, double);
      double Fsig_gp(double, double);
      double Fbkg_gp(double, double);

      string name;
      double lx, lmass, gnorm;
      double lbx, rbx;
      vector<double> ptrain;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      void SetGPopts();
      void TrainGP( map< string, map<string, TH1D*> >& );
      double GPkern(double, double, double, double, double, double, double);

};

#endif

