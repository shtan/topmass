#include "TH1D.h"
#include "TVectorD.h"

#include <map>

#ifndef SHAPES_H
#define SHAPES_H

using namespace std;

class Shapes{

   public:

<<<<<<< HEAD
      Shapes( map< string, map<string,TH1D*> >& );
=======
      //Shapes( map< string, map<string,TH1D*> >& );
      Shapes( TH1D*& );
>>>>>>> 06531be
      ~Shapes();

      double Fmbl_tot(double*, double*);
      double Fmbl_sig_param(double, double);
      double Fmbl_sig_gp(double, double);
<<<<<<< HEAD
      double Fmbl_bkg(double);


      double lmbl, lmass;
      vector<double> ptrain;
      TVectorD aGP;
      void SetGPopts();
      void TrainGP();
      double GPkern(double, double, double, double, double, double);

      map< string, map<string, TH1D*> > hists_; // copy from TopMass class
=======
      double Fmbl_bkg_gp(double, double);
      double Fmbl_bkg(double);


      double lmbl, lmass, gnorm, gnoise;
      vector<double> ptrain;
      TVectorD aGPsig;
      TVectorD aGPbkg;
      void SetGPopts();
      void TrainGP( map< string, map<string, TH1D*> >& );
      double GPkern(double, double, double, double, double, double, double, double);

      //map< string, map<string, TH1D*> > hists_; // copy from TopMass class
>>>>>>> 06531be

      TH1D* hmbl_bkg;
      double norm_mbl_bkg;

};

#endif

