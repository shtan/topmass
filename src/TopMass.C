#include "TopMass.h"
#include "Shapes.h"
#include "Mt2Calculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

using namespace std;


//
// constructor and destructor
//

Fitter::Fitter(){

   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // compute weighted errors
   TH1::SetDefaultSumw2();

   // plot formatting
   gROOT->ProcessLineSync(".L scripts/tdrstyle.C");
   gROOT->ProcessLineSync("setTDRStyle()");


}

Fitter::~Fitter(){
      if (gMinuit) delete gMinuit;
      if (fFunc) delete fFunc;
}


//
// member definitions
//

void Fitter::LoadDatasets( map<string, Dataset>& datasets ){

   // file path
   string path = "root://cmseos:1094//eos/uscms/store/user/nmirman/Ntuples/TopMass/20140322/";

   // filenames
   datasets[ "data" ]      = Dataset( path, "ntuple_data.root" );
   datasets[ "t_tw" ]      = Dataset( path, "ntuple_T_tW.root" );
   datasets[ "tbar_tw" ]   = Dataset( path, "ntuple_Tbar_tW.root" );
   datasets[ "dy" ]        = Dataset( path, "ntuple_DYJetsToLL.root" );
   datasets[ "wjets" ]     = Dataset( path, "ntuple_WJetsToLNu.root" );
   datasets[ "ww" ]        = Dataset( path, "ntuple_WW.root" );
   datasets[ "wz" ]        = Dataset( path, "ntuple_WZ.root" );
   datasets[ "zz" ]        = Dataset( path, "ntuple_ZZ.root" );
   datasets[ "ttbar161" ]  = Dataset( path, "ntuple_TTJets_mass161_5.root" ); 
   datasets[ "ttbar163" ]  = Dataset( path, "ntuple_TTJets_mass163_5.root" ); 
   datasets[ "ttbar166" ]  = Dataset( path, "ntuple_TTJets_mass166_5.root" ); 
   datasets[ "ttbar169" ]  = Dataset( path, "ntuple_TTJets_mass169_5.root" ); 
   datasets[ "ttbar172" ]  = Dataset( path, "ntuple_TTJets.root" );
   datasets[ "ttbar175" ]  = Dataset( path, "ntuple_TTJets_mass175_5.root" ); 
   datasets[ "ttbar178" ]  = Dataset( path, "ntuple_TTJets_mass178_5.root" ); 
   datasets[ "ttbar181" ]  = Dataset( path, "ntuple_TTJets_mass181_5.root" ); 

   // for mc weights
   datasets[ "dy" ].mc_nevts         = 30459503;
   datasets[ "ww" ].mc_nevts         = 10000431;
   datasets[ "wz" ].mc_nevts         = 10000283;
   datasets[ "zz" ].mc_nevts         = 9799908;
   datasets[ "tbar_tw" ].mc_nevts    = 493460;
   datasets[ "t_tw" ].mc_nevts       = 497658;
   datasets[ "wjets" ].mc_nevts      = 57709905;
   datasets[ "ttbar161" ].mc_nevts   = 6923750;
   datasets[ "ttbar163" ].mc_nevts   = 6923750;
   datasets[ "ttbar166" ].mc_nevts   = 6923750;
   datasets[ "ttbar169" ].mc_nevts   = 6923750;
   datasets[ "ttbar172" ].mc_nevts   = 6923750;
   datasets[ "ttbar175" ].mc_nevts   = 6923750;
   datasets[ "ttbar178" ].mc_nevts   = 6923750;
   datasets[ "ttbar181" ].mc_nevts   = 6923750;

   datasets[ "dy" ].mc_xsec          = 3351.97;
   datasets[ "ww" ].mc_xsec          = 54.838;
   datasets[ "wz" ].mc_xsec          = 33.21;
   datasets[ "zz" ].mc_xsec          = 8.059;
   datasets[ "tbar_tw" ].mc_xsec     = 11.1;
   datasets[ "t_tw" ].mc_xsec        = 11.1;
   datasets[ "wjets" ].mc_xsec       = 37509.0;
   datasets[ "ttbar161" ].mc_xsec    = 234;
   datasets[ "ttbar163" ].mc_xsec    = 234;
   datasets[ "ttbar166" ].mc_xsec    = 234;
   datasets[ "ttbar169" ].mc_xsec    = 234;
   datasets[ "ttbar172" ].mc_xsec    = 234;
   datasets[ "ttbar175" ].mc_xsec    = 234;
   datasets[ "ttbar178" ].mc_xsec    = 234;
   datasets[ "ttbar181" ].mc_xsec    = 234;

}

void Fitter::ReadNtuple( string path, string process, double mcweight, 
      string selection, vector<Event>& eventvec ){
   cout << "... " << process << endl;
   
   // declare variables
   TLorentzVector *jet1 = new TLorentzVector();
   TLorentzVector *jet2 = new TLorentzVector();
   TLorentzVector *lep1 = new TLorentzVector();
   TLorentzVector *lep2 = new TLorentzVector();
   TLorentzVector *met = new TLorentzVector();

   int nvert;
   double jet1PtRes, jet1PhiRes, jet1EtaRes, jet2PtRes, jet2PhiRes, jet2EtaRes;
   int lpPdgIdGEN, lmPdgIdGEN, nPdgIdGEN, nbPdgIdGEN;
   int jet1GenId, jet2GenId;
   float puMyWeight = 1.0;

   // open ntuple
   TFile file( path.c_str() );
   TTree *tree = (TTree*)file.Get(selection.c_str());

   tree->SetBranchAddress("jet1FourVector", &jet1);
   tree->SetBranchAddress("jet2FourVector", &jet2);
   tree->SetBranchAddress("lepton1FourVector", &lep1);
   tree->SetBranchAddress("lepton2FourVector", &lep2);
   tree->SetBranchAddress("metFourVector", &met);
   tree->SetBranchAddress("jet1PtResolution", &jet1PtRes);
   tree->SetBranchAddress("jet1PhiResolution", &jet1PhiRes);
   tree->SetBranchAddress("jet1EtaResolution", &jet1EtaRes);
   tree->SetBranchAddress("jet2PtResolution", &jet2PtRes);
   tree->SetBranchAddress("jet2PhiResolution", &jet2PhiRes);
   tree->SetBranchAddress("jet2EtaResolution", &jet2EtaRes);
   tree->SetBranchAddress("vertices", &nvert);

   if( process.find("ttbar") != string::npos ){
      tree->SetBranchAddress("lpPdgIdGEN", &lpPdgIdGEN);
      tree->SetBranchAddress("lmPdgIdGEN", &lmPdgIdGEN);
      tree->SetBranchAddress("nPdgIdGEN", &nPdgIdGEN);
      tree->SetBranchAddress("nbPdgIdGEN", &nbPdgIdGEN);
      tree->SetBranchAddress("jet1GenId", &jet1GenId);
      tree->SetBranchAddress("jet2GenId", &jet2GenId);
   }

   if( !(process.find("data") != string::npos) ){
      tree->SetBranchAddress("puMyWeight", &puMyWeight);
   }

   // fill event vector
   for( int ev=0; ev < tree->GetEntries(); ev++ ){

      tree->GetEntry(ev);

      Event evtemp;

      // global quantities
      evtemp.process = process;
      evtemp.weight = mcweight * puMyWeight;
      evtemp.nvertices = nvert;

      // jets, leptons, met
      evtemp.jet1 = *jet1;
      evtemp.jet2 = *jet2;

      evtemp.lep1 = *lep1;
      evtemp.lep2 = *lep2;

      evtemp.met = *met;

      //
      // classify events
      // 
      if( process.find("data") != string::npos ){
         evtemp.type = process;
      }
      // physics backgrounds
      else if( process.find("ttbar") == string::npos ){
         evtemp.type = "other";
      }
      // hadronic decays
      else if( fabs(lpPdgIdGEN) < 5 or fabs(lmPdgIdGEN) < 5
            or fabs(nPdgIdGEN) < 5 or fabs(nbPdgIdGEN) < 5 ){
         evtemp.type = process+"_hadronic";
      }
      // tau decays
      else if( fabs(lpPdgIdGEN) == 15 or fabs(lmPdgIdGEN) == 15
            or fabs(nPdgIdGEN) == 15 or fabs(nbPdgIdGEN) == 15 ){
         evtemp.type = process+"_taus";
      }
      // mistag bkg
      else if( fabs(jet1GenId) != 5 or fabs(jet2GenId) != 5 ){
         evtemp.type = process+"_mistag";
      }
      // signal ttbar
      else if( (lpPdgIdGEN == -11 or lpPdgIdGEN == -13)
            and (lmPdgIdGEN == 11 or lmPdgIdGEN == 13)
            and (nPdgIdGEN == 12 or nPdgIdGEN == 14)
            and (nbPdgIdGEN == -12 or nbPdgIdGEN == -14)
            and (fabs(jet1GenId) == 5 and fabs(jet2GenId) == 5) ){
         evtemp.type = process+"_signal";
      }
      // check for unclassified events
      else{
         cout << "ERROR IN TTBAR EVENT CLASSIFICATION" << endl;
         break;
      }

      // push back event
      eventvec.push_back( evtemp );

   }

   return;
}


void Fitter::GetVariables( vector<Event>& eventvec ){

   Mt2Calculator::Calculator Calc;
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){

      Calc.SetParticles( ev->jet1, ev->jet2, ev->lep1, ev->lep2, ev->met );

      ev->mt2_221 = Calc.GetMt2(2,1);
      ev->mt2_220 = Calc.GetMt2(2,0);
      ev->mt2_210 = Calc.GetMt2(1,0);
      ev->mbls = Calc.GetBlInvariantMasses();

   }

}


void Fitter::RunMinimizer( vector<Event>& eventvec ){

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetTolerance(0.001);
   gMinuit->SetPrintLevel(3);

   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 2 );
   gMinuit->SetFunction( *fFunc );
   gMinuit->SetVariable(0, "topMass", 173.0, 0.01);
   gMinuit->SetLimitedVariable(1, "norm", 0.5, 0.01, 0, 1.0);

   // set event vector and minimize
   eventvec_fit = &eventvec;
   gMinuit->Minimize();
   gMinuit->Hesse();

   const double *xmin = gMinuit->X();

   // plot fit results
   TFile *fileout = new TFile( "results/plotsFitResults.root", "RECREATE" );
   fileout->cd();
   
   TCanvas *cfit_mbl = new TCanvas("cfit_mbl","cfit_mbl",800,800);
   cfit_mbl->cd();

   hists_["mbl_r1000"]["data"]->Draw();

   TH1D *hbkg = (TH1D*)hists_["mbl_r1000"]["ttbar172_mistag"]->Clone("hbkg2");
   hbkg->Add( hists_["mbl_r1000"]["ttbar172_taus"] );
   hbkg->Add( hists_["mbl_r1000"]["ttbar172_hadronic"] );
   hbkg->Add( hists_["mbl_r1000"]["other"] );
   Shapes * fptr = new Shapes( hists_["mbl_r1000"]["data_bkgcontrol"] );
   //Shapes * fptr = new Shapes( hists_["mbl_r1000"]["ttbar172_mistag"] );
   //Shapes * fptr = new Shapes( hbkg );
   TF1 *fmbl_tot = new TF1("fmbl_tot", fptr, &Shapes::Fmbl_tot, 0, 1000, 4);
   // evaluate integral
   fmbl_tot->SetParameters( xmin[0], 1.0, 1.0, 1.0 );
   double integral = fmbl_tot->Integral(0,1000);
   fmbl_tot->SetParameters( xmin[0], xmin[1], hists_["mbl_r1000"]["data"]->Integral("width"), integral );
   fmbl_tot->DrawCopy("same");

   cfit_mbl->Write();

   fileout->Close();

}

double Fitter::Min2LL(const double *x){

   // temp -- evaluate integral
   Shapes * fptr = new Shapes( hists_["mbl_r1000"]["data_bkgcontrol"] );
   TF1 *fmbl_tot = new TF1("fmbl_tot", fptr, &Shapes::Fmbl_tot, 0, 1000, 4);
   fmbl_tot->SetParameters( x[0], 1.0, 1.0, 1.0 );
   double integral = fmbl_tot->Integral(0,1000);
   delete fmbl_tot;
   delete fptr;

   TH1D *hbkg = (TH1D*)hists_["mbl_r1000"]["ttbar172_mistag"]->Clone("hbkg");
   hbkg->Add( hists_["mbl_r1000"]["ttbar172_taus"] );
   hbkg->Add( hists_["mbl_r1000"]["ttbar172_hadronic"] );
   hbkg->Add( hists_["mbl_r1000"]["other"] );
   Shapes shape( hists_["mbl_r1000"]["data_bkgcontrol"] );
   //Shapes shape( hists_["mbl"]["ttbar172_mistag"] );
   //Shapes shape( hbkg );
   // evaluate likelihood
   double pmbl [] = {x[0], x[1], 1.0, integral};
   double m2ll = 0;
   for( vector<Event>::iterator ev = eventvec_fit->begin(); ev < eventvec_fit->end(); ev++){
      if( ev->process.compare("data") != 0 ) continue;
      for( unsigned int i=0; i < ev->mbls.size(); i++ ){

         double Lmbl = shape.Fmbl_tot( &(ev->mbls[i]), pmbl );
         m2ll -= 2.0*ev->weight*log( Lmbl );

      }
   }

   return m2ll;
}

