#include "TopMass.h"
#include "Shapes.h"
#include "Mt2Calculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
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

   return;
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

   return;
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

   return;
}

double Fitter::Min2LL(const double *x){

   // normalization inside likelihood function (temp)
   Shapes * fptr = new Shapes( hists_["mbl_fit"]["data_bkgcontrol"] );
   TF1 *fmbl_tot = new TF1("fmbl_tot", fptr, &Shapes::Fmbl_tot, 0, 1000, 4);
   fmbl_tot->SetParameters( x[0], 1.0, 1.0, 1.0 );
   double integral = fmbl_tot->Integral(0,1000);
   delete fmbl_tot;
   delete fptr;

   Shapes shape( hists_["mbl_fit"]["data_bkgcontrol"] );

   double pmbl [] = {x[0], x[1], 1.0, integral};
   double m2ll = 0;
   // evaluate likelihood
   for( vector<Event>::iterator ev = eventvec_fit->begin(); ev < eventvec_fit->end(); ev++){
      if( ev->process.compare("data") != 0 ) continue;

      for( unsigned int i=0; i < ev->mbls.size(); i++ ){
         m2ll -= 2.0*ev->weight*log( shape.Fmbl_tot( &(ev->mbls[i]), pmbl ) );
      }

   }

   return m2ll;
}

void Fitter::PlotResults(){

   const double *xmin = gMinuit->X();
   const double *xerr = gMinuit->Errors();
   double minvalue = gMinuit->MinValue();

   // plot fit results
   TFile *fileout = new TFile( "results/plotsFitResults.root", "RECREATE" );
   fileout->cd();

   string names [] = {"mbl","mbl_fit"};
   for(unsigned int i=0; i < sizeof(names)/sizeof(names[0]); i++){

      TCanvas *canvas = new TCanvas( ("c_"+names[i]).c_str(), ("c_"+names[i]).c_str(), 800, 800);
      canvas->SetFillColor(0);
      canvas->cd();

      TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
      pad1->SetTopMargin(0.1);
      pad1->SetBottomMargin(0.01);
      pad1->SetRightMargin(0.1);
      pad1->SetFillColor(0);
      pad2->SetTopMargin(0.01);
      pad2->SetBottomMargin(0.3);
      pad2->SetRightMargin(0.1);
      pad2->SetFillColor(0);
      pad1->Draw();
      pad2->Draw();

      // line for ratio plot
      TF1 *func = new TF1("func","[0]",-10E6,10E6);
      func->SetParameter(0,1.0);
      func->SetLineWidth(1);
      func->SetLineStyle(7);
      func->SetLineColor(1);

      // pad 1
      pad1->cd();

      TH1D *hdata = (TH1D*)hists_[names[i]]["data"]->Clone("hdata");
      if( names[i].compare("mbl_fit") == 0 ){
         hdata->Rebin(4);
         hdata->GetYaxis()->SetTitle("Events/10 GeV");
      }

      hdata->GetXaxis()->SetTitleSize(0.00);
      hdata->GetYaxis()->SetLabelSize(0.07);
      hdata->GetYaxis()->SetTitleSize(0.08);
      hdata->GetYaxis()->SetTitleOffset(1.0);
      hdata->GetXaxis()->SetLabelFont(42);
      hdata->GetYaxis()->SetLabelFont(42);
      hdata->GetXaxis()->SetTitleFont(42);
      hdata->GetYaxis()->SetTitleFont(42);

      hdata->SetMarkerStyle(20);
      hdata->Draw();

      Shapes * fptr = new Shapes( hists_["mbl_fit"]["data_bkgcontrol"] );
      TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Fmbl_tot, 0, 1000, 4);

      // normalization inside likelihood function (temp)
      ftemplate->SetParameters( xmin[0], 1.0, 1.0, 1.0 );
      double integral = ftemplate->Integral(0,1000);
      ftemplate->SetParameters( xmin[0], xmin[1],
            hdata->Integral("width"), integral );

      ftemplate->SetLineWidth(2);
      ftemplate->DrawCopy("same");
      hdata->DrawCopy("same"); // redraw points

      // pad 2
      pad2->cd();
      TH1D *hratio = (TH1D*)hists_[names[i]]["data"]->Clone("hratio");
      if( names[i].compare("mbl_fit") == 0 ) hratio->Rebin(4);
      hratio->Divide( ftemplate );

      hratio->SetTitle("");
      hratio->GetYaxis()->SetTitle("data/mc");
      hratio->GetYaxis()->CenterTitle();
      hratio->SetStats(0);

      hratio->GetXaxis()->SetTitleSize(0.14);
      hratio->GetXaxis()->SetLabelSize(0.14);
      hratio->GetYaxis()->SetLabelSize(0.11);
      hratio->GetYaxis()->SetTitleSize(0.14);
      hratio->GetYaxis()->SetTitleOffset(0.28);
      hratio->GetXaxis()->SetLabelFont(42);
      hratio->GetYaxis()->SetLabelFont(42);
      hratio->GetXaxis()->SetTitleFont(42);
      hratio->GetYaxis()->SetTitleFont(42);
      hratio->SetMaximum( 1.6 );
      hratio->SetMinimum( 0.4 );
      hratio->GetYaxis()->SetNdivisions(505);

      hratio->Draw("EP");
      func->Draw("same");

      canvas->Write();

      delete canvas;
      delete func;
      delete fptr;
      delete ftemplate;
   }

   //
   // plot likelihood near minimum
   //
   unsigned int npnts_mt = 30;
   unsigned int npnts_kmbl = 30;
   double mt_lrange = xmin[0]-3*xerr[0];
   double mt_rrange = xmin[0]+3*xerr[0];
   double kmbl_lrange = xmin[1]-3*xerr[1];
   double kmbl_rrange = xmin[1]+3*xerr[1];
   
   // mt profile
   TGraph *gLmt = new TGraph();
   for(unsigned int i=0; i <= npnts_mt; i++){
      cout << "mt profile, pnt " << i << endl;
      double mt = mt_lrange + (mt_rrange-mt_lrange)*i/npnts_mt;
      const double par [] = {mt, xmin[1], xmin[2], xmin[3]};
      gLmt->SetPoint(i, mt, Min2LL(par) - minvalue);
   }
   // kmbl profile
   TGraph *gLkmbl = new TGraph();
   for(unsigned int i=0; i <= npnts_kmbl; i++){
      cout << "kmbl profile, pnt " << i << endl;
      double kmbl = kmbl_lrange + (kmbl_rrange-kmbl_lrange)*i/npnts_kmbl;
      const double par [] = {xmin[0], kmbl, xmin[2], xmin[3]};
      gLkmbl->SetPoint(i, kmbl, Min2LL(par) - minvalue);
   }
   // kmbl vs mt
   TH2D *hLmbl = new TH2D("hLmbl", "hLmbl", npnts_mt, mt_lrange, mt_rrange,
         npnts_kmbl, kmbl_lrange, kmbl_rrange);
   for(unsigned int i=0; i <= npnts_mt; i++){
      for(unsigned int j=0; j <= npnts_kmbl; j++){
         cout << "2d profile, pnt " << i << ", " << j << endl;
         double mt = hLmbl->GetXaxis()->GetBinCenter(i);
         double kmbl = hLmbl->GetYaxis()->GetBinCenter(j);
         const double par [] = {mt, kmbl, xmin[2], xmin[3]};
         hLmbl->SetBinContent(i, j, Min2LL(par) - minvalue);
      }
   }

   TCanvas *cLmt = new TCanvas("cLmt", "cLmt", 800, 800);
   cLmt->cd();
   gLmt->SetMarkerStyle(20);
   gLmt->Draw("ACP");
   cLmt->Write("cLmt");

   TCanvas *cLkmbl = new TCanvas("cLkmbl", "cLkmbl", 800, 800);
   cLkmbl->cd();
   gLkmbl->SetMarkerStyle(20);
   gLkmbl->Draw("ACP");
   cLkmbl->Write("cLkmbl");

   TCanvas *cLmbl = new TCanvas("cLmbl", "cLmbl", 800, 800);
   cLmbl->cd();
   hLmbl->Draw("colz");
   cLmbl->Write("cLmbl");

   /*
   // 2 sigma contour
   unsigned int npnts_contour = 100;
   gMinuit->SetErrorDef(4); // note 4 and not 2
   double *x2s = new double[npnts_contour];
   double *y2s = new double[npnts_contour];
   gMinuit->Contour(0, 1, npnts_contour, x2s, y2s); // par1, par2, numpnts, x, y
   TGraph *gs2 = new TGraph(npnts_contour, x2s, y2s);
   gs2->SetLineColor(5);
   gs2->Draw("same C");

   // 1 sigma contour
   gMinuit->SetErrorDef(1); // note 4 and not 2
   double *xs = new double[100];
   double *ys = new double[100];
   gMinuit->Contour(0, 1, npnts_contour, xs, ys); // par1, par2, numpnts, x, y
   TGraph *gs = new TGraph(npnts_contour, xs, ys);
   gs->SetLineColor(5);
   gs->Draw("same C");
   */

   fileout->Write();
   fileout->Close();

   return;
}
