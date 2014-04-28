#include "TopMass.h"
#include "Shapes.h"
#include "Mt2Calculator.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
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
#include <sstream>

using namespace std;


void Fitter::DeclareHists( map< string, map<string, TH1D*> >& hists_, string label ){

   vector<string> type;
   type.push_back("data");
   type.push_back("data_bkgcontrol");
   type.push_back("other");
   type.push_back("other_bkgcontrol");
   string masspnts [] = {"161","163","166","169","172","175","178","181"};
   for(int i=0; i < 8; i++){
      type.push_back("ttbar"+masspnts[i]+"_signal");
      type.push_back("ttbar"+masspnts[i]+"_mistag");
      type.push_back("ttbar"+masspnts[i]+"_taus");
      type.push_back("ttbar"+masspnts[i]+"_hadronic");
      type.push_back("ttbar"+masspnts[i]+"_bkgcontrol_signal");
      type.push_back("ttbar"+masspnts[i]+"_bkgcontrol_mistag");
      type.push_back("ttbar"+masspnts[i]+"_bkgcontrol_taus");
      type.push_back("ttbar"+masspnts[i]+"_bkgcontrol_hadronic");
   }
   type.push_back("fitevts"); // for displaying fit results (could be mc or data)

   for(vector<string>::iterator t = type.begin(); t < type.end(); t++){
         
      string name = *t;
      string namel = *t+label;

      //
      // variables
      //

      hists_["mt2_220"][name] = new TH1D( ("hmt2_220_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["mt2_221"][name] = new TH1D( ("hmt2_221_"+namel).c_str(),
            "M_{T2} 221;M_{T2} 221 (GeV);Events/2 GeV", 100, 50, 250 );
      hists_["mbl"][name] = new TH1D( ("hmbl_"+namel).c_str(),
            "M_{bl};M_{bl} (GeV);Events/3 GeV", 100, 0, 300 );
      hists_["mbl_fit"][name] = new TH1D( ("hmbl_fit_"+namel).c_str(),
            "M_{bl};M_{bl} (GeV);Events/2.5 GeV", 400, 0, 1000 );
      hists_["mt2_210"][name] = new TH1D( ("hmt2_210_"+namel).c_str(),
            "M_{T2} 210;M_{T2} 210 (GeV);Events/1.5 GeV", 100, 0, 150 );
      hists_["mt2_220_matchmbl"][name] = new TH1D( ("hmt2_220_matchmbl_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["mt2_220_nomatchmbl"][name] = new TH1D( ("hmt2_220_nomatchmbl_"+namel).c_str(),
            "M_{T2} 220;M_{T2} 220 (GeV);Events/2.5 GeV", 100, 0, 250 );

      //hists2d_["mblV220"][name] = new TH2D( ("hmblV220_"+namel).c_str(),
      //      "220 vs. mbl", 50, 0, 250, 50, 0, 250 );

      //
      // kinematic distributions
      //

      // pt, eta, phi
      hists_["b_pt"][name] = new TH1D( ("hb_pt_"+namel).c_str(),
            "b quark p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["b_eta"][name] = new TH1D( ("hb_eta_"+namel).c_str(),
            "b quark #eta;#eta;Events", 100, -3, 3 );
      hists_["b_phi"][name] = new TH1D( ("hb_phi_"+namel).c_str(),
            "b quark #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["l_pt"][name] = new TH1D( ("hl_pt_"+namel).c_str(),
            "lepton p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["l_eta"][name] = new TH1D( ("hl_eta_"+namel).c_str(),
            "lepton #eta;#eta;Events", 100, -3, 3 );
      hists_["l_phi"][name] = new TH1D( ("hl_phi_"+namel).c_str(),
            "lepton #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["bb_pt"][name] = new TH1D( ("hbb_pt_"+namel).c_str(),
            "bb p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["bb_eta"][name] = new TH1D( ("hbb_eta_"+namel).c_str(),
            "bb Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["bb_phi"][name] = new TH1D( ("hbb_phi_"+namel).c_str(),
            "bb #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["ll_pt"][name] = new TH1D( ("hll_pt_"+namel).c_str(),
            "ll p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["ll_eta"][name] = new TH1D( ("hll_eta_"+namel).c_str(),
            "ll Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["ll_phi"][name] = new TH1D( ("hll_phi_"+namel).c_str(),
            "ll #phi;#phi;Events", 100, -3.5, 3.5 );
      hists_["bbll_pt"][name] = new TH1D( ("hbbll_pt_"+namel).c_str(),
            "bbll p_{T};p_{T} (GeV);Events/2.5 GeV", 100, 0, 250 );
      hists_["bbll_eta"][name] = new TH1D( ("hbbll_eta_"+namel).c_str(),
            "bbll Rapidity;Rapidity;Events", 100, -3, 3 );
      hists_["bbll_phi"][name] = new TH1D( ("hbbll_phi_"+namel).c_str(),
            "bbll #phi;#phi;Events", 100, -3.5, 3.5 );

      // invariant mass
      hists_["bb_m"][name] = new TH1D( ("hbb_m_"+namel).c_str(),
            "bb invariant mass;m_{bb} (GeV);Events/4 GeV", 100, 0, 400 );
      hists_["ll_m"][name] = new TH1D( ("hll_m_"+namel).c_str(),
            "ll invariant mass;m_{ll} (GeV);Events/4 GeV", 100, 0, 400 );
      hists_["bbl_m"][name] = new TH1D( ("hbbl_m_"+namel).c_str(),
            "bbl invariant mass;m_{bbl} (GeV);Events/6 GeV", 100, 0, 600 );
      hists_["bll_m"][name] = new TH1D( ("hbll_m_"+namel).c_str(),
            "bll invariant mass;m_{bll} (GeV);Events/6 GeV", 100, 0, 600 );
      hists_["bbll_m"][name] = new TH1D( ("hbbll_m_"+namel).c_str(),
            "bbll invariant mass;m_{bbll} (GeV);Events/10 GeV", 100, 0, 1000 );

      // angles
      hists_["bb_dR"][name] = new TH1D( ("hbb_dR_"+namel).c_str(),
            "#DeltaR between bb;#DeltaR;Events", 100, 0, 5.0 );
      hists_["bl_dR"][name] = new TH1D( ("hbl_dR_"+namel).c_str(),
            "#DeltaR between bl;#DeltaR;Events", 100, 0, 5.0 );
      hists_["ll_dR"][name] = new TH1D( ("hll_dR_"+namel).c_str(),
            "#DeltaR between ll;#DeltaR;Events", 100, 0, 5.0 );

      // met
      hists_["met_x"][name] = new TH1D( ("hmet_x_"+namel).c_str(),
            "MET_{x};MET_{x} (GeV);Events/4 GeV", 100, -200, 200 );
      hists_["met_y"][name] = new TH1D( ("hmet_y_"+namel).c_str(),
            "MET_{y};MET_{y} (GeV);Events/4 GeV", 100, -200, 200 );
      hists_["met"][name] = new TH1D( ("hmet_"+namel).c_str(),
            "total MET;MET (GeV);Events/3 GeV", 100, 0, 300 );
   }

   return;
}

void Fitter::DeleteHists( map< string, map<string, TH1D*> >& hists_ ){

   typedef map<string, TH1D*> tmap;
   typedef map<string, tmap> hmap;

   for( hmap::iterator h = hists_.begin(); h != hists_.end(); h++){
      for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
         //t->second->Delete();
         delete (t->second);
      }
   }
   

}

void Fitter::FillHists( map< string, map<string, TH1D*> >& hists_,
      vector<Event>& eventvec, bool fit_events ){

   // event loop
   for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){

      string type = fit_events ? "fitevts" : ev->type;

      if( fit_events and !(ev->fit_event) ){
         cout << "ERROR: FILL WITH FIT EVENTS" << endl;
         return;
      }

      // define objects
      TLorentzVector jet1 = ev->jet1;
      TLorentzVector jet2 = ev->jet2;
      TLorentzVector lep1 = ev->lep1;
      TLorentzVector lep2 = ev->lep2;
      TLorentzVector met = ev->met;
      TLorentzVector up221 = -((jet1)+(jet2)+(lep1)+(lep2)+(met));
      TLorentzVector bb = jet1 + jet2;
      TLorentzVector ll = lep1 + lep2;
      TLorentzVector bbll = jet1 + jet2 + lep1 + lep2;
      TLorentzVector bll1 = jet1 + lep1 + lep2;
      TLorentzVector bll2 = jet2 + lep1 + lep2;
      TLorentzVector bbl1 = jet1 + jet2 + lep1;
      TLorentzVector bbl2 = jet1 + jet2 + lep2;

      if ( !(jet1.M() < 40 and jet2.M() < 40) ) continue;

      if (sin((jet1).DeltaPhi(up221))*sin((jet2).DeltaPhi(up221)) > 0){
         hists_["mt2_221"][type]->Fill( ev->mt2_221, ev->weight );
      }

      bool matchmbl = false;
      hists_["mt2_220"][type]->Fill( ev->mt2_220, ev->weight );
      if( ev->mt2_210 > 1 ) hists_["mt2_210"][type]->Fill( ev->mt2_210, ev->weight );

      for( unsigned int m=0; m < ev->mbls.size(); m++ ){
         hists_["mbl"][type]->Fill( ev->mbls[m], ev->weight );
         hists_["mbl_fit"][type]->Fill( ev->mbls[m], ev->weight );
         //hists2d_["mblV220"][type]->Fill( ev->mt2_220, ev->mbls[m], ev->weight );
         if( ev->mbls[m] == ev->mt2_220 ) matchmbl = true;
      }

      if( matchmbl){
         hists_["mt2_220_matchmbl"][type]->Fill( ev->mt2_220, ev->weight ); 
      }
      else{
         hists_["mt2_220_nomatchmbl"][type]->Fill( ev->mt2_220, ev->weight );
      }

      //
      // kinematic distributions
      //

      double bb_dR = fabs(jet1.DeltaR(jet2));
      double ll_dR = fabs(lep1.DeltaR(lep2));
      double bl_dR1 = fabs(jet1.DeltaR(lep1));
      double bl_dR2 = fabs(jet2.DeltaR(lep1));
      double bl_dR3 = fabs(jet1.DeltaR(lep2));
      double bl_dR4 = fabs(jet2.DeltaR(lep2));

      // pt, eta, phi
      hists_["b_pt"][type]->Fill( jet1.Pt(), ev->weight );
      hists_["b_pt"][type]->Fill( jet2.Pt(), ev->weight );
      hists_["b_eta"][type]->Fill( jet1.Eta(), ev->weight );
      hists_["b_eta"][type]->Fill( jet2.Eta(), ev->weight );
      hists_["b_phi"][type]->Fill( jet1.Phi(), ev->weight );
      hists_["b_phi"][type]->Fill( jet2.Phi(), ev->weight );
      hists_["l_pt"][type]->Fill( lep1.Pt(), ev->weight );
      hists_["l_pt"][type]->Fill( lep2.Pt(), ev->weight );
      hists_["l_eta"][type]->Fill( lep1.Eta(), ev->weight );
      hists_["l_eta"][type]->Fill( lep2.Eta(), ev->weight );
      hists_["l_phi"][type]->Fill( lep1.Phi(), ev->weight );
      hists_["l_phi"][type]->Fill( lep2.Phi(), ev->weight );
      hists_["bb_pt"][type]->Fill( bb.Pt(), ev->weight );
      hists_["bb_eta"][type]->Fill( bb.Rapidity(), ev->weight );
      hists_["bb_phi"][type]->Fill( bb.Phi(), ev->weight );
      hists_["ll_pt"][type]->Fill( ll.Pt(), ev->weight );
      hists_["ll_eta"][type]->Fill( ll.Rapidity(), ev->weight );
      hists_["ll_phi"][type]->Fill( ll.Phi(), ev->weight );
      hists_["bbll_pt"][type]->Fill( bbll.Pt(), ev->weight );
      hists_["bbll_eta"][type]->Fill( bbll.Rapidity(), ev->weight );
      hists_["bbll_phi"][type]->Fill( bbll.Phi(), ev->weight );

      // invariant mass
      hists_["bb_m"][type]->Fill( bb.M(), ev->weight );
      hists_["ll_m"][type]->Fill( ll.M(), ev->weight );
      hists_["bbll_m"][type]->Fill( bbll.M(), ev->weight );
      hists_["bll_m"][type]->Fill( bll1.M(), ev->weight );
      hists_["bll_m"][type]->Fill( bll2.M(), ev->weight );
      hists_["bbl_m"][type]->Fill( bbl1.M(), ev->weight );
      hists_["bbl_m"][type]->Fill( bbl2.M(), ev->weight );

      // dR
      hists_["bb_dR"][type]->Fill( bb_dR, ev->weight );
      hists_["ll_dR"][type]->Fill( ll_dR, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR1, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR2, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR3, ev->weight );
      hists_["bl_dR"][type]->Fill( bl_dR4, ev->weight );

      // met
      hists_["met"][type]->Fill( met.Pt(), ev->weight );
      hists_["met_x"][type]->Fill( met.Px(), ev->weight );
      hists_["met_y"][type]->Fill( met.Py(), ev->weight );

   }

   return;
}

void Fitter::PrintHists( map< string, map<string, TH1D*> >& hists_ ){
   cout << "Printing data/mc plots." << endl;

   typedef map<string, TH1D*> tmap;
   typedef map<string, tmap> hmap;

   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }

   TFile *fileout = new TFile( (pathstr+"/plotsDataMC.root").c_str(), "RECREATE" );
   fileout->cd();

   TDirectory *dall = fileout->mkdir( "all" );
   TDirectory *dkin = fileout->mkdir( "kinematics" );
   TDirectory *dmass = fileout->mkdir( "masses" );

   
   // print all hists, all channels
   dall->cd();
   for( hmap::iterator h = hists_.begin(); h != hists_.end(); h++){
      TDirectory *dtemp = dall->mkdir( (h->first).c_str() );
      dtemp->cd();
      for( tmap::iterator t = h->second.begin(); t != h->second.end(); t++){
         t->second->Write();
      }
   }
   

   // data/mc stack plots w/ mt = 172.5
   for(hmap::iterator h = hists_.begin(); h != hists_.end(); h++){

      string name = h->first;
      tmap hist = h->second;
      TH1D* hdata = (TH1D*)hist["data"]->Clone("data_cl");

      // formatting
      TH1D* httbar_signal = (TH1D*)hist["ttbar172_signal"]->Clone("ttbar172_signal_cl");
      TH1D* httbar_mistag = (TH1D*)hist["ttbar172_mistag"]->Clone("ttbar172_mistag_cl");
      TH1D* httbar_taus = (TH1D*)hist["ttbar172_taus"]->Clone("ttbar172_taus_cl");
      TH1D* httbar_hadronic = (TH1D*)hist["ttbar172_hadronic"]->Clone("ttbar172_hadronic_cl");
      TH1D* hother = (TH1D*)hist["other"]->Clone("other_cl");

      httbar_signal->SetFillColor( kAzure-9 );
      httbar_mistag->SetFillColor( kYellow-9 );
      httbar_taus->SetFillColor( kRed-10 );
      httbar_hadronic->SetFillColor( kCyan-10 );
      hother->SetFillColor( kGreen-10 );

      TH1D* hallmc = (TH1D*)httbar_signal->Clone("hallmc");
      hallmc->Add( httbar_mistag );
      hallmc->Add( httbar_taus );
      hallmc->Add( httbar_hadronic );
      hallmc->Add( hother );

      double norm = hdata->Integral("width") / hallmc->Integral("width");

      // renormalize distribution in mc for all processes
      httbar_signal->Scale( norm );
      httbar_mistag->Scale( norm );
      httbar_taus->Scale( norm );
      httbar_hadronic->Scale( norm );
      hother->Scale( norm );
      hallmc->Scale( norm );

      if( name.compare("mbl_fit") == 0 ){
         httbar_signal->Rebin(4);
         httbar_mistag->Rebin(4);
         httbar_taus->Rebin(4);
         httbar_hadronic->Rebin(4);
         hother->Rebin(4);
         hallmc->Rebin(4);
         hdata->Rebin(4);
         hdata->GetYaxis()->SetTitle("Events/10 GeV");
      }

      TCanvas * canvas = new TCanvas( ("c"+name).c_str(), ("c"+name).c_str(), 800, 800 );
      canvas->SetFillColor(0);
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

      pad1->cd();

      THStack * hstack = new THStack( ("hs"+name).c_str(), ("hs"+name).c_str() );
      hstack->Add( hother );
      hstack->Add( httbar_hadronic );
      hstack->Add( httbar_taus );
      hstack->Add( httbar_mistag );
      hstack->Add( httbar_signal );

      hstack->Draw("HIST");

      hstack->GetYaxis()->SetTitle( hdata->GetYaxis()->GetTitle() );
      hstack->SetMinimum( 0.5 );
      hstack->GetXaxis()->SetTitleSize(0.00);
      hstack->GetYaxis()->SetLabelSize(0.07);
      hstack->GetYaxis()->SetTitleSize(0.08);
      hstack->GetYaxis()->SetTitleOffset(1.0);
      hstack->GetXaxis()->SetLabelFont(42);
      hstack->GetYaxis()->SetLabelFont(42);
      hstack->GetXaxis()->SetTitleFont(42);
      hstack->GetYaxis()->SetTitleFont(42);

      hdata->SetMarkerStyle(20);
      hdata->Draw( "same EP" );

      TLegend * legend = new TLegend(0.717,0.650,0.874,0.870);
      legend->AddEntry( hdata, "data" );
      legend->AddEntry( httbar_signal, "signal", "f" );
      legend->AddEntry( httbar_mistag, "mistag bkg", "f" );
      legend->AddEntry( httbar_taus, "tau decays", "f" );
      legend->AddEntry( httbar_hadronic, "hadronic decays", "f" );
      legend->AddEntry( hother, "non-ttbar bkg", "f" );

      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->Draw( "same" );

      pad2->cd();

      TH1D *hratio = (TH1D*)hdata->Clone("hratio");
      hratio->Divide( hallmc );
      hratio->SetTitle(";"+TString(hdata->GetXaxis()->GetTitle())+";data/mc");
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

      TF1 *func = new TF1("func","[0]",-10E6,10E6);
      func->SetParameter(0,1.0);
      func->SetLineWidth(1);
      func->SetLineStyle(7);
      func->SetLineColor(1);
      func->Draw("same");

      fileout->cd();
      if( name.find("220")==string::npos and name.find("221")==string::npos
            and name.find("mbl")==string::npos ){
         dkin->cd();
      }
      canvas->Write();

      pad1->cd();
      TLatex *text = new TLatex();
      text->SetTextSize(0.05);
      text->SetNDC();
      text->DrawLatex(0.4, 0.9, hdata->GetTitle() );
      canvas->Print( ("results/plots/stack_"+name+".pdf").c_str() );

      delete canvas;
      delete legend;
      delete func;
      delete text;
   }

   TCanvas *c220_mblstack = new TCanvas("c220_mblstack","c220_mblstack",700,700);
   c220_mblstack->cd();

   TH1D* h220_matchmbl = (TH1D*)hists_["mt2_220_matchmbl"]["data"]->Clone("mt2_220_matchmbl_cl");
   TH1D* h220_nomatchmbl = (TH1D*)hists_["mt2_220_nomatchmbl"]["data"]->Clone("mt2_220_nomatchmbl_cl");

   h220_matchmbl->Add( h220_nomatchmbl);

   h220_matchmbl->SetFillColor( kRed-9 );
   h220_nomatchmbl->SetFillColor( kAzure-9 );

   h220_matchmbl->Draw( "HIST" );
   h220_nomatchmbl->Draw( "same HIST" );

   TLegend *l220_mblstack = new TLegend(0.605528,0.655866,0.866834,0.816333);
   l220_mblstack->AddEntry( h220_matchmbl, "220 = mbl", "f" );
   l220_mblstack->AddEntry( h220_nomatchmbl, "220 #neq mbl", "f" );
   l220_mblstack->SetFillStyle(0);
   l220_mblstack->SetBorderSize(0);
   l220_mblstack->Draw("same");

   c220_mblstack->Write();

   //
   // chi2 plots
   //
   gErrorIgnoreLevel = 2000;
   for(hmap::iterator h = hists_.begin(); h != hists_.end(); h++){

      string hname = h->first;
      tmap hist = h->second;
      TH1D* hdata = (TH1D*)hist["data"]->Clone("data_cl");

      TGraph *gchi2 = new TGraph();

      // do this for all top masses
      double masspnts [] = {161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5};
      for(int i=0; i < 8; i++){

         stringstream dstr;
         dstr << floor(masspnts[i]);
         string dname = "ttbar"+dstr.str();

         TH1D *hallmc = (TH1D*)hist[dname+"_signal"]->Clone("hallmc");
         hallmc->Add( hist[dname+"_mistag"] );
         hallmc->Add( hist[dname+"_taus"] );
         hallmc->Add( hist[dname+"_hadronic"] );
         hallmc->Add( hist["other"] );

         double chi2 = hdata->Chi2Test(hallmc,"UW CHI2");

         gchi2->SetPoint(i, masspnts[i], chi2);

         // write out mass hists
         if( hname.find("mbl") != string::npos or hname.find("mt2") != string::npos ){
            dmass->cd();
            hallmc->SetName( (hname+"_"+dname).c_str() );
            hallmc->Write();
            fileout->cd();
         }

      }

      TCanvas * canvas = new TCanvas( ("g"+hname).c_str(), ("g"+hname).c_str(), 800, 800 );
      canvas->SetFillColor(0);
      canvas->cd();

      gchi2->SetName(("g"+hname).c_str());
      gchi2->SetMarkerStyle(20);
      gchi2->SetTitle(hname.c_str());
      gchi2->GetXaxis()->SetTitle("Top Mass (GeV)");
      gchi2->GetYaxis()->SetTitle("#chi^{2}");

      bool yrange = false;
      string names [] = {"b_pt","bb_pt","bbl_m", "bbll_m", "mbl", 
         "mt2_220","mt2_221","mt2_220_nomatchmbl"};
      for(unsigned int i=0; i < sizeof(names)/sizeof(names[0]); i++){
         if( hname.compare(names[i]) == 0 ) yrange = true;
      }

      if( yrange ){
         gchi2->SetMaximum(4000);
         gchi2->SetMinimum(0);
      }

      gchi2->Draw("AP");

      TF1 *func = new TF1("func","pol2",160,182);
      gchi2->Fit(func,"Q");

      gStyle->SetOptFit();

      double pb = func->GetParameter(1);
      double pa = func->GetParameter(2);
      double minimum = -pb/(2*pa);

      cout << hname << " minimum = " << minimum << endl;

      stringstream minstr;
      minstr << minimum;
      TLatex *text = new TLatex();
      text->SetTextSize(0.04);
      text->SetNDC();
      text->DrawLatex(0.5, 0.7, ("min @ "+minstr.str()+" GeV").c_str() );

      fileout->cd();
      if( hname.find("220")==string::npos and hname.find("221")==string::npos
            and hname.find("mbl")==string::npos ){
         dkin->cd();
      }
      canvas->Write();

      text->SetTextSize(0.05);
      text->DrawLatex(0.4, 0.9, hdata->GetTitle() );
      canvas->Print( ("results/plots/chi2_"+hname+".pdf").c_str() );


      delete gchi2;
      delete func;
      delete canvas;
      delete text;
   }
   fileout->Close();

   return;
}

void Fitter::PlotTemplates( map< string, map<string, TH1D*> >& hists_ ){

   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }

   TFile *fileout = new TFile( (pathstr+"/plotsTemplates.root").c_str(), "RECREATE" );
   fileout->cd();

   // templates for all masses
   string names [] = {"mbl","mbl_fit"};
   string sb[] = {"sig","bkg"};
   double masspnts [] = {161.5,163.5,166.5,169.5,172.5,175.5,178.5,181.5};
   for(unsigned int i=0; i < sizeof(names)/sizeof(names[0]); i++){ // distributions

      TDirectory *dir = fileout->mkdir( names[i].c_str() );
      dir->cd();

      for(unsigned int k=0; k < sizeof(sb)/sizeof(sb[0]); k++){ // sig,bkg
         for(int j=0; j < 8; j++){ // masses

            stringstream ssmass;
            ssmass << floor(masspnts[j]);
            string smass = ssmass.str();

            TCanvas *canvas = new TCanvas( ("c"+sb[k]+"_"+names[i]+smass).c_str(),
                  ("c"+sb[k]+"_"+names[i]+smass).c_str(), 800, 800);
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

            TH1D *hmc;
            if( sb[k] == "sig" ){
               hmc = (TH1D*)hists_[names[i]]["ttbar"+smass+"_signal"]->Clone("hmc");
            }else{
               hmc = (TH1D*)hists_[names[i]]["ttbar"+smass+"_mistag"]->Clone("hmc");
               hmc->Add( hists_[names[i]]["ttbar"+smass+"_taus"] );
               hmc->Add( hists_[names[i]]["ttbar"+smass+"_hadronic"] );
               hmc->Add( hists_[names[i]]["other"] );
            }

            // set bkg control sample
            TH1D *hmc_bkg;
            hmc_bkg = (TH1D*)hists_[names[i]]["ttbar"+smass+"_bkgcontrol_signal"]->Clone("hmc_bkg");
            hmc_bkg->Add( hists_[names[i]]["ttbar"+smass+"_bkgcontrol_mistag"] );
            hmc_bkg->Add( hists_[names[i]]["ttbar"+smass+"_bkgcontrol_taus"] );
            hmc_bkg->Add( hists_[names[i]]["ttbar"+smass+"_bkgcontrol_hadronic"] );
            hmc_bkg->Add( hists_[names[i]]["other_bkgcontrol"] );

            if( names[i].compare("mbl_fit") == 0 ){
               hmc->Rebin(4);
               hmc->GetYaxis()->SetTitle("Events/10 GeV");
            }

            hmc->SetTitle( hmc->GetTitle()+TString(" "+sb[k]+" shape @ "+smass+".5") );
            hmc->GetXaxis()->SetTitleSize(0.00);
            hmc->GetYaxis()->SetLabelSize(0.07);
            hmc->GetYaxis()->SetTitleSize(0.08);
            hmc->GetYaxis()->SetTitleOffset(1.0);
            hmc->GetXaxis()->SetLabelFont(42);
            hmc->GetYaxis()->SetLabelFont(42);
            hmc->GetXaxis()->SetTitleFont(42);
            hmc->GetYaxis()->SetTitleFont(42);

            hmc->SetMarkerStyle(20);
            hmc->DrawCopy();

            Shapes * fptr = new Shapes( hmc_bkg );
            fptr->aGPsig.ResizeTo( aGPsig.GetNoElements() );
            fptr->aGPsig = aGPsig;
            fptr->aGPbkg.ResizeTo( aGPbkg.GetNoElements() );
            fptr->aGPbkg = aGPbkg;
            TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Fmbl_tot, 0, rangembl, 5);

            // normalization inside likelihood function (temp)
            ftemplate->SetParameters( masspnts[j], 1-k, 1.0, 1.0, 1.0 );
            double integralsig = (sb[k] == "sig") ? ftemplate->Integral(0,rangembl) : 1.0;
            double integralbkg = (sb[k] == "bkg") ? ftemplate->Integral(0,rangembl) : 1.0;
            ftemplate->SetParameters( masspnts[j], 1-k,
                  hmc->Integral("width"), integralsig, integralbkg );

            ftemplate->SetLineWidth(2);
            ftemplate->DrawCopy("same");
            hmc->DrawCopy("same"); // redraw points

            // pad 2
            pad2->cd();
            TH1D *hratio = (TH1D*)hmc->Clone("hratio");
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
            delete ftemplate;
            delete fptr;
         }
      }
   }
   fileout->cd();
   
   // plot template as a function of top mass
   /*
   TDirectory *dir = fileout->mkdir( "mtshape" );
   dir->cd();
   for(unsigned int k=0; k < sizeof(sb)/sizeof(sb[0]); k++){ // sig,bkg
      for(double x=0; x <= rangembl; x+=10){ // bin of mbl

         stringstream ssmbl;
         ssmbl << x;
         string smbl = ssmbl.str();

         TCanvas *canvas = new TCanvas( ("c"+sb[k]+"_mbl"+smbl).c_str(),
               ("c"+sb[k]+"_mbl"+smbl).c_str(), 800, 800);
         canvas->SetFillColor(0);
         canvas->cd();

         // set bkg control sample (dummy)
         TH1D *hmc_bkg;
         hmc_bkg = (TH1D*)hists_["mbl"]["ttbar172_bkgcontrol_signal"]->Clone("hmc_bkg");
         hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_mistag"] );
         hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_taus"] );
         hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_hadronic"] );
         hmc_bkg->Add( hists_["mbl"]["other_bkgcontrol"] );

         // graph with template value at mbl = x
         TGraph *gtemplate = new TGraph();
         Shapes * fptr = new Shapes( hmc_bkg );
         fptr->aGPsig.ResizeTo( aGPsig.GetNoElements() );
         fptr->aGPsig = aGPsig;
         fptr->aGPbkg.ResizeTo( aGPbkg.GetNoElements() );
         fptr->aGPbkg = aGPbkg;
         TF1 *ftemplate = new TF1("ftemplate", fptr, &Shapes::Fmbl_tot, 0, rangembl, 5);
         int count=0;
         for(double m=160.0; m <= 183.0; m+=1){ // value of mt
            // normalization inside likelihood function (temp)
            ftemplate->SetParameters( m, 1-k, 1.0, 1.0, 1.0 );
            double integralsig = (sb[k] == "sig") ? ftemplate->Integral(0,rangembl) : 1.0;
            double integralbkg = (sb[k] == "bkg") ? ftemplate->Integral(0,rangembl) : 1.0;
            ftemplate->SetParameters( m, 1-k, 1.0, integralsig, integralbkg );

            gtemplate->SetPoint(count, m, ftemplate->Eval(x));
            count++;
         }
         gtemplate->SetTitle( TString("M_{bl} "+sb[k]+" shape @ mbl = "+smbl) );
         gtemplate->SetLineColor(2);
         gtemplate->SetLineWidth(2);

         // now do the same at mc masspoints
         TGraphErrors *gmc = new TGraphErrors();
         count=0;
         for(int j=0; j < 8; j++){
            stringstream ssmass;
            ssmass << floor(masspnts[j]);
            string smass = ssmass.str();

            TH1D *hmc;
            if( sb[k] == "sig" ){
               hmc = (TH1D*)hists_["mbl"]["ttbar"+smass+"_signal"]->Clone("hmc");
            }else{
               hmc = (TH1D*)hists_["mbl"]["ttbar"+smass+"_mistag"]->Clone("hmc");
               hmc->Add( hists_["mbl"]["ttbar"+smass+"_taus"] );
               hmc->Add( hists_["mbl"]["ttbar"+smass+"_hadronic"] );
               hmc->Add( hists_["mbl"]["other"] );
            }
            hmc->Scale( 1.0/hmc->Integral("width") );

            gmc->SetPoint(count, masspnts[j], hmc->GetBinContent(hmc->FindBin(x)) );
            gmc->SetPointError(count, 0.0, hmc->GetBinError(hmc->FindBin(x)) );
            count++;
         }

         gmc->SetMarkerStyle(20);

         gtemplate->SetMinimum( min(gtemplate->GetMinimum(),gmc->GetMinimum()) );
         gtemplate->SetMaximum( max(gtemplate->GetMaximum(),gmc->GetMaximum()) );
         gtemplate->Draw("AC");
         gmc->Draw("P");

         canvas->Write();

         delete canvas;
         delete ftemplate;
         delete fptr;
      }
   }
*/
   fileout->cd();

   
   TCanvas *cmbl_signal = new TCanvas("cmbl_signal","M_{bl} Template",800,800);
   cmbl_signal->cd();

   // mass points
   TH1D* mbl161 = (TH1D*)hists_["mbl"]["ttbar161_signal"]->Clone("mbl161");
   TH1D* mbl172 = (TH1D*)hists_["mbl"]["ttbar172_signal"]->Clone("mbl172");
   TH1D* mbl181 = (TH1D*)hists_["mbl"]["ttbar181_signal"]->Clone("mbl181");

   mbl161->Scale( 1.0/mbl161->Integral("width") );
   mbl172->Scale( 1.0/mbl172->Integral("width") );
   mbl181->Scale( 1.0/mbl181->Integral("width") );

   mbl161->SetLineColor(2);
   mbl172->SetLineColor(1);
   mbl181->SetLineColor(3);
   mbl161->SetMarkerColor(2);
   mbl172->SetMarkerColor(1);
   mbl181->SetMarkerColor(3);

   mbl161->DrawCopy();
   mbl172->DrawCopy("same");
   mbl181->DrawCopy("same");

   // mbl likelihood

   // set bkg control sample (dummy)
   TH1D *hmc_bkg;
   hmc_bkg = (TH1D*)hists_["mbl"]["ttbar172_bkgcontrol_signal"]->Clone("hmc_bkg");
   hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_mistag"] );
   hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_taus"] );
   hmc_bkg->Add( hists_["mbl"]["ttbar172_bkgcontrol_hadronic"] );
   hmc_bkg->Add( hists_["mbl"]["other_bkgcontrol"] );

   Shapes * fptr = new Shapes( hmc_bkg );
   fptr->aGPsig.ResizeTo( aGPsig.GetNoElements() );
   fptr->aGPsig = aGPsig;
   fptr->aGPbkg.ResizeTo( aGPbkg.GetNoElements() );
   fptr->aGPbkg = aGPbkg;
   TF1 *fmbl_tot = new TF1("fmbl_tot", fptr, &Shapes::Fmbl_tot, 0, rangembl, 5);

   fmbl_tot->SetParameters( 161.5, 1.0, 1.0, 1.0, 1.0 );
   fmbl_tot->SetParameters( 161.5, 1.0, 1.0, fmbl_tot->Integral(0,rangembl), 1.0 );
   fmbl_tot->SetLineColor(2);
   fmbl_tot->DrawCopy("same");

   fmbl_tot->SetParameters( 172.5, 1.0, 1.0, 1.0, 1.0 );
   fmbl_tot->SetParameters( 172.5, 1.0, 1.0, fmbl_tot->Integral(0,rangembl), 1.0 );
   fmbl_tot->SetLineColor(1);
   fmbl_tot->DrawCopy("same");

   fmbl_tot->SetParameters( 181.5, 1.0, 1.0, 1.0, 1.0 );
   fmbl_tot->SetParameters( 181.5, 1.0, 1.0, fmbl_tot->Integral(0,rangembl), 1.0 );
   fmbl_tot->SetLineColor(3);
   fmbl_tot->DrawCopy("same");

   TLegend *lm = new TLegend(0.61,0.59,0.83,0.77);
   lm->SetFillStyle(0);
   lm->SetBorderSize(0);
   lm->AddEntry( mbl161, "161.5" );
   lm->AddEntry( mbl172, "172.5" );
   lm->AddEntry( mbl181, "181.5" );
   lm->Draw("same");

   cmbl_signal->Write();

   fileout->Close();

   delete fptr;

   return;
}

