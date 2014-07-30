#include "Plotting.h"
#include "TopMass.h"
#include "Mt2Calculator.h"

#include "TNtuple.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

Plotting::Plotting(){
}

Plotting::~Plotting(){
}

void Plotting::MaosNan( vector<Event> eventvec){

	TFile *fileout = new TFile( "results/maosnans.root", "RECREATE" );

	string names [] = {"220","210"};

	for (unsigned int j=0; j < sizeof(names)/sizeof(names[0]); j++){

		TNtuple *metdif = new TNtuple("metdifPlot", ("MAOS neutrinos from M_{T2}" + names[j] + " that have Pz of Nan").c_str(), "metdif:mtransdif");
		TH1F *metdifhist = new TH1F("metdifhist", ("MAOS neutrinos from M_{T2} " + names[j] + ", red=Nan, black=non-Nan;abs(MET - nGEN - nbarGEN, pT) (GeV);Events/1.5GeV,normalised").c_str(), 100, 0, 350);
		TH1F *metdifhistnan = new TH1F("metdifhistnan", ("MAOS neutrinos from M_{T2} " + names[j] + ";abs(MET - nGEN - nbarGEN, pT) (GeV);Events/1.5GeV").c_str(), 100, 0, 350);
		vector <float> metdifx;
		vector <float> metdify;
		metdifx.clear();
		metdify.clear();

		int nancount = 0;
		int nonnancount = 0;

		for( vector<Event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){

			vector<TLorentzVector> neuarray;
			neuarray.clear();
			if( j==1 ){
				neuarray.push_back( ev->maos210_neutrino1p);
				neuarray.push_back( ev->maos210_neutrino1m);
				neuarray.push_back( ev->maos210_neutrino2p);
				neuarray.push_back( ev->maos210_neutrino2m);
			}
			else if( j==0 ){
				neuarray.push_back( ev->maos220_neutrino1ap);
				neuarray.push_back( ev->maos220_neutrino1am);
				neuarray.push_back( ev->maos220_neutrino2ap);
				neuarray.push_back( ev->maos220_neutrino2am);
				neuarray.push_back( ev->maos220_neutrino1bp);
				neuarray.push_back( ev->maos220_neutrino1bm);
				neuarray.push_back( ev->maos220_neutrino2bp);
				neuarray.push_back( ev->maos220_neutrino2bm);
			}

			for (unsigned int i=0; i < neuarray.size(); i++){

				TLorentzVector lepton;
				if (j==1){
					if (i==0 or i==1 or i==4 or i==5){ lepton = ev->lep1; }
					else { lepton = ev->lep2; }
				} else if (j==0){
					if (i==0 or i==1){ lepton = ev->lep1; }
					else { lepton = ev->lep2; }
				}

				if ( isnan(neuarray[i].Pz()) ){

					Mt2Calculator::Calculator calc;
					double mTrans = calc.CalcTransverseMass(lepton, neuarray[i], 0);

					metdifx.push_back( abs(( ev->met - ev->nGEN - ev->nbarGEN ).Pt()) );
					metdify.push_back( abs(( mTrans - 80.4 )) );

					metdif->Fill( abs(( ev->met - ev->nGEN - ev->nbarGEN ).Pt()) , abs(( mTrans - 80.4 )) );

					metdifhistnan->Fill( abs(( ev->met - ev->nGEN - ev->nbarGEN ).Pt()) );

					nancount++;

				} else {

					metdifhist->Fill( abs(( ev->met - ev->nGEN - ev->nbarGEN ).Pt()) );

					nonnancount++;

				}
			}

		}

/*		float metdifxarr [metdifx.size()];
		float metdifyarr [metdify.size()];

		for (unsigned int i=0; i < metdifx.size(); i++){
			metdifxarr[i] = metdifx[i];
			metdiryarr[i] = metdify[i];
		}
*/
		TGraph *graphMetdif = new TGraph( metdifx.size(), &(metdifx[0]), &(metdify[0]) );
		graphMetdif->SetMarkerColor(1);
		graphMetdif->SetMarkerStyle(2);
		graphMetdif->SetTitle(("MAOS neutrinos from M_{T2}" + names[j] + " that have Pz of Nan").c_str());
		graphMetdif->GetXaxis()->SetTitle("abs(MET-nGEN-nbarGEN, pT)");
		graphMetdif->GetYaxis()->SetTitle("abs( (lep+MaosNeu, mTrans) - 80.4)"); 
//		graphMetdif->GetXaxis()->SetTitleSize(100.0);
//		graphMetdif->GetYaxis()->SetTitleSize(100.0);

		TLine *line = new TLine(0,0,400,400);
		line->SetLineColor(2);
		line->SetLineWidth(1);
		line->SetLineStyle(2);

		metdifhist->SetLineColor(1);
		metdifhistnan->SetLineColor(2);

		metdifhist->Scale( (metdifhistnan->Integral() )/(metdifhist->Integral() ) );

		TCanvas *canvasMetdif = new TCanvas ( ("metdif_"+names[j]).c_str(),("metdif_"+names[j]).c_str(),700,700);
		TCanvas *canvasMetdifhist = new TCanvas ( ("metdifhist_"+names[j]).c_str(),("metdifhist_"+names[j]).c_str(),700,700);
		canvasMetdifhist->SetTitle(("MAOS neutrinos from M_{T2} " + names[j] + ", red = Nan, black = non-Nan").c_str()); 

		metdif->SetMarkerColor(1);
		metdif->SetMarkerStyle(2);

		fileout->cd();
		canvasMetdif->cd();
//		metdif->Draw("mtransdif:metdif");
		graphMetdif->Draw("AP");
		line->Draw("SAME");
		canvasMetdif->Write();
		canvasMetdif->ls();

		canvasMetdifhist->cd();
		metdifhist->Draw("HIST");
		metdifhistnan->Draw("HIST SAME");
		canvasMetdifhist->Write();
		canvasMetdifhist->ls();

		cout << "nancount = " << nancount << endl;
		cout << "nonnancount = " << nonnancount << endl;

		delete canvasMetdif;
		delete metdif;
		delete canvasMetdifhist;
		delete metdifhist;
		delete graphMetdif;
		delete line;

	}

	fileout->Close();

}
