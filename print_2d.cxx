#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>

#include "constants.h"

using namespace std;
using namespace constants;

void print_2d() {

	//	----------------------------------------//

	TH2D::SetDefaultSumw2();
	
	double TextSize = 0.07;
	const Int_t NCont = 999;	
	gStyle->SetPalette(55); 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);

	//	----------------------------------------//

	vector<TString> PlotNames;

	PlotNames.push_back("TrueTheta_vsThetaVisPlot"); 
	
	const int N2DPlots = PlotNames.size();
	cout << "Number of 2D Plots = " << N2DPlots << endl;

	//	----------------------------------------//
	
	vector<TString> NameOfSamples;
	NameOfSamples.push_back("atmo_AR23_20i_00_000");
	NameOfSamples.push_back("bnb_AR23_20i_00_000");
	const int NSamples = NameOfSamples.size();
		
	//	----------------------------------------//

	vector< vector<TH2D*> > Plots; Plots.clear(); Plots.resize(NSamples);

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		Plots[WhichSample].resize(N2DPlots);
		TFile* f = TFile::Open("output_files/analyzer_output_"+NameOfSamples[WhichSample]+".root","readonly");

		for (int WhichPlot = 0; WhichPlot < N2DPlots; WhichPlot ++) {

			//	----------------------------------------//

			Plots[WhichSample][WhichPlot] = (TH2D*)(f->Get(PlotNames[WhichPlot]));

			//	----------------------------------------//			
	
			TString PlotCanvasName = PlotNames[WhichPlot]+NameOfSamples[WhichSample];
			TCanvas* PlotCanvas = new TCanvas(PlotCanvasName,PlotCanvasName,205,34,1024,768);
			PlotCanvas->cd();
			PlotCanvas->SetBottomMargin(0.16);
			PlotCanvas->SetLeftMargin(0.15);
			PlotCanvas->SetRightMargin(0.15);				
					
			gStyle->SetMarkerSize(1.5);
			gStyle->SetPaintTextFormat("4.2f");				
				
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);				
			Plots[WhichSample][WhichPlot]->GetXaxis()->CenterTitle();
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetNdivisions(6);
			Plots[WhichSample][WhichPlot]->GetXaxis()->SetLabelOffset(0.01);				
					
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);				
			Plots[WhichSample][WhichPlot]->GetYaxis()->CenterTitle();
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetNdivisions(6);
			Plots[WhichSample][WhichPlot]->GetYaxis()->SetTitleOffset(1.);				
									
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelFont(FontStyle);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetLabelSize(TextSize);
			Plots[WhichSample][WhichPlot]->GetZaxis()->SetNdivisions(6);				

			Plots[WhichSample][WhichPlot]->SetMarkerColor(kWhite);				
			Plots[WhichSample][WhichPlot]->SetMarkerSize(0.9);
			Plots[WhichSample][WhichPlot]->Draw("colz"); 				

			PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/PeLEETuples_Atmospherics/dune_FlatTTreePlots/"+PlotCanvasName+".pdf");

			delete PlotCanvas;		
					
		} // End of the loop over the plots
		
		f->Close();

	} // End of the loop over the samples

} // End of the program