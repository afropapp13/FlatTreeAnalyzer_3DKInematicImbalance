#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

void create_spline() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			

	//------------------------------//

	vector<TString> plot_name; vector<TString> label;
	
	plot_name.push_back("TrueFineBinEvPlot"); label.push_back("");
	plot_name.push_back("TrueFineBinECalPlot"); label.push_back("ECal");

	int nplots = plot_name.size();

	//------------------------------//

	// input / output files

	TFile* f_bnb = new TFile("OutputFiles/FlatTreeAnalyzerOutput_GENIE_v3_0_6.root","readonly");
	TFile* f_honda = new TFile("OutputFiles/FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_Honda.root","readonly");
	TFile* f_spline = new TFile("OutputFiles/spline.root","recreate");

	//------------------------------//

	for (int iplot = 0; iplot < nplots; iplot++) {

		TH1D* h_bnb = (TH1D*)(f_bnb->Get(plot_name[iplot]));
		TH1D* h_honda = (TH1D*)(f_honda->Get(plot_name[iplot]));
	
		TH1D* spline = (TH1D*)(h_honda->Clone());
		spline->Divide(h_bnb);

		TString CanvasName = "HondaSpline" + label[iplot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.15);
		PlotCanvas->SetLeftMargin(0.17);
		PlotCanvas->SetRightMargin(0.05);
		PlotCanvas->SetBottomMargin(0.16);		
		PlotCanvas->Draw();	

		spline->SetLineWidth(3);
		spline->SetLineColor( kRed + 1 );	
		spline->SetLineStyle( kSolid );	

		spline->GetXaxis()->SetTitleFont(FontStyle);
		spline->GetXaxis()->SetLabelFont(FontStyle);
		spline->GetXaxis()->SetNdivisions(8);
		spline->GetXaxis()->SetLabelSize(TextSize);
		spline->GetXaxis()->SetTitleSize(TextSize);	
		spline->GetXaxis()->SetTitleOffset(1.1);					
		spline->GetXaxis()->SetLabelOffset(0.01);					
		spline->GetXaxis()->CenterTitle();						

		spline->GetYaxis()->SetTitleFont(FontStyle);
		spline->GetYaxis()->SetLabelFont(FontStyle);
		spline->GetYaxis()->SetNdivisions(6);
		spline->GetYaxis()->SetLabelSize(TextSize);
		spline->GetYaxis()->SetTitle("Honda/BNB");
		spline->GetYaxis()->SetTitleSize(TextSize);
		spline->GetYaxis()->SetTitleOffset(1.2);
		spline->GetYaxis()->CenterTitle();	

		PlotCanvas->cd();
		spline->Draw("hist");

		PlotCanvas->SaveAs("/exp/uboone/data/users/apapadop/PeLEETuples_Atmospherics/FlatTTreePlots/"+CanvasName+".pdf");
		delete PlotCanvas;
		
		f_spline->cd();
		spline->Write();
	
	}

	//------------------------------//

f_spline->cd();
	f_spline->Close();

	//------------------------------//


} // End of the program

