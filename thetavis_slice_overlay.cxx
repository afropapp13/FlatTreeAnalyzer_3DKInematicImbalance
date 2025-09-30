#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
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

void thetavis_slice_overlay(TString Tag = "") {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			

	TString OutFilePath = "OutputFiles/";

	//------------------------------//

	// input file

	TString name = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; 
	TFile* f = new TFile(name,"readonly");;

	//------------------------------//

	// plots

	// TH1D* lowe = (TH1D*)(f->Get("TrueFineBinLowEThetaVisPlot"));
	// TH1D* mide = (TH1D*)(f->Get("TrueFineBinMidEThetaVisPlot"));
	// TH1D* highe = (TH1D*)(f->Get("TrueFineBinHighEThetaVisPlot"));	
	
	TH1D* lowe = (TH1D*)(f->Get("TrueFineBinLowErecoThetaVisPlot"));
	TH1D* mide = (TH1D*)(f->Get("TrueFineBinMidErecoThetaVisPlot"));
	TH1D* highe = (TH1D*)(f->Get("TrueFineBinHighErecoThetaVisPlot"));		

	//------------------------------//

	// canvas

	TString CanvasName = "thetavis_overlay";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetTopMargin(0.15);
	PlotCanvas->SetLeftMargin(0.17);
	PlotCanvas->SetRightMargin(0.05);
	PlotCanvas->SetBottomMargin(0.16);		
	PlotCanvas->Draw();	

	//------------------------------//
	
	// legend

	TLegend* leg = new TLegend(0.55,0.6,0.75,0.8);
	leg->SetBorderSize(0);
	leg->SetNColumns(1);
	leg->SetTextSize(TextSize);	
	leg->SetTextFont(FontStyle);						
	leg->SetMargin(0.1);						

	//------------------------------//

	lowe->SetLineWidth(3);
	lowe->SetLineColor( kAzure+7 );	
	lowe->SetLineStyle( kSolid );	

	lowe->GetXaxis()->SetTitleFont(FontStyle);
	lowe->GetXaxis()->SetLabelFont(FontStyle);
	lowe->GetXaxis()->SetNdivisions(8);
	lowe->GetXaxis()->SetLabelSize(TextSize);
	lowe->GetXaxis()->SetTitleSize(TextSize);	
	lowe->GetXaxis()->SetTitleOffset(1.1);					
	lowe->GetXaxis()->CenterTitle();						

	lowe->GetYaxis()->SetTitleFont(FontStyle);
	lowe->GetYaxis()->SetLabelFont(FontStyle);
	lowe->GetYaxis()->SetNdivisions(6);
	lowe->GetYaxis()->SetLabelSize(TextSize);
	lowe->GetYaxis()->SetTitle("probability density");
	lowe->GetYaxis()->SetTitleSize(TextSize);
	lowe->GetYaxis()->SetTitleOffset(1.2);
	lowe->GetYaxis()->CenterTitle();		

	lowe->Scale(1./lowe->Integral("width"));
	mide->Scale(1./mide->Integral("width"));
	highe->Scale(1./highe->Integral("width"));

	lowe->GetYaxis()->SetRangeUser(0,1.2*highe->GetMaximum());
	lowe->Draw("hist same");

	TLegendEntry* lowlegColor = leg->AddEntry(lowe,"E_{reco} < 0.5 GeV","");
	lowlegColor->SetTextColor( kAzure+7 ); 

	//------------------------------//	

	mide->SetLineWidth(3);
	mide->SetLineColor(kOrange+7);
	mide->Draw("hist same");

	TLegendEntry* midlegColor = leg->AddEntry(mide,"0.5 < E_{reco} < 0.8 GeV","");
	midlegColor->SetTextColor( kOrange+7 ); 	

	//------------------------------//	

	highe->SetLineWidth(3);
	highe->SetLineColor(kGreen+1);
	highe->Draw("hist same");	

	TLegendEntry* highlegColor = leg->AddEntry(highe,"E_{reco} > 0.8 GeV","");
	highlegColor->SetTextColor( kGreen+1 ); 

	//------------------------------//	

	leg->Draw();

	PlotCanvas->SaveAs(CanvasName+".pdf");
	delete PlotCanvas;

	//------------------------------//

} // End of the program
