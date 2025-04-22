#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <THStack.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

#include "../../myClasses/Constants.h"

using namespace std;
using namespace Constants;

int DUNEColor = kBlack;
int T2KColor = kRed+1;

// -------------------------------------------------------------------------------------------------------------------------------------

void PrettyPlot(TH1D* h, int Color,bool shift) {

	h->SetLineColor(Color);
	h->SetLineWidth(3);
	//h->Scale(1./h->GetMaximum());

	h->GetXaxis()->SetRangeUser(0.,2.);
	
	h->GetXaxis()->SetNdivisions(8);
	h->GetXaxis()->CenterTitle();
	h->GetXaxis()->SetTitleFont(FontStyle);
	h->GetXaxis()->SetLabelFont(FontStyle);
	h->GetXaxis()->SetTitleSize(TextSize);
	h->GetXaxis()->SetLabelSize(TextSize);
	h->GetXaxis()->SetTitleOffset(0.9);
	h->GetXaxis()->SetTickSize(0.02);
	h->GetXaxis()->SetTitle("E_{#nu} [GeV]");

	h->GetYaxis()->SetRangeUser(0.,17.9);
	h->GetYaxis()->SetNdivisions(5);
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleFont(FontStyle);
	h->GetYaxis()->SetLabelFont(FontStyle);
	h->GetYaxis()->SetTitleSize(TextSize);
	h->GetYaxis()->SetLabelSize(TextSize);
	h->GetYaxis()->SetTickSize(0.02);
	h->GetYaxis()->SetTitleOffset(0.65);
	h->GetYaxis()->SetTitle("#nu flux pdf");

}

// -------------------------------------------------------------------------------------------------------------------------------------

void fluxes() {

	// -------------------------------------------------------------------------------------------------------------------------------------

	gStyle->SetOptStat(0);
	TGaxis::SetMaxDigits(3);

	gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(TextSize,"t"); gStyle->SetTitleFont(FontStyle,"t");

	// --------------------------------------------------------------------------------------------------------------------------------	

	TString CanvasName = "FluxPdfCanvas";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->SetBottomMargin(0.15);

	TString SplineCanvasName = "SplineCanvas";
	TCanvas* SplinePlotCanvas = new TCanvas(SplineCanvasName,SplineCanvasName,205,34,1024,768);
	SplinePlotCanvas->SetBottomMargin(0.15);	

	// ---------------------------------------------------------------------------------------------------------------------------

	TFile* bnb_file = TFile::Open("MCC9_FluxHist_volTPCActive.root");
	TH1D* bnb_flux = (TH1D*)( bnb_file->Get("hEnumu_cv") );

	//TFile* honda_file = TFile::Open("histogram.root","readonly");
	TFile* honda_file = TFile::Open("honda_flux_solmax.root","readonly");
	TH1D* honda_flux = (TH1D*)(honda_file->Get("h"));

	TFile* min_honda_file = TFile::Open("honda_flux_solmin.root","readonly");
	TH1D* min_honda_flux = (TH1D*)(min_honda_file->Get("h"));	

	TFile* old_honda_file = TFile::Open("histogram.root","readonly");
	TH1D* old_honda_flux = (TH1D*)(old_honda_file->Get("h"));	

	TFile* f_spline = new TFile("honda_spline.root","recreate");

	// --------------------------------------------------------------------------------------------------

	PrettyPlot(bnb_flux,DUNEColor,false);
	PrettyPlot(honda_flux,T2KColor,false);
	PrettyPlot(min_honda_flux,kAzure+7,false);
	PrettyPlot(old_honda_flux,kMagenta,false);	

	//--------------------//

	// BNB pdf
	double bnb_int_flux = bnb_flux->Integral("width");
	bnb_flux->Scale(1./bnb_int_flux);

	cout << bnb_flux->GetBinWidth(1) << endl;
	
	// max Honda pdf
	double honda_int_flux = honda_flux->Integral("width");	
	honda_flux->Scale(1./honda_int_flux);

	// min Honda pdf
	double honda_int_flux_min = min_honda_flux->Integral("width");	
	min_honda_flux->Scale(1./honda_int_flux_min);	

	// old Honda pdf
	double honda_int_flux_old = old_honda_flux->Integral("width");	
	old_honda_flux->Scale(1./honda_int_flux_old);		

	TH1D* average = (TH1D*)(honda_flux->Clone());
	average->Add(min_honda_flux);
	average->Scale(0.5);

	// --------------------------------------------------------------------------------------------------

	PlotCanvas->cd();
	average->Draw("hist same");
	bnb_flux->Draw("hist same");	

	TLatex *bnb = new TLatex(); 
	bnb->SetTextFont(FontStyle); 
	bnb->SetTextColor(DUNEColor); 
	bnb->SetTextSize(TextSize);
	bnb->DrawLatexNDC(0.6,0.82,"BNB");

	/*TLatex *Honda = new TLatex(); 
	Honda->SetTextFont(FontStyle); 
	Honda->SetTextColor(T2KColor); 
	Honda->SetTextSize(TextSize);
	Honda->DrawLatexNDC(0.6,0.75,"Honda max");*/

	TLatex *average_latex = new TLatex(); 
	average_latex->SetTextFont(FontStyle); 
	average_latex->SetTextColor(T2KColor); 
	average_latex->SetTextSize(TextSize);
	average_latex->DrawLatexNDC(0.6,0.75,"Honda");	

	/*TLatex *Honda_min = new TLatex(); 
	Honda_min->SetTextFont(FontStyle); 
	Honda_min->SetTextColor(kAzure+7); 
	Honda_min->SetTextSize(TextSize);
	Honda_min->DrawLatexNDC(0.6,0.68,"Honda min");	

	TLatex *Honda_old = new TLatex(); 
	Honda_old->SetTextFont(FontStyle); 
	Honda_old->SetTextColor(kMagenta); 
	Honda_old->SetTextSize(TextSize);
	Honda_old->DrawLatexNDC(0.6,0.61,"Honda old");*/
	
	/*honda_flux->Draw("hist same");
	min_honda_flux->Draw("hist same");	
	bnb_flux->Draw("hist same");
	old_honda_flux->Draw("hist same");*/	

	PlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/PeLEETuples_Atmospherics/FlatTTreePlots/"+CanvasName+".pdf");
	//delete PlotCanvas;

	//---------------------//

	// evaluate the spline

	double min_e = 0.;
	double max_e = 2.;
	double step = 0.05;
	int nbins = (max_e - min_e) / step;

	TH1D* spline = new TH1D("spline",";E_{#nu} [GeV];",nbins, min_e, max_e);
	//TH1D* spline = new TH1D("spline",";E_{#nu} [GeV];",NBinsEv,ArrayNBinEv);	

	for (int ibin = 0; ibin < nbins; ibin++){

		double e = min_e + ibin * step;

		int bnb_bin = bnb_flux->FindBin(e);
		int honda_bin = honda_flux->FindBin(e);
		int min_honda_bin = min_honda_flux->FindBin(e);
		int old_honda_bin = old_honda_flux->FindBin(e);				

		double bnb_entry = bnb_flux->GetBinContent(bnb_bin);
		double honda_entry = honda_flux->GetBinContent(honda_bin);	
		double min_honda_entry = min_honda_flux->GetBinContent(min_honda_bin);	
		double old_honda_entry = old_honda_flux->GetBinContent(old_honda_bin);			
		
		double weight = 0.5*(honda_entry + min_honda_entry) / bnb_entry;

		//cout << "weight = " << weight << endl;

		spline->SetBinContent(ibin+1,weight);

	}

	// make sure that we have non zero entries

	f_spline->cd();
	spline->Write();

	PrettyPlot(spline,T2KColor,false);

	spline->GetYaxis()->SetTitle("Honda/BNB");
	SplinePlotCanvas->cd();
	spline->Draw("hist");	

	SplinePlotCanvas->SaveAs("/exp/uboone/data/users/"+UserID+"/PeLEETuples_Atmospherics/FlatTTreePlots/"+SplineCanvasName+".pdf");
	//delete SplinePlotCanvas;	

	//---------------------//	

} // End of the program
