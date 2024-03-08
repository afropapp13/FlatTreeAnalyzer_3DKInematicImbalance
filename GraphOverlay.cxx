#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

using namespace std;

void GraphOverlay(TString Tag = "G18") {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.07;			
	int LineWidth = 5;

	TString OutFilePath = "OutputFiles/";

	//------------------------------//

	TString Names = "";

	// Event generators

	if (Tag == "G18") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }

	//if (Tag == "Gv2") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }
	if (Tag == "NEUT") { Names = OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"; }
	if (Tag == "NuWro") { Names = OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"; }
	if (Tag == "GiBUU") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021.root"; }
	//if (Tag == "G18 hN2018") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }
	//if (Tag == "G18 NoCoulomb") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }
	//if (Tag == "G18 RFG") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }
	//if (Tag == "G18 EffSF") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }
	//if (Tag == "G18 NoRPA") { Names = OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"; }

	TFile* File = new TFile(Names,"readonly");

	//------------------------------//

	TString CanvasName = "RejectionPower";
	TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	PlotCanvas->cd();
	PlotCanvas->SetTopMargin(0.12);
	PlotCanvas->SetLeftMargin(0.14);
	PlotCanvas->SetRightMargin(0.12);
	PlotCanvas->SetBottomMargin(0.15);		
	PlotCanvas->Draw();	

	TLegend* leg = new TLegend(0.7,0.65,0.8,0.85);
	leg->SetBorderSize(0);
	leg->SetNColumns(1);
	leg->SetTextSize(TextSize);	
	leg->SetTextFont(FontStyle);						
	leg->SetMargin(0.3);						
	
        TGraph* gPT = (TGraph*)(File->Get("gDeltaPT"));
        TGraph* gPn = (TGraph*)(File->Get("gDeltaPn"));

	gPT->SetTitle("");
	gPT->SetLineWidth(LineWidth);
	gPT->SetLineColor( kAzure-4 );	
	gPT->SetLineStyle( kSolid );	

	gPT->GetXaxis()->SetTitleFont(FontStyle);
	gPT->GetXaxis()->SetLabelFont(FontStyle);
	gPT->GetXaxis()->SetNdivisions(8);
	gPT->GetXaxis()->SetLabelSize(TextSize);
	gPT->GetXaxis()->SetTitleSize(TextSize);	
	gPT->GetXaxis()->SetTitleOffset(0.99);					
	gPT->GetXaxis()->CenterTitle();						
	gPT->GetXaxis()->SetTitle("signal acceptance");					

	gPT->GetYaxis()->SetTitleFont(FontStyle);
	gPT->GetYaxis()->SetLabelFont(FontStyle);
	gPT->GetYaxis()->SetNdivisions(6);
	gPT->GetYaxis()->SetLabelSize(TextSize);
	gPT->GetYaxis()->SetTitle("background rejection");
	gPT->GetYaxis()->SetTitleSize(TextSize);
	gPT->GetYaxis()->SetTitleOffset(1.);
	gPT->GetYaxis()->CenterTitle();	
	gPT->GetYaxis()->SetRangeUser(0.,1.15);	

	leg->AddEntry(gPT,Tag,"");
	TLegendEntry* legColorPT = leg->AddEntry(gPT,"#deltap_{T}","l");
	legColorPT->SetTextColor( kAzure-4 ); 
	
	gPT->Draw("AC");

	gPn->SetLineWidth(LineWidth);
	gPn->SetLineColor( kOrange+7 );	
	gPn->SetLineStyle( kSolid );	

	gPn->Draw("C");

	TLegendEntry* legColorPn = leg->AddEntry(gPn,"p_{n}","l");
	legColorPn->SetTextColor( kOrange+7 ); 

	leg->Draw();

	TLatex *textPanel = new TLatex();
	textPanel->SetTextFont(FontStyle);
	textPanel->SetTextSize(TextSize);
	//textPanel->DrawLatexNDC(0.17, 0.83, "(a)");

	PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+Tag+CanvasName+".pdf");
	delete PlotCanvas;

	//------------------------------//

	TString CutCanvasName = "CutPower";
	TCanvas* CutPlotCanvas = new TCanvas(CutCanvasName,CutCanvasName,205,34,1024,768);
	CutPlotCanvas->cd();
	CutPlotCanvas->SetTopMargin(0.12);
	CutPlotCanvas->SetLeftMargin(0.14);
	CutPlotCanvas->SetRightMargin(0.12);
	CutPlotCanvas->SetBottomMargin(0.15);		
	CutPlotCanvas->Draw();	

	TLegend* Cutleg = new TLegend(0.7,0.65,0.8,0.85);
	Cutleg->SetBorderSize(0);
	Cutleg->SetNColumns(1);
	Cutleg->SetTextSize(TextSize);	
	Cutleg->SetTextFont(FontStyle);						
	Cutleg->SetMargin(0.3);						
	
        TGraph* CutgPT = (TGraph*)(File->Get("gCutDeltaPT"));
        TGraph* CutgPn = (TGraph*)(File->Get("gCutDeltaPn"));

	CutgPT->SetTitle("");
	CutgPT->SetLineWidth(LineWidth);
	CutgPT->SetLineColor( kAzure-4 );	
	CutgPT->SetLineStyle( kSolid );	

	CutgPT->GetXaxis()->SetTitleFont(FontStyle);
	CutgPT->GetXaxis()->SetLabelFont(FontStyle);
	CutgPT->GetXaxis()->SetNdivisions(8);
	CutgPT->GetXaxis()->SetLabelSize(TextSize);
	CutgPT->GetXaxis()->SetTitleSize(TextSize);	
	CutgPT->GetXaxis()->SetTitleOffset(0.99);					
	CutgPT->GetXaxis()->CenterTitle();						
	CutgPT->GetXaxis()->SetTitle("minimum requirement [GeV/c]");					

	CutgPT->GetYaxis()->SetTitleFont(FontStyle);
	CutgPT->GetYaxis()->SetLabelFont(FontStyle);
	CutgPT->GetYaxis()->SetNdivisions(6);
	CutgPT->GetYaxis()->SetLabelSize(TextSize);
	CutgPT->GetYaxis()->SetTitle("rejection #times acceptance");
	CutgPT->GetYaxis()->SetTitleSize(TextSize);
	CutgPT->GetYaxis()->SetTitleOffset(1.);
	CutgPT->GetYaxis()->CenterTitle();	
	CutgPT->GetYaxis()->SetRangeUser(0.,0.69);	

	Cutleg->AddEntry(gPT,Tag,"");
	TLegendEntry* CutlegColorPT = Cutleg->AddEntry(CutgPT,"#deltap_{T}","l");
	CutlegColorPT->SetTextColor( kAzure-4 ); 
	
	CutgPT->Draw("AC");

	CutgPn->SetLineWidth(LineWidth);
	CutgPn->SetLineColor( kOrange+7 );	
	CutgPn->SetLineStyle( kSolid );	

	CutgPn->Draw("C");

	TLegendEntry* CutlegColorPn = Cutleg->AddEntry(CutgPn,"p_{n}","l");
	CutlegColorPn->SetTextColor( kOrange+7 ); 

	Cutleg->Draw();

	TLatex *textPanelb = new TLatex();
	textPanelb->SetTextFont(FontStyle);
	textPanelb->SetTextSize(TextSize);
	//textPanelb->DrawLatexNDC(0.17, 0.83, "(b)");

	CutPlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+Tag+CutCanvasName+".pdf");
	delete CutPlotCanvas;
	
	//------------------------------//

} // End of the program
