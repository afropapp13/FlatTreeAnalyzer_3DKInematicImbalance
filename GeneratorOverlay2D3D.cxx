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

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void GeneratorOverlay2D3D(TString Tag = "") {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			

	TString OutFilePath = "OutputFiles/";

	static std::map<TString,TString> Mapping = {

	  { "TrueFineBinDeltaPTPlot", "TrueFineBinDeltaPnPlot" },
	  { "TrueFineBinDeltaAlphaTPlot", "TrueFineBinDeltaAlpha3DqPlot" },
	  { "TrueFineBinDeltaPhiTPlot", "TrueFineBinDeltaPhi3DPlot" }
	};

	static std::map<TString,TString> TKILabels = {

	  { "TrueFineBinDeltaPTPlot", "#delta p_{T}" },
	  { "TrueFineBinDeltaAlphaTPlot", "#delta#alpha_{T}" },
	  { "TrueFineBinDeltaPhiTPlot", "#delta#phi_{T}" }
	};

	static std::map<TString,TString> GKILabels = {

	  { "TrueFineBinDeltaPnPlot","p_{n}" },
	  { "TrueFineBinDeltaAlpha3DqPlot", "#alpha_{3D}" },
	  { "TrueFineBinDeltaPhi3DPlot", "#phi_{3D}" }
	};



	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; 

	vector<int> Colors; 

	if (Tag == "") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18"); Colors.push_back(kBlack);
	  //Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT"); Colors.push_back(kMagenta-9);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT"); Colors.push_back(kYellow-6);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"); Labels.push_back("NuWro"); Colors.push_back(NEUTColor);
	
	}


	const int NSamples = Names.size();
	const int NColors = Colors.size();

	// Sanity check
	if (NColors < NSamples) { cout << "Give me some more colors!" << endl; return; }

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;
	vector<TString> PanelLabel;
	vector<TString> XaxisLabel;

	//------------------------------//

	// Post FSI

	// 1D

	//PlotNames.push_back("TrueFineBinDeltaPTPlot"); PanelLabel.push_back("(a)"); XaxisLabel.push_back("Missing Momentum [GeV/c]");

	//PlotNames.push_back("TrueFineBinDeltaAlphaTPlot"); PanelLabel.push_back("(a)"); XaxisLabel.push_back("Struck Nucleon - Missing Momentum Angle [deg]");
	//PlotNames.push_back("TrueFineBinDeltaPhiTPlot"); PanelLabel.push_back("(b)"); XaxisLabel.push_back("Proton - Missing Momentum Angle [deg]");

	PlotNames.push_back("TrueFineBinDeltaAlphaTPlot"); PanelLabel.push_back(""); XaxisLabel.push_back("Struck Nucleon - Missing Momentum Angle [deg]");
	PlotNames.push_back("TrueFineBinDeltaPhiTPlot"); PanelLabel.push_back(""); XaxisLabel.push_back("Proton - Missing Momentum Angle [deg]");


	//------------------------------//

	const int NPlots = PlotNames.size();

	//------------------------------//	

	// Loop over the samples to open the files and the TTree

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples

	//------------------------------//

	// Loop over the plots to be compared

	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		TString CanvasName = "GeneratorOverlay2D3D_" + PlotNames[iPlot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.16);
		PlotCanvas->SetLeftMargin(0.1);
		PlotCanvas->SetRightMargin(0.05);
		PlotCanvas->SetBottomMargin(0.16);		
		PlotCanvas->Draw();	

		TLegend* leg = new TLegend(0.15,0.925,0.99,0.99);
		leg->SetBorderSize(0);
		leg->SetNColumns(3);
		leg->SetTextSize(TextSize);	
		leg->SetTextFont(FontStyle);						
		leg->SetMargin(0.1);						
		leg->SetFillStyle(0);

		TLegend* legTKI = new TLegend(0.15,0.855,0.99,0.925);
		legTKI->SetBorderSize(0);
		legTKI->SetNColumns(3);
		legTKI->SetTextSize(TextSize);	
		legTKI->SetTextFont(FontStyle);						
		legTKI->SetMargin(0.1);	
		legTKI->SetFillStyle(0);

		// Loop over the samples to open the files and to get the corresponding plot

		std::vector<TH1D*> Histos; Histos.resize(NSamples);
		std::vector<TH1D*> NoFSIHistos; NoFSIHistos.resize(NSamples);

		std::vector<TH1D*> GKIHistos; GKIHistos.resize(NSamples);
		std::vector<TH1D*> NoFSIGKIHistos; NoFSIGKIHistos.resize(NSamples);

		for (int iSample = 0; iSample < NSamples; iSample++) {	

		        Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]));
		        GKIHistos[iSample] = (TH1D*)(Files[iSample]->Get(Mapping[ PlotNames[iPlot] ] ));


			// Special case, GiBUU needs a separate production with no fsi                                                                                                                 
			// otherwise the _vert branches make no sense                                                                                                                                   

			TFile* file = Files[iSample];

			if (Labels[iSample] == "GiBUU") {

			  file = new TFile(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021_NoFSI.root","readonly");


			}

		        NoFSIHistos[iSample] = (TH1D*)(file->Get("NoFSI" + PlotNames[iPlot]));
		        NoFSIGKIHistos[iSample] = (TH1D*)(file->Get( "NoFSI" + Mapping[ PlotNames[iPlot] ] ));


			Histos[iSample]->Divide(NoFSIHistos[iSample]);
			GKIHistos[iSample]->Divide(NoFSIGKIHistos[iSample]);

			Histos[iSample]->SetLineWidth(3);
			Histos[iSample]->SetLineColor( Colors.at(iSample) );	
			Histos[iSample]->SetLineStyle( kDashed );	

			GKIHistos[iSample]->SetLineWidth(3);
			GKIHistos[iSample]->SetLineColor( Colors.at(iSample) );	
			GKIHistos[iSample]->SetLineStyle( kSolid );	

			Histos[iSample]->GetXaxis()->SetTitle(XaxisLabel[iPlot]);
			Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
			Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
			Histos[iSample]->GetXaxis()->SetNdivisions(8);
			Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
			Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
			Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
			Histos[iSample]->GetXaxis()->CenterTitle();						

			Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
			Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
			Histos[iSample]->GetYaxis()->SetNdivisions(6);
			Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitle("Ratio FSI / no-FSI");
			Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitleOffset(0.8);
			//Histos[iSample]->GetYaxis()->SetTickSize(0);
			Histos[iSample]->GetYaxis()->CenterTitle();	

			double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());
			imax = TMath::Max(imax,GKIHistos[0]->GetMaximum());
			imax = TMath::Max(imax,GKIHistos[iSample]->GetMaximum());

			double MaxScale = 1.05;
			Histos[iSample]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);
			Histos[0]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);			

			PlotCanvas->cd();
			Histos[iSample]->Draw("hist same");
			GKIHistos[iSample]->Draw("hist same");
			Histos[0]->Draw("hist same");	

			TLegendEntry* legColor = legTKI->AddEntry(Histos[iSample],Labels[iSample] + " " + TKILabels[PlotNames[iPlot]],"l");
			legColor->SetTextColor( Colors.at(iSample) ); 

			TLegendEntry* legColorGKI = leg->AddEntry(GKIHistos[iSample],Labels[iSample] + " " + GKILabels[Mapping[PlotNames[iPlot]]],"l");
			legColorGKI->SetTextColor( Colors.at(iSample) ); 

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(TextSize);
			TString PlotNameDuplicate = PlotNames[iPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("TrueFineBin","") ;
			textSlice->DrawLatexNDC(0.2, 0.8, LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

			TLatex *textPanel = new TLatex();
			textPanel->SetTextFont(FontStyle);
			textPanel->SetTextSize(TextSize);

			textPanel->DrawLatexNDC(0.87, 0.78, PanelLabel[iPlot]);
			
			//----------------------------------------//					

		} // End of the loop over the samples grabing the plots	

		PlotCanvas->cd();
		leg->Draw("same");
		legTKI->Draw("same");

		PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+Tag+CanvasName+".pdf");
		delete PlotCanvas;

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
