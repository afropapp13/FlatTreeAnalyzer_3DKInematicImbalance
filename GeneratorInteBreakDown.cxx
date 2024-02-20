#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
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

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

void GeneratorInteBreakDown() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			

	TString OutFilePath = "OutputFiles/";

	//------------------------------//

	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; 

	std::vector<int> Colors{kBlack,kBlue,kRed+1,kOrange+7,kGreen+1, kMagenta+1};
	std::vector<TString> Process{"","QE","MEC","RES","DIS"};

	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10_MEC.root"); Labels.push_back("GENIE v2.12.10.MEC");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT 5.4.0.1");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"); Labels.push_back("NuWro 19.02.1");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021.root"); Labels.push_back("GiBUU");

	const int NSamples = Names.size();
	const int NColors = Colors.size();
	const int NProcesses = Process.size();

	// Sanity check
	if (NColors < NProcesses) { cout << "Give me some more colors for all the processes!" << endl; return; }

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;
	std::vector<TString> YAxisLabel;

	//------------------------------//

	// Post FSI

	// 1D

	PlotNames.push_back("InclusiveMuonCosThetaSingleBinPlot"); YAxisLabel.push_back("#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");	
	PlotNames.push_back("InclusiveMuonCosThetaPlot"); YAxisLabel.push_back("#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("MuonCosThetaSingleBinPlot"); YAxisLabel.push_back("#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("InclusiveMuonEnergyPlot"); YAxisLabel.push_back("#sigma #left[10^{-38} #frac{cm^{2}}{Ar GeV}#right]");
	PlotNames.push_back("InclusiveNeutronMultiPlot"); YAxisLabel.push_back("#sigma #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	
	//------------------------------//

	const int NPlots = PlotNames.size();
	const int NLabels = YAxisLabel.size();

	// sanity check
	if ( NPlots != NLabels) { cout << "Inconsistent number of plots and labels! Aborting !" << endl; return; }

	//------------------------------//	

	// Loop over the samples to open the files and the TTree

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples

	//------------------------------//

	// Loop over the plots to be compared

	for (int iPlot = 0; iPlot < NPlots; iPlot++) {

		for (int iSample = 0; iSample < NSamples; iSample++) {	

		  TString LabelCopy = Labels[iSample];
		  TString CanvasName = "ThreeDKI_"+LabelCopy.ReplaceAll(" ","_")+"_InteBreakDown_" + PlotNames[iPlot];
		  TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		  PlotCanvas->cd();
		  PlotCanvas->SetTopMargin(0.13);
		  PlotCanvas->SetLeftMargin(0.17);
		  PlotCanvas->SetRightMargin(0.05);
		  PlotCanvas->SetBottomMargin(0.16);		
		  PlotCanvas->Draw();	

		  TLegend* leg = new TLegend(0.01,0.88,0.99,0.99);
		  leg->SetBorderSize(0);
		  leg->SetNColumns(5);
		  leg->SetTextSize(TextSize-0.01);	
		  leg->SetTextFont(FontStyle);						
		  leg->SetMargin(0.2);				
		  leg->SetFillColor(0);				

		  // Loop over the interaction processes

		  std::vector<TH1D*> Histos; Histos.resize(NProcesses);

		  for (int iprocess = 0; iprocess < NProcesses; iprocess++) {	

		        Histos[iprocess] = (TH1D*)(Files[iSample]->Get(Process[iprocess]+PlotNames[iPlot]));

			Histos[iprocess]->SetLineWidth(4);
			Histos[iprocess]->SetLineColor( Colors.at(iprocess) );	

			Histos[iprocess]->GetXaxis()->SetTitleFont(FontStyle);
			Histos[iprocess]->GetXaxis()->SetLabelFont(FontStyle);
			Histos[iprocess]->GetXaxis()->SetNdivisions(8);
			Histos[iprocess]->GetXaxis()->SetLabelSize(TextSize);
			Histos[iprocess]->GetXaxis()->SetTitleSize(TextSize);	
			Histos[iprocess]->GetXaxis()->SetTitleOffset(1.1);					
			Histos[iprocess]->GetXaxis()->CenterTitle();						

			Histos[iprocess]->GetYaxis()->SetTitleFont(FontStyle);
			Histos[iprocess]->GetYaxis()->SetLabelFont(FontStyle);
			Histos[iprocess]->GetYaxis()->SetNdivisions(6);
			Histos[iprocess]->GetYaxis()->SetLabelSize(TextSize);
			Histos[iprocess]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
			Histos[iprocess]->GetYaxis()->SetTitleSize(TextSize);
			Histos[iprocess]->GetYaxis()->SetTitleOffset(1.3);
			//Histos[iprocess]->GetYaxis()->SetTickSize(0);
			Histos[iprocess]->GetYaxis()->CenterTitle();	
			Histos[iprocess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[0]->GetMaximum());

			Histos[iprocess]->Draw("hist same");
			Histos[0]->Draw("hist same");	

			TString LegLabel = Process[iprocess];
			if (LegLabel == "") { LegLabel = "Total"; }
			int frac = Histos[iprocess]->Integral()/Histos[0]->Integral()*100.;
			TLegendEntry* legColor = leg->AddEntry(Histos[iprocess],LegLabel + " (" + TString(std::to_string(frac)) + "%)","l");
			legColor->SetTextColor( Colors.at(iprocess) ); 

		  } // End of the loop over the processes



		  PlotCanvas->cd();
		  leg->Draw();

		  TLatex *textSlice = new TLatex();
		  textSlice->SetTextFont(FontStyle);
		  textSlice->SetTextSize(TextSize);
		  TString PlotNameDuplicate = PlotNames[iPlot];
		  TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("TrueFineBin","") ;
		  textSlice->DrawLatexNDC(0.2, 0.81, Labels[iSample] + "      " + LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

		  gPad->RedrawAxis();

		  PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+CanvasName+".pdf");
		  delete PlotCanvas;

		} // End of the loop over the samples grabing the plots	

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
