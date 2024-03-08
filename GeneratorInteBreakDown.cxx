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

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"

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

	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18");
	Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("G18T");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10.root"); Labels.push_back("Gv2");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10_MEC.root"); Labels.push_back("GENIE v2.12.10.MEC");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT 5.4.0.1");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"); Labels.push_back("NuWro 19.02.1");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021.root"); Labels.push_back("GiBUU");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021_NoFSI.root"); Labels.push_back("GiBUU no-FSI");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_hN2018.root"); Labels.push_back("v3.0.6 hN2018");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoRPA.root"); Labels.push_back("v3.0.6 NoRPA");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoCoulomb.root"); Labels.push_back("v3.0.6 NoCoulomb");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_RFG.root"); Labels.push_back("v3.0.6 RFG");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_EffSF.root"); Labels.push_back("v3.0.6 EffSF");

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

	PlotNames.push_back("TrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DMuPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D,#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPtxPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPtxPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPtyPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPtyPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");

	//------------------------------//

	// Pre FSI

	// 1D
	/*
	PlotNames.push_back("NoFSITrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} [10^{-38} #frac{cm^{2}}{Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlpha3DMuPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D,#mu}} [10^{-38} #frac{cm^{2}}{deg Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} [10^{-38} #frac{cm^{2}}{deg Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} [10^{-38} #frac{cm^{2}}{deg Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} [10^{-38} #frac{cm^{2}}{deg Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} [10^{-38} #frac{cm^{2}}{deg Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,x}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,y}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}]");
	*/

	// 2D Post FSI

	PlotNames.push_back("TrueFineBinDeltaAlphaT_DeltaPT_0_00To0_20Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlphaT_DeltaPT_0_20To0_40Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlphaT_DeltaPT_0_40To1_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

	PlotNames.push_back("TrueFineBinDeltaAlpha3Dq_DeltaPn_0_00To0_20Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3Dq_DeltaPn_0_20To0_40Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

	PlotNames.push_back("TrueFineBinDeltaAlpha3DMu_DeltaPn_0_00To0_20Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DMu_DeltaPn_0_20To0_40Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DMu_DeltaPn_0_40To1_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

	PlotNames.push_back("TrueFineBinDeltaPT_DeltaAlphaT_0_00To45_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPT_DeltaAlphaT_45_00To90_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPT_DeltaAlphaT_90_00To135_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPT_DeltaAlphaT_135_00To180_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3Dq_0_00To45_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3Dq_45_00To90_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3Dq_90_00To135_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3Dq_135_00To180_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3DMu_0_00To45_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3DMu_45_00To90_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3DMu_90_00To135_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPn_DeltaAlpha3DMu_135_00To180_00Plot"); YAxisLabel.push_back("#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]");

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

		  TLegend* leg = new TLegend(0.2,0.88,0.8,0.99);
		  leg->SetBorderSize(0);
		  leg->SetNColumns(5);
		  leg->SetTextSize(TextSize);	
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
			Histos[iprocess]->GetYaxis()->SetTitleOffset(1.25);
			//Histos[iprocess]->GetYaxis()->SetTickSize(0);
			Histos[iprocess]->GetYaxis()->CenterTitle();	
			Histos[iprocess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[0]->GetMaximum());

			Histos[iprocess]->Draw("hist same");
			Histos[0]->Draw("hist same");	

			TString LegLabel = Process[iprocess];
			if (LegLabel == "") { LegLabel = "Total"; }
			TLegendEntry* legColor = leg->AddEntry(Histos[iprocess],LegLabel,"l");
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

		  //TLatex *textLabel = new TLatex();
		  //textLabel->SetTextFont(FontStyle);
		  //textLabel->SetTextSize(TextSize);
		  //textLabel->DrawLatexNDC(0.2, 0.81, Labels[iSample]);

		  if ( 
		      PlotNames[iPlot] == "TrueFineBinDeltaPnPlot" || 
		      PlotNames[iPlot] == "TrueFineBinDeltaPTPlot"  || 
		      PlotNames[iPlot] == "TrueFineBinDeltaPhi3DPlot"  || 
		      PlotNames[iPlot] == "TrueFineBinDeltaPhiTPlot" ||
		      PlotNames[iPlot] == "TrueFineBinDeltaPnParPlot" ||
		      PlotNames[iPlot] == "NoFSITrueFineBinDeltaPnParPlot" ||
		      PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" ||
		      PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" ||
		      PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_00To0_20Plot" ||
		      PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
		       ) {

		    TLatex *textPanel = new TLatex();
		    textPanel->SetTextFont(FontStyle);
		    textPanel->SetTextSize(TextSize);

		    if ( 
			PlotNames[iPlot] == "TrueFineBinDeltaPnPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3DqPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaPhi3DPlot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaPnParPlot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaPnPerpPlot" ||
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPnPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaAlpha3DqPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPhi3DPlot" ||
			//PlotNames[iPlot] == "NoFSITrueFineBinDeltaPnParPlot" ||
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPnPerpPlot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_00To0_20Plot"
			 ) { 

			//	textPanel->DrawLatexNDC(0.87, 0.8, "(a)"); 

			}

		    if ( 
			PlotNames[iPlot] == "TrueFineBinDeltaPTPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaAlphaTPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaPhiTPlot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaPtxPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaPtyPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPTPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaAlphaTPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPhiTPlot" ||
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPtxPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPtyPlot" || 
			PlotNames[iPlot] == "NoFSITrueFineBinDeltaPnParPlot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
			 ) { 

				//textPanel->DrawLatexNDC(0.87, 0.8, "(b)"); 

			}

		  }


		  gPad->RedrawAxis();

		  PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+CanvasName+".pdf");
		  delete PlotCanvas;

		} // End of the loop over the samples grabing the plots	

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
