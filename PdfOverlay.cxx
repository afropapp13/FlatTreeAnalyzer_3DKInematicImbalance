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

void PdfOverlay(TString Tag = "") {

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

	vector<int> Colors;  vector<int> LineStyle;

	if (Tag == "") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18"); Colors.push_back(kBlack); LineStyle.push_back(kSolid);
	  //Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10.root"); Labels.push_back("Gv2"); Colors.push_back(kBlue); LineStyle.push_back(Gv2LineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10_MEC.root"); Labels.push_back("Gv2"); Colors.push_back(kBlue); LineStyle.push_back(Gv2LineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT"); Colors.push_back(kMagenta); LineStyle.push_back(NEUTLineStyle); // kMagenta - 9
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"); Labels.push_back("NuWro"); Colors.push_back(NEUTColor); LineStyle.push_back(NuWroLineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021.root"); Labels.push_back("GiBUU"); Colors.push_back(GiBUUColor); LineStyle.push_back(GiBUULineStyle);

	}

	if (Tag == "GENIE") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18");Colors.push_back(kBlack); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_hN2018.root"); Labels.push_back("hN"); Colors.push_back(kBlue); LineStyle.push_back(Gv2LineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoRPA.root"); Labels.push_back("NoRPA"); Colors.push_back(kMagenta-9); LineStyle.push_back(NEUTLineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoCoulomb.root"); Labels.push_back("NoC"); Colors.push_back(NEUTColor); LineStyle.push_back(NuWroLineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_RFG.root"); Labels.push_back("RFG"); Colors.push_back(GiBUUColor); LineStyle.push_back(GiBUULineStyle);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_EffSF.root"); Labels.push_back("SF");Colors.push_back(kRed); LineStyle.push_back(kSolid);

	}

	if (Tag == "Closure") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18 v3.0.6");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("G18 v3.4.0");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);

	}

	if (Tag == "Weights") {

	  //Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18a");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("G18");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11b.root"); Labels.push_back("G18b");Colors.push_back(kMagenta-9); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_03_320.root"); Labels.push_back("03_320"); Colors.push_back(kGreen+2); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_03_330.root"); Labels.push_back("03_330"); Colors.push_back(kBlue); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_RS.root"); Labels.push_back("RS");Colors.push_back(kRed-1); LineStyle.push_back(kSolid);

	}

	if (Tag == "WeightsClosure") {

	  // Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("T2K G18 v3.0.6");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"WeightsFlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("T2K G18 v3.4.0");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);

	}

	if (Tag == "DUNE") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_AR23_20i_00_000.root"); Labels.push_back("AR23_20i_00_000");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("G18_10a_02_11a");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);

	}

	if (Tag == "EnergyIndependence") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("uB Flux");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_0_3GeV.root"); Labels.push_back("0.3 GeV");Colors.push_back(kBlue); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_0_5GeV.root"); Labels.push_back("0.5 GeV");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_0_8GeV.root"); Labels.push_back("0.8 GeV");Colors.push_back(kGreen+1); LineStyle.push_back(kSolid);

	  //Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_1GeV.root"); Labels.push_back("G18 1 GeV");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  //Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a_3GeV.root"); Labels.push_back("G18 3 GeV");Colors.push_back(kGreen+1); LineStyle.push_back(kSolid);

	}

	if (Tag == "uBDUNEFDND") {

	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_G18_DUNE_FD.root"); Labels.push_back("DUNE FD");Colors.push_back(kGreen+1); LineStyle.push_back(kSolid);	
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_G18_DUNE_ND.root"); Labels.push_back("DUNE ND");Colors.push_back(kOrange+7); LineStyle.push_back(kSolid);
	  Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_4_0_G18_10a_02_11a.root"); Labels.push_back("MicroBooNE");Colors.push_back(kAzure-4); LineStyle.push_back(kSolid);

	}


	const int NSamples = Names.size();
	const int NColors = Colors.size();

	// Sanity check
	if (NColors < NSamples) { cout << "Give me some more colors!" << endl; return; }

	std::vector<TFile*> Files; Files.resize(NSamples);

	//------------------------------//

	// Plots to overlay

	std::vector<TString> PlotNames;
	std::vector<TString> YAxisLabel;

	//------------------------------//

	// Post FSI

	// 1D

	PlotNames.push_back("DISTrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("RESTrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("MECTrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("QETrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("TrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DMuPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D,#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("RESTrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}}  #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp ,x}}  #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp ,y}}  #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}}  #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");

	//------------------------------//

	// Pre FSI

	// 1D
	/*
	PlotNames.push_back("NoFSITrueFineBinMuonCosThetaPlot"); YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlpha3DMuPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D,#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("NoFSITrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}} [10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
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

		TString CanvasName = "ThreeDKI_PdfOverlay_" + PlotNames[iPlot];
		TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetTopMargin(0.15);
		PlotCanvas->SetLeftMargin(0.17);
		PlotCanvas->SetRightMargin(0.05);
		PlotCanvas->SetBottomMargin(0.16);		
		PlotCanvas->Draw();	

		TLegend* leg = new TLegend(0.17,0.855,0.99,0.985);
		leg->SetBorderSize(0);
		leg->SetNColumns(6);
		leg->SetTextSize(TextSize);	
		leg->SetTextFont(FontStyle);						
		leg->SetMargin(0.3);						

		// Loop over the samples to open the files and to get the corresponding plot

		std::vector<TH1D*> Histos; Histos.resize(NSamples);

		for (int iSample = 0; iSample < NSamples; iSample++) {	

		        Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]));

			Histos[iSample]->SetLineWidth(3);
			Histos[iSample]->SetLineColor( Colors.at(iSample) );	
			Histos[iSample]->SetLineStyle( LineStyle.at(iSample) );	

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
			//Histos[iSample]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
			Histos[iSample]->GetYaxis()->SetTitle("PDF");
			Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
			Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
			//Histos[iSample]->GetYaxis()->SetTickSize(0);
			Histos[iSample]->GetYaxis()->CenterTitle();	
			Histos[iSample]->Scale(1./Histos[iSample]->Integral());	

			double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());
			double MaxScale = 1.05;
			Histos[iSample]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);
			Histos[0]->GetYaxis()->SetRangeUser(0.,MaxScale*imax);			

			PlotCanvas->cd();
			Histos[iSample]->Draw("hist same");
			//Histos[0]->Draw("hist same");	

			TLegendEntry* legColor = leg->AddEntry(Histos[iSample],Labels[iSample],"l");
			legColor->SetTextColor( Colors.at(iSample) ); 

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(TextSize);
			TString PlotNameDuplicate = PlotNames[iPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("TrueFineBin","") ;
			textSlice->DrawLatexNDC(0.2, 0.8, LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

			if ( 
			    PlotNames[iPlot] == "TrueFineBinDeltaPnPlot" ||
			    PlotNames[iPlot] == "TrueFineBinDeltaAlpha3DqPlot" || 
			    PlotNames[iPlot] == "TrueFineBinDeltaPhi3DPlot"
			     ) {

			  TLatex *textPanel = new TLatex();
			  textPanel->SetTextFont(FontStyle);
			  textPanel->SetTextSize(TextSize);

			  if ( PlotNames[iPlot] == "TrueFineBinDeltaPnPlot") { textPanel->DrawLatexNDC(0.87, 0.8, "(a)"); }
			  if ( PlotNames[iPlot] == "TrueFineBinDeltaAlpha3DqPlot") { textPanel->DrawLatexNDC(0.87, 0.8, "(b)"); }
			  if ( PlotNames[iPlot] == "TrueFineBinDeltaPhi3DPlot") { textPanel->DrawLatexNDC(0.87, 0.8, "(c)"); }
			
			}

			//----------------------------------------//					

		} // End of the loop over the samples grabing the plots	

		PlotCanvas->cd();
		leg->Draw();

		PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+Tag+CanvasName+".pdf");
		delete PlotCanvas;

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
