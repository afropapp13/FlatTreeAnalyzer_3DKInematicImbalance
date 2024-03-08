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

void GeneratorFSIStudy() {

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

	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6.root"); Labels.push_back("G18");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10.root"); Labels.push_back("GENIE v2.12.10");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v2_12_10_MEC.root"); Labels.push_back("GENIE v2.12.10.MEC");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT_5_4_0_1.root"); Labels.push_back("NEUT");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro_19_02_1.root"); Labels.push_back("NuWro");
	Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021.root"); Labels.push_back("GiBUU");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021_NoFSI.root"); Labels.push_back("GiBUU 2021 NoFSI");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_hN2018.root"); Labels.push_back("v3.0.6 hN2018");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoRPA.root"); Labels.push_back("v3.0.6 NoRPA");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_NoCoulomb.root"); Labels.push_back("v3.0.6 NoCoulomb");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_RFG.root"); Labels.push_back("v3.0.6 RFG");
	//Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_v3_0_6_EffSF.root"); Labels.push_back("v3.0.6 EffSF");

	const int NSamples = Names.size();

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
	PlotNames.push_back("TrueFineBinDeltaPnPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DMuPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D,#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlpha3DqPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaAlphaTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhiTPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPhi3DPlot"); YAxisLabel.push_back("#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpxPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,x}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnPerpyPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#perp,y}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");
	PlotNames.push_back("TrueFineBinDeltaPnParPlot"); YAxisLabel.push_back("#frac{d#sigma}{dp_{n#parallel}} #left[10^{-38} #frac{cm^{2}}{(GeV/c) Ar}#right]");

	//------------------------------//

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
		  TString CanvasName = "ThreeDKI_"+LabelCopy.ReplaceAll(" ","_")+"_FSIStudy_" + PlotNames[iPlot];
		  TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
		  PlotCanvas->cd();
		  PlotCanvas->SetTopMargin(0.13);
		  PlotCanvas->SetLeftMargin(0.19);
		  PlotCanvas->SetRightMargin(0.05);
		  PlotCanvas->SetBottomMargin(0.16);		
		  PlotCanvas->Draw();	

		  TLegend* leg = new TLegend(0.2,0.88,0.6,0.99);
		  leg->SetBorderSize(0);
		  leg->SetNColumns(3);
		  leg->SetTextSize(TextSize);	
		  leg->SetTextFont(FontStyle);						
		  leg->SetMargin(0.2);				

		  TH1D* Histos = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]));
		  Histos->SetLineWidth(4);
		  Histos->SetLineColor(kBlack);	

		  Histos->GetXaxis()->SetTitleFont(FontStyle);
		  Histos->GetXaxis()->SetLabelFont(FontStyle);
		  Histos->GetXaxis()->SetNdivisions(8);
		  Histos->GetXaxis()->SetLabelSize(TextSize);
		  Histos->GetXaxis()->SetTitleSize(TextSize);	
		  Histos->GetXaxis()->SetTitleOffset(1.1);					
		  Histos->GetXaxis()->CenterTitle();						

		  Histos->GetYaxis()->SetTitleFont(FontStyle);
		  Histos->GetYaxis()->SetLabelFont(FontStyle);
		  Histos->GetYaxis()->SetNdivisions(6);
		  Histos->GetYaxis()->SetLabelSize(TextSize);
		  Histos->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
		  Histos->GetYaxis()->SetTitleSize(TextSize);
		  Histos->GetYaxis()->SetTitleOffset(1.4);
		  //Histos->GetYaxis()->SetTickSize(0);
		  Histos->GetYaxis()->CenterTitle();	

		  Histos->Draw("hist same");
		  TLegendEntry* legColor = leg->AddEntry(Histos,"FSI","l");
		  legColor->SetTextColor(kBlack); 

		  //-------------------//

		  // Special case, GiBUU needs a separate production with no fsi
		  // otherwise the _vert branches make no sense

		  TFile* file = Files[iSample];

		  if (Labels[iSample] == "GiBUU 2021") {

		    file = new TFile(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_2021_NoFSI.root","readonly");

		  }

		  TH1D* NoFSIHistos = (TH1D*)(file->Get("NoFSI"+PlotNames[iPlot]));

		  NoFSIHistos->SetLineWidth(4);
		  NoFSIHistos->SetLineColor(kBlack);	
		  NoFSIHistos->SetLineStyle(kDashed);	

		  NoFSIHistos->Draw("hist same");
		  TLegendEntry* NoFSIlegColor = leg->AddEntry(NoFSIHistos,"No FSI","l");
		  NoFSIlegColor->SetTextColor(kBlack); 

		  //-------------------//

		  double max = TMath::Max(Histos->GetMaximum(),NoFSIHistos->GetMaximum());
		  Histos->GetYaxis()->SetRangeUser(0.,1.15*max);
		  Histos->Draw("same hist");

		  //-------------------//

		  PlotCanvas->cd();
		  leg->Draw();

		  TLatex *textSlice = new TLatex();
		  textSlice->SetTextFont(FontStyle);
		  textSlice->SetTextSize(TextSize);
		  TString PlotNameDuplicate = PlotNames[iPlot];
		  TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("TrueFineBin","") ;
		  textSlice->DrawLatexNDC(0.25, 0.82, LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

                  TLatex *textLabel = new TLatex();
                  textLabel->SetTextFont(FontStyle);
                  textLabel->SetTextSize(TextSize);
                  textLabel->DrawLatexNDC(0.8, 0.91, Labels[iSample]);

                  if ( 
		      PlotNames[iPlot] == "TrueFineBinDeltaAlpha3DqPlot" || PlotNames[iPlot] == "TrueFineBinDeltaAlphaTPlot" ||
                      PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" ||
                      PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" ||
                      PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_00To0_20Plot" ||
                      PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
		       ) {

                    TLatex *textPanel = new TLatex();
                    textPanel->SetTextFont(FontStyle);
                    textPanel->SetTextSize(TextSize);
                    if ( 
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3DqPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_00To0_20Plot"
			 ) { 
			
				//textPanel->DrawLatexNDC(0.87, 0.81, "(a)"); 

			}

                    if ( 
			PlotNames[iPlot] == "TrueFineBinDeltaAlphaTPlot" || 
			PlotNames[iPlot] == "TrueFineBinDeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" ||
			PlotNames[iPlot] == "TrueFineBinDeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"

			 ) { 

			//	textPanel->DrawLatexNDC(0.87, 0.81, "(b)"); 

			}

                  }


		  PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/"+CanvasName+".pdf");
		  delete PlotCanvas;

		} // End of the loop over the samples grabing the plots	

	} // End of the loop over the plots

	//------------------------------//

} // End of the program
