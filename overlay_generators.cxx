#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../myClasses/Constants.h"
#include "../myClasses/myFunctions.cpp"

using namespace std;
using namespace Constants;

//----------------------------------------//

void overlay_generators() {

	//----------------------------------------//

	TH1D::SetDefaultSumw2();
	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.07,"t");
	gStyle->SetOptStat(0);
	
	//----------------------------------------//

	vector<TString> PlotNames;
	//PlotNames.push_back("ThetaVisPlot");
	PlotNames.push_back("MuonCosThetaPlot");
	//PlotNames.push_back("MuonCosThetaSingleBinPlot");
	//PlotNames.push_back("SerialThetaVis_ECalPlot");
	//PlotNames.push_back("SerialThetaVis_DeltaPnPlot");

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
	//Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();

		NameOfSamples.push_back("Overlay9"); Colors.push_back(kAzure-4); Labels.push_back("G18T ");
		NameOfSamples.push_back("GENIE_v3_0_6"); Colors.push_back(kBlack); Labels.push_back("G18 ");
		NameOfSamples.push_back("NoTuneOverlay9"); Colors.push_back(kOrange+7); Labels.push_back("G18D ");

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToExtractedXSec+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TString TrueString = "NoSmearTrue";

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

				}

			}

			else if (NameOfSamples[WhichSample] == "NoTuneOverlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToExtractedXSec+"/NoTuneOverlay9WienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TString TrueString = "NoSmearAltTrue";

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					rm_bin_width(histTrue);
					CurrentPlotsTrue.push_back(histTrue);

				}

			}

			else if (NameOfSamples[WhichSample] == "Overlay9NuWro") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToExtractedXSec+"/Overlay9NuWroWienerSVD_ExtractedXSec_Overlay9_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TString TrueString = "NoSmearAltTrue";

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					rm_bin_width(histTrue);
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}


			else {
		
			  FileSample.push_back(TFile::Open("OutputFiles/FlatTreeAnalyzerOutput_"+NameOfSamples[WhichSample]+".root")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TString TrueString = "True";
					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);
		
				}

			}

			PlotsTrue.push_back(CurrentPlotsTrue);

		}

		//----------------------------------------//

		// Loop over the plots

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			// -----------------------------------------------------------------------------------------------------------------------------			

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.16);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.19);
			midPad->SetRightMargin(0.03);			
			midPad->Draw();

			TLegend* leg = new TLegend(0.39,0.69,0.72,0.85);
			TLegend* legMC = new TLegend(0.3,0.69,0.4,0.85);
			
			leg->SetBorderSize(0);
			leg->SetTextSize(0.05);
			leg->SetTextFont(FontStyle);
			leg->SetNColumns(1);
			leg->SetMargin(0.15);

			legMC->SetBorderSize(0);
			legMC->SetTextSize(0.05);
			legMC->SetTextFont(FontStyle);
			legMC->SetNColumns(1);
			legMC->SetMargin(0.3);				

			midPad->cd();
			
			//------------------------------//

			// The N-dimensional analysis has been developed based on the bin number, not the actual range

			if (string(PlotNames[WhichPlot]).find("Serial") != std::string::npos) {	

				TString XaxisTitle = PlotsTrue[0][WhichPlot]->GetXaxis()->GetTitle();
				XaxisTitle.ReplaceAll("deg","bin #");
				XaxisTitle.ReplaceAll("GeV/c","bin #");
				XaxisTitle.ReplaceAll("GeV","bin #");				
				PlotsTrue[0][WhichPlot]->GetXaxis()->SetTitle(XaxisTitle);

				TString YaxisTitle = VarLabel[PlotNames[WhichPlot]];
				YaxisTitle.ReplaceAll("deg","");
				YaxisTitle.ReplaceAll("GeV/c","");
				YaxisTitle.ReplaceAll("GeV","");
				YaxisTitle.ReplaceAll("/c","");
				PlotsTrue[0][WhichPlot]->GetYaxis()->SetTitle(YaxisTitle);				

			}			

			//------------------------------//

			// Overlay GENIE v3 + Tune

			PrettyPlot(PlotsTrue[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetLineWidth(3);
			PlotsTrue[0][WhichPlot]->Draw("hist same");

			TLegendEntry* lOv = legMC->AddEntry(PlotsTrue[0][WhichPlot],Labels[0],"l");
			lOv->SetTextColor(Colors[0]); 										

			// Clones for the NSamples-1 model predictions
			// index 0 corresponds to nominal overlay / CV

			TH1D* Clone[NSamples-1];			

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				// Apply the additional smearing matrix Ac
				Clone[WhichSample-1] = PlotsTrue[WhichSample][WhichPlot];
				
				// Divide by the bin width
				Reweight(Clone[WhichSample-1],1.);
				Clone[WhichSample-1]->SetLineColor(Colors[WhichSample]);
				Clone[WhichSample-1]->SetMarkerColor(Colors[WhichSample]);

				PrettyPlot(Clone[WhichSample-1]); // includes scaling factor for multi dimensional analysis
				Clone[WhichSample-1]->SetLineWidth(3);		
				Clone[WhichSample-1]->Draw("hist same");		

				TLegendEntry* lGenie = legMC->AddEntry(Clone[WhichSample-1],Labels[WhichSample],"l");
				lGenie->SetTextColor(Colors[WhichSample]); 										


			}

			legMC->Draw();			

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			TString PlotNameDuplicate = PlotNames[WhichPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
			textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ ReducedPlotName ].ReplaceAll("All events","") );

			//----------------------------------------//

			// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

		//	PlotCanvas->SaveAs("/exp/uboone/data/users/apapadop/FlatTTreePlots/Atmospherics/"+Extra+"XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
		//	delete PlotCanvas;

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
