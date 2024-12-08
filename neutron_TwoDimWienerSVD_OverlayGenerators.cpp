#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../myClasses/myFunctions.cpp"
#include "../myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "../myClasses/Util.h"

//----------------------------------------//

void neutron_TwoDimWienerSVD_OverlayGenerators() {

	//----------------------------------------//

	Tools tools;

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		
	gStyle->SetOptStat(0);

        TString PathToFiles = "/exp/uboone/data/users/"+UserID+"/Atmospherics/myXSec/v08_00_00_70/";

	//----------------------------------------//

	// CV Flux File

	TFile* FluxFile = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

	//----------------------------------------//

	vector<TString> PlotNames;
	vector< vector<double> > SliceDiscriminators;
	vector< vector< vector<double> > > SliceBinning;

	//----------------------------------------//		

	// 2D analysis

	PlotNames.push_back("ThetaVis_PMissPlot"); 
	
	//----------------------------------------//	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		double DataPOT = PeLEE_ReturnBeamOnRunPOT(Runs[WhichRun]);						
		double IntegratedFlux = (HistoFlux->Integral() * DataPOT / POTPerSpill / Nominal_UB_XY_Surface);

		vector<vector<TH1D*> > PlotsFullUncReco; PlotsFullUncReco.clear();
		vector<vector<TH1D*> > PlotsXSecReco; PlotsXSecReco.clear();
		vector<vector<TH1D*> > PlotsTotalReco; PlotsTotalReco.clear();
		vector<vector<TH1D*> > PlotsNormOnly; PlotsNormOnly.clear();		
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > zeroneutronPlotsTrue; zeroneutronPlotsTrue.clear();
		vector<vector<TH1D*> > oneneutronPlotsTrue; oneneutronPlotsTrue.clear();
		vector<vector<TH1D*> > twoneutronsPlotsTrue; twoneutronsPlotsTrue.clear();
		vector<vector<TH1D*> > threeneutronsPlotsTrue; threeneutronsPlotsTrue.clear();	

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t");

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();
		vector<int> LineStyle; LineStyle.clear();

		vector<TString> weighted; weighted.clear();

		// CV

		NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("G18T "); LineStyle.push_back(G18LineStyle); weighted.push_back("");

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();
		vector<TFile*> MCFileSample; MCFileSample.clear();

		//----------------------------------------//

		// Open the files and grap the relevant plots

		for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

			vector<TH1D*> CurrentPlotsFullUncReco; CurrentPlotsFullUncReco.clear();
			vector<TH1D*> CurrentPlotsXSecReco; CurrentPlotsXSecReco.clear();
			vector<TH1D*> CurrentPlotsTotalReco; CurrentPlotsTotalReco.clear();
			vector<TH1D*> CurrentPlotsNormOnly; CurrentPlotsNormOnly.clear();			
			vector<TH1D*> CurrentPlotsReco; CurrentPlotsReco.clear();
			vector<TH1D*> CurrentPlotsCC1pReco; CurrentPlotsCC1pReco.clear();
			vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();
			vector<TH1D*> zeroneutronCurrentPlotsTrue; zeroneutronCurrentPlotsTrue.clear();
			vector<TH1D*> oneneutronCurrentPlotsTrue; oneneutronCurrentPlotsTrue.clear();	
			vector<TH1D*> twoneutronsCurrentPlotsTrue; twoneutronsCurrentPlotsTrue.clear();
			vector<TH1D*> threeneutronsCurrentPlotsTrue; threeneutronsCurrentPlotsTrue.clear();	

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToFiles+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

				TString MCFileSampleName = "/exp/uboone/data/users/apapadop/Atmospherics/OutputFiles/"+UBCodeVersion+"/TruthSTVAnalysis_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root";
				MCFileSample.push_back(TFile::Open(MCFileSampleName,"readonly")); 

				for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatRecoSerial"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histXSecReco = (TH1D*)(FileSample[WhichSample]->Get("XSecRecoSerial"+PlotNames[WhichPlot]));
					CurrentPlotsXSecReco.push_back(histTotalReco);

					TH1D* histFullUncReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUncSerial"+PlotNames[WhichPlot]));
					CurrentPlotsFullUncReco.push_back(histFullUncReco);										

					TH1D* histNormOnly = (TH1D*)(FileSample[WhichSample]->Get("NormOnlyRecoSerial"+PlotNames[WhichPlot]));
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("RecoSerial"+PlotNames[WhichPlot]));
					CurrentPlotsReco.push_back(histReco);

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("TrueSerial"+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* zeroneutronhistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("0n_TrueSerial"+PlotNames[WhichPlot]));
					zeroneutronCurrentPlotsTrue.push_back(zeroneutronhistTrue);
					TH1D* oneneutronhistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("1n_TrueSerial"+PlotNames[WhichPlot]));
					oneneutronCurrentPlotsTrue.push_back(oneneutronhistTrue);
					TH1D* twoneutronshistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("2n_TrueSerial"+PlotNames[WhichPlot]));
					twoneutronsCurrentPlotsTrue.push_back(twoneutronshistTrue);
					TH1D* threeneutronshistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("3n_TrueSerial"+PlotNames[WhichPlot]));
					threeneutronsCurrentPlotsTrue.push_back(threeneutronshistTrue);
		
				}

			} 

			PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);
			PlotsXSecReco.push_back(CurrentPlotsXSecReco);
			PlotsTotalReco.push_back(CurrentPlotsTotalReco);
			PlotsNormOnly.push_back(CurrentPlotsNormOnly);					
			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);	

			zeroneutronPlotsTrue.push_back(zeroneutronCurrentPlotsTrue);			
			oneneutronPlotsTrue.push_back(oneneutronCurrentPlotsTrue);
			twoneutronsPlotsTrue.push_back(twoneutronsCurrentPlotsTrue);
			threeneutronsPlotsTrue.push_back(threeneutronsCurrentPlotsTrue);

		}

		//----------------------------------------//

		// Loop over the plots

		vector< vector<TH1D*> > BeamOnFullUnc;
		vector< vector<TH1D*> > BeamOnXSec;
		vector< vector<TH1D*> > BeamOnStatShape;
		vector< vector<TH1D*> > BeamOnStatOnly;
		vector< vector<TH1D*> > BeamOnNormOnly;
		vector< vector< vector<TH1D*> > > MC;
		vector< vector< vector<TH1D*> > > zeroneutronMC;
		vector< vector< vector<TH1D*> > > oneneutronMC;
		vector< vector< vector<TH1D*> > > twoneutronsMC;
		vector< vector< vector<TH1D*> > > threeneutronsMC;

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {		

			//----------------------------------------//

			// Setting up the relevant discriminators

			SliceDiscriminators.clear();
			SliceBinning.clear();

			if (PlotNames[WhichPlot] == "ThetaVis_ECalPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsECal); 
				SliceBinning.push_back(TwoDArrayNBinsThetaVisInECalSlices);

			}

			if (PlotNames[WhichPlot] == "ThetaVis_DeltaPnPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPn); 
				SliceBinning.push_back(TwoDArrayNBinsThetaVisInDeltaPnSlices);

			}

			if (PlotNames[WhichPlot] == "ThetaVis_PMissPlot") {

				SliceDiscriminators.push_back(TwoDArrayNBinsPMiss); 
				SliceBinning.push_back(TwoDArrayNBinsThetaVisInPMissSlices);

			}


			//----------------------------------------//

			BeamOnFullUnc.resize(N1DPlots);
			BeamOnXSec.resize(N1DPlots);
			BeamOnStatShape.resize(N1DPlots);
			BeamOnStatOnly.resize(N1DPlots);
			BeamOnNormOnly.resize(N1DPlots);							
			MC.resize(N1DPlots);
			zeroneutronMC.resize(N1DPlots);
			oneneutronMC.resize(N1DPlots);
			twoneutronsMC.resize(N1DPlots);
			threeneutronsMC.resize(N1DPlots);

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("AcSerial"+PlotNames[WhichPlot]);

			TH2D* Cov = (TH2D*)FileSample[0]->Get("UnfCovSerial"+PlotNames[WhichPlot]);	
			Cov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial" + PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* NormCov = (TH2D*)FileSample[0]->Get("NormUnfCovSerial"+PlotNames[WhichPlot]);	
			NormCov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial"+PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* ShapeCov = (TH2D*)FileSample[0]->Get("ShapeUnfCovSerial"+PlotNames[WhichPlot]);	
			ShapeCov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial"+PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			//----------------------------------------//

			// The covariance matrix needs to be scaled by the 2D bin width

			TH2D* CovClone = (TH2D*)(Cov->Clone()); 
			TH2D* NormCovClone = (TH2D*)(NormCov->Clone());	
			TH2D* ShapeCovClone = (TH2D*)(ShapeCov->Clone());			

			int n = Cov->GetXaxis()->GetNbins();

			for (int ix = 1; ix <= n; ix++) {

				for (int iy = 1; iy <= n; iy++) {

					double WidthX = Cov->GetXaxis()->GetBinWidth(ix);
					double WidthY = Cov->GetYaxis()->GetBinWidth(iy);

					double TwoDWidth = WidthX * WidthY;

					double BinContent = Cov->GetBinContent(ix,iy);
					// Division by bin width already included
					// Still need to include scaling due to slice range
					// That is done in Tools::Get2DHistoBins
					double NewBinContent = BinContent/TwoDWidth;

					double NormBinContent = NormCov->GetBinContent(ix,iy);
					double NormNewBinContent = NormBinContent/TwoDWidth;

					double ShapeBinContent = ShapeCov->GetBinContent(ix,iy);
					double ShapeNewBinContent = ShapeBinContent/TwoDWidth;					

					// Only for the diagonal elements
					// Add the unfolding uncertainty
					// On top of everything else
					// That is done both for the final xsec result and for the unfolded covariance
					if (ix == iy) { 
						
						// unfolded covariance matrix
//						double UnfUncBin = UncHist->GetBinContent(ix);
						double UnfUncBin = 0.;

						NewBinContent = NewBinContent + TMath::Power(UnfUncBin,2.) ;
						ShapeNewBinContent = ShapeNewBinContent + TMath::Power(UnfUncBin,2.) ;						 

						// xsec uncertainty
						double CurrentUnc = PlotsReco[0][WhichPlot]->GetBinError(ix);
						double NewError = TMath::Sqrt( TMath::Power(CurrentUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsReco[0][WhichPlot]->SetBinError(ix,NewError);

						double CurrentFullUnc = PlotsFullUncReco[0][WhichPlot]->GetBinError(ix);
						double NewFullError = TMath::Sqrt( TMath::Power(CurrentFullUnc,2.) + TMath::Power(UnfUncBin,2.) ) ;
						PlotsFullUncReco[0][WhichPlot]->SetBinError(ix,NewFullError);						
						
					}

					CovClone->SetBinContent(ix,iy,NewBinContent);
					ShapeCovClone->SetBinContent(ix,iy,ShapeNewBinContent);
					NormCovClone->SetBinContent(ix,iy,NormNewBinContent);					

				}					

			}	

			//CovClone->Draw("coltz text");

			//------------------------------------//

			// Number of N-dimensional slices

			int NSlices = 1;
			vector<double> SerialVectorRanges;
			vector<int> SerialVectorBins;
			vector<int> SerialVectorLowBin;
			vector<int> SerialVectorHighBin;						

			int BinCounter = 1;

			// vector< vector<double> > SliceDiscriminators;
			// 1st index = how many disciminators
			// 2nd index = values / ranges of discriminators

			// How many discriminators do we have ? e.g. cos theta mu, Pmu et al
			for (int islice = 0; islice < (int)(SliceDiscriminators.size()); islice++) { 
				
				// For a given discriminator, how many slices do we have ? SliceDiscrimSize - 1
				int SliceDiscrimSize = SliceDiscriminators.at(islice).size()-1;
				NSlices *= SliceDiscrimSize; 

				for (int iSliceDiscrimSize = 0; iSliceDiscrimSize < SliceDiscrimSize; iSliceDiscrimSize++) {

					// Accessing the vector<double> with the bin ranges
					int SliceDiscrimValue = SliceBinning.at(0).at(iSliceDiscrimSize).size();

					// Storing the number of bins for a specific slice					
					SerialVectorBins.push_back(SliceDiscrimValue-1);
					for (int iBin = 0; iBin < SliceDiscrimValue; iBin++) {

						double BinValue = SliceBinning.at(0).at(iSliceDiscrimSize).at(iBin);
						// First bin number for a given slice
						if (iBin == 0) { SerialVectorLowBin.push_back(BinCounter); }
						// Last bin number for a given slice
						if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(BinCounter); }	

						// Storing the binning for a specific slice
						SerialVectorRanges.push_back(BinValue);

						// Increase the global bin counter
						// But not for the last bin
						if (iBin != SliceDiscrimValue-1) { BinCounter++; }

					} // End of the loop over the bins of a given slice
					//cout << "End of the loop over the bins of a given slice" << endl;

				} // End of the loop over the slices for a given discriminator

			} // End of the loop over the number of discriminators	

			//------------------------------------//

			BeamOnFullUnc[WhichPlot].resize(NSlices);
			BeamOnXSec[WhichPlot].resize(NSlices);
			BeamOnStatShape[WhichPlot].resize(NSlices);
			BeamOnStatOnly[WhichPlot].resize(NSlices);
			BeamOnNormOnly[WhichPlot].resize(NSlices);
			MC[WhichPlot].resize(NSlices);
			zeroneutronMC[WhichPlot].resize(NSlices);
			oneneutronMC[WhichPlot].resize(NSlices);
			twoneutronsMC[WhichPlot].resize(NSlices);
			threeneutronsMC[WhichPlot].resize(NSlices);	

			//------------------------------------//
		       
			int StartIndex = 0;
			int BinStartIndex = 0;			

			//------------------------------------//

			// Loop over the N-dimensional slices
			
			for (int NDimSlice = 0; NDimSlice < NSlices; NDimSlice++) {	

				//------------------------------------//

				MC[WhichPlot][NDimSlice].resize(NSamples);
				zeroneutronMC[WhichPlot][NDimSlice].resize(NSamples);
				oneneutronMC[WhichPlot][NDimSlice].resize(NSamples);	
				twoneutronsMC[WhichPlot][NDimSlice].resize(NSamples);
				threeneutronsMC[WhichPlot][NDimSlice].resize(NSamples);	

				//------------------------------------//

				TString NameCopy = PlotNames[WhichPlot];

				NameCopy.ReplaceAll("unf_","");
				NameCopy.ReplaceAll("TrueUnf_","");
				NameCopy.ReplaceAll("True_","");
				NameCopy.ReplaceAll("True","");
				NameCopy.ReplaceAll("NoSmearAlt","");			

				NameCopy.ReplaceAll("_Run1","");
				NameCopy.ReplaceAll("_Run2","");
				NameCopy.ReplaceAll("_Run3","");
				NameCopy.ReplaceAll("_Run4","");
				NameCopy.ReplaceAll("_Run4a","");	
				NameCopy.ReplaceAll("_Run5","");
				NameCopy.ReplaceAll("_Combined","");	

				NameCopy = "Serial" + NameCopy + "_" + TString(std::to_string(NDimSlice));	

				//------------------------------------//		
				
				// Get the number of bins and the bin ranges for the specific slice	

				int SliceNBins = SerialVectorBins.at(NDimSlice);
				std::vector<double> SerialSliceBinning;		

				for (int iBin = 0; iBin < SliceNBins+1; iBin++) { 

					double value = SerialVectorRanges.at(StartIndex+iBin);
					SerialSliceBinning.push_back(value);

				} // End of the number of bins and the bin ranges declaration	

				//------------------------------------//

				// Canvas, pads & legend				
				
				TString CanvasName = "Serial"+PlotNames[WhichPlot]+"_Slice_"+TString(std::to_string(NDimSlice));
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.14);
				PlotCanvas->SetTopMargin(0.12);
				PlotCanvas->SetLeftMargin(0.19);
				PlotCanvas->SetRightMargin(0.05);				

				TLegend* leg = new TLegend(0.55,0.67,0.9,0.85);
				TLegend* legMC = new TLegend(0.55,0.49,0.90,.67);
			       
				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(1);
				leg->SetMargin(0.15);
				leg->SetFillStyle(0);

				legMC->SetBorderSize(0);
				legMC->SetTextSize(0.05);
				legMC->SetTextFont(FontStyle);
				legMC->SetNColumns(1);
				legMC->SetMargin(0.15);
				legMC->SetFillStyle(0);

				//------------------------------------//

				// Corresponding covariance matrix

				TH2D* SliceCovMatrix = tools.Get2DHistoBins(CovClone,SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning);
				TH2D* SliceAc = tools.Get2DHistoBins(Ac,SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), 1., SerialSliceBinning, false);				 
				
				//------------------------------------//

				// BeamOn Total Uncertainty																		

				BeamOnStatShape[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatShape");	// includes scaling for 2D analysis									
				PrettyPlot(BeamOnStatShape[WhichPlot][NDimSlice]);
				
				double MaxValue = BeamOnStatShape[WhichPlot][NDimSlice]->GetMaximum();
				int MaxValueBin = LocateBinWithValue(BeamOnStatShape[WhichPlot][NDimSlice],MaxValue);
				double MaxValueError = BeamOnStatShape[WhichPlot][NDimSlice]->GetBinError(MaxValueBin);

				double MinValue = BeamOnStatShape[WhichPlot][NDimSlice]->GetMinimum();

				BeamOnStatShape[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerSize(1.);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetMarkerStyle(20);
				BeamOnStatShape[WhichPlot][NDimSlice]->SetLineWidth(1);	

				BeamOnStatShape[WhichPlot][NDimSlice]->GetXaxis()->CenterTitle();	

				BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->SetRangeUser(XSecRange[ MapUncorCor[ NameCopy ] ].first,XSecRange[ MapUncorCor[ NameCopy ] ].second);

				BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->SetTitle(VarLabel["Serial"+PlotNames[WhichPlot]]);							
				BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->CenterTitle();	
				
				BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Total Unc (Shape + Stat)
								
				//------------------------------//

				// arrays for MC NSamples

				double Chi2[NSamples];			
				int Ndof[NSamples];
				double pval[NSamples];
				double sigma[NSamples];

				//------------------------------//									
	
				// tools.GetHistoBins scales by both the bin width and the slice width
	
				// Overlay
	
				//PlotsTrue[0][WhichPlot] = Multiply(PlotsTrue[0][WhichPlot],Ac); // has already been done
				zeroneutronPlotsTrue[0][WhichPlot] = Multiply(zeroneutronPlotsTrue[0][WhichPlot],Ac);
				oneneutronPlotsTrue[0][WhichPlot] = Multiply(oneneutronPlotsTrue[0][WhichPlot],Ac);
				twoneutronsPlotsTrue[0][WhichPlot] = Multiply(twoneutronsPlotsTrue[0][WhichPlot],Ac);
				threeneutronsPlotsTrue[0][WhichPlot] = Multiply(threeneutronsPlotsTrue[0][WhichPlot],Ac);
				
				MC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(PlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				PrettyPlot(MC[WhichPlot][NDimSlice][0]);
				MC[WhichPlot][NDimSlice][0]->SetLineColor(kMagenta);
				MC[WhichPlot][NDimSlice][0]->SetMarkerColor(Colors[0]);	
				MC[WhichPlot][NDimSlice][0]->SetLineWidth(3);	
	
			
				zeroneutronMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(zeroneutronPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay_0n");
				oneneutronMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(oneneutronPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay_1n");								
				twoneutronsMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(twoneutronsPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay_2n");
				threeneutronsMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(threeneutronsPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay_3n");

				// Sacel by the integated flux and number of targets
				
				zeroneutronMC[WhichPlot][NDimSlice][0]->Scale(Units/(IntegratedFlux*NTargets) );
				oneneutronMC[WhichPlot][NDimSlice][0]->Scale(Units/(IntegratedFlux*NTargets));
				twoneutronsMC[WhichPlot][NDimSlice][0]->Scale(Units/(IntegratedFlux*NTargets));
				threeneutronsMC[WhichPlot][NDimSlice][0]->Scale(Units/(IntegratedFlux*NTargets));
		
				CalcChiSquared(MC[WhichPlot][NDimSlice][0],BeamOnStatShape[WhichPlot][NDimSlice],SliceCovMatrix,Chi2[0],Ndof[0],pval[0],sigma[0]);
				TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";
				TLegendEntry* lGenie = legMC->AddEntry(MC[WhichPlot][NDimSlice][0],Labels[0] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(Colors[0]); 										

				//------------------------------//	
				
				// Stat Unc Only	
				
				BeamOnStatOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsTotalReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatOnly");											
				PrettyPlot(BeamOnStatOnly[WhichPlot][NDimSlice]);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineWidth(1);			
				//BeamOnStatOnly[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Stat Only

				// Plot again on top
				//BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Total Unc (Shape + Stat)				
				
				//------------------------------//
				
				// Norm Unc Only

				BeamOnNormOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsNormOnly[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"NormOnly");
				PrettyPlot(BeamOnNormOnly[WhichPlot][NDimSlice]); // includes scaling factor for multi dimensional analysis			
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetFillColorAlpha(kGray+1, 0.45);	
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetLineColor(kGray+1);
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetMarkerColor(kGray+1);
				//BeamOnNormOnly[WhichPlot][NDimSlice]->Draw("e2 same");		

				//------------------------------//

				// XSec Only
				
				BeamOnXSec[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsXSecReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"XSecOnly");
				PrettyPlot(BeamOnXSec[WhichPlot][NDimSlice]); // includes scaling factor for multi dimensional analysis		
				BeamOnXSec[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnXSec[WhichPlot][NDimSlice]->SetLineWidth(1);		
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerSize(1.);
				BeamOnXSec[WhichPlot][NDimSlice]->SetMarkerStyle(20);	
				BeamOnXSec[WhichPlot][NDimSlice]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);		
				BeamOnXSec[WhichPlot][NDimSlice]->GetYaxis()->SetRangeUser(XSecRange[ MapUncorCor[ NameCopy ] ].first,XSecRange[ MapUncorCor[ NameCopy ] ].second);		     									

				//------------------------------//

				// Full Unc
				
				BeamOnFullUnc[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsFullUncReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"FullUnc");

				// -----------------------------------------------------------------------------------------------------------------			

				// Legend & Run / POT
				
				double tor860_wcut = -99.;
				if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
				if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
				if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
				if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
				if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }

				TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times10").ReplaceAll("+","^{")+"} POT";	
				if (Runs[WhichRun] == "Combined") { Label = "1.30 #times 10^{21} POT"; }
			
				// ---------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------
				
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"MicroBooNE Data","");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],Label,"");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"Stat #oplus Shape","ep");
				leg->AddEntry(BeamOnNormOnly[WhichPlot][NDimSlice],"Norm Unc","f");
				leg->Draw();			

				legMC->Draw();
								
				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(0.06);
				TString PlotNameDuplicate = NameCopy;
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				textSlice->DrawLatexNDC(0.24, 0.94, LatexLabel[ MapUncorCor[ReducedPlotName] ]);	

				//----------------------------------------//
				//----------------------------------------//

				// interaction breakdown canvas			
				// start of the loop over the mc samples

				for (int WhichSample = 0; WhichSample < NSamples; WhichSample++) {
			
					// Canvas with interaction breakdown for each generator
			
					TCanvas* inte_can = new TCanvas(NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun],NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
					inte_can->cd();
					inte_can->SetBottomMargin(0.14);
					inte_can->SetLeftMargin(0.18);

					//----------------------------------------//

					TLegend* ilegmc = new TLegend(0.39,0.59,0.83,0.89);
					ilegmc->SetBorderSize(0);
					ilegmc->SetTextSize(0.05);
					ilegmc->SetTextFont(FontStyle);
					ilegmc->SetNColumns(2);
					ilegmc->SetMargin(0.18);
					ilegmc->SetFillStyle(0);

					//----------------------------------------//

					TString THStackName = NameOfSamples[WhichSample] + "THStack_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun];
					THStack* thstack = new THStack(THStackName,"");	

					//----------------------------------------//

					// plot data for the first time

					BeamOnStatShape[WhichPlot][NDimSlice]->GetYaxis()->SetTitleOffset(1.3);
					BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // BeamOn Stat Total

					//----------------------------------------//

					// zeroneutron

					zeroneutronMC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(OverlayColor);
					zeroneutronMC[WhichPlot][NDimSlice][WhichSample]->SetFillColor(OverlayColor);
					thstack->Add(zeroneutronMC[WhichPlot][NDimSlice][WhichSample],"hist");
					thstack->Draw("same");

					//----------------------------------------//

					// oneneutron

					oneneutronMC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(kOrange-3);
					oneneutronMC[WhichPlot][NDimSlice][WhichSample]->SetFillColor(kOrange-3);
					thstack->Add(oneneutronMC[WhichPlot][NDimSlice][WhichSample],"hist");
					thstack->Draw("same");

					//----------------------------------------//

					// twoneutrons

					twoneutronsMC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(kGreen+1);
					twoneutronsMC[WhichPlot][NDimSlice][WhichSample]->SetFillColor(kGreen+1);
					thstack->Add(twoneutronsMC[WhichPlot][NDimSlice][WhichSample],"hist");
					thstack->Draw("same");

					//----------------------------------------//

					// threeneutrons

					threeneutronsMC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(kRed+1);
					threeneutronsMC[WhichPlot][NDimSlice][WhichSample]->SetFillColor(kRed+1);
					thstack->Add(threeneutronsMC[WhichPlot][NDimSlice][WhichSample],"hist");
					thstack->Draw("same");

					//----------------------------------------//
					
					// plot the data points again so that they can be on top

					BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // BeamOn Stat Total
					BeamOnStatOnly[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Stat Only
					BeamOnNormOnly[WhichPlot][NDimSlice]->Draw("e2 hist same"); // norm only	
 MC[WhichPlot][NDimSlice][0]->Draw("same");

					//----------------------------------------//
				
					TH1D* hstack = (TH1D*)(thstack->GetStack()->Last());
		
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"MicroBooNE Data","");	
					
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"","");	
					
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"1.30 #times 10^{21} POT","");
					
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"","");	
		
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"Stat#oplusShape","ep");
				
					ilegmc->AddEntry(BeamOnNormOnly[WhichPlot][NDimSlice],"Norm","f");

					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"#chi^{2}/ndf = " + to_string_with_precision(Chi2[WhichSample],1.) + "/"+ to_string_with_precision(Ndof[WhichSample],0) ,"");
					
					ilegmc->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"p-value = " + to_string_with_precision(pval[WhichSample],2.),"");
		
  					double qe_frac = zeroneutronMC[WhichPlot][NDimSlice][WhichSample]->Integral() / hstack->Integral() * 100.;
					TLegendEntry* lqe = ilegmc->AddEntry(zeroneutronMC[WhichPlot][NDimSlice][WhichSample],"0n (" + to_string_with_precision(qe_frac,1.) + "%)","f");
					lqe->SetTextColor(OverlayColor);	

					double mec_frac = oneneutronMC[WhichPlot][NDimSlice][WhichSample]->Integral() / hstack->Integral() * 100.;
					TLegendEntry* lmec = ilegmc->AddEntry(oneneutronMC[WhichPlot][NDimSlice][WhichSample],"1n (" + to_string_with_precision(mec_frac,1.) + "%)","f");
					lmec->SetTextColor(kOrange-3);	

					double res_frac = twoneutronsMC[WhichPlot][NDimSlice][WhichSample]->Integral() / hstack->Integral() * 100.;
					TLegendEntry* lres = ilegmc->AddEntry(twoneutronsMC[WhichPlot][NDimSlice][WhichSample],"2n (" + to_string_with_precision(res_frac,1.) + "%)","f");
					lres->SetTextColor(kGreen+1);	

					double dis_frac = threeneutronsMC[WhichPlot][NDimSlice][WhichSample]->Integral() / hstack->Integral() * 100.;
					TLegendEntry* ldis = ilegmc->AddEntry(threeneutronsMC[WhichPlot][NDimSlice][WhichSample],"3(+)n (" + to_string_with_precision(dis_frac,1.) + "%)","f");
					ldis->SetTextColor(kRed+1);	
				
					ilegmc->Draw();

					textSlice->DrawLatexNDC(0.19, 0.94, LatexLabel[ MapUncorCor[ReducedPlotName] ] );
					
					//----------------------------------------//

					// Save canvas with interaction breakdown

					gPad->RedrawAxis();
					inte_can->SaveAs("/exp/uboone/data/users/"+UserID+"/FlatTTreePlots/Atmospherics/neutronbreak_" + NameOfSamples[WhichSample] + "_XSections_"+CanvasName+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
					delete inte_can;

				} // end of the loop over the mc samples

				//------------------------------------//

				// Update the starting index to move to the next slice

				StartIndex += (SliceNBins+1);
				BinStartIndex += SliceNBins;
				
				//------------------------------------//
				
			} // End of the loop over the discriminators slices
			
			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//

	} // End of the loop over the runs	

} // End of the program 
