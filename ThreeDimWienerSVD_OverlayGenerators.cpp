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

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/AnalysisCode/Secondary_Code/CenterAxisTitle.cpp"
#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/AnalysisCode/Secondary_Code/SetOffsetAndSize.cpp"
#include "myFunctions.cpp"
#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/AnalysisCode/Secondary_Code/MakeMyPlotPretty.cpp"

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"

using namespace std;
using namespace Constants;

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Util.h"

//----------------------------------------//

void ThreeDimWienerSVD_OverlayGenerators(TString PlotName = "", int FirstDiscrIndex = 0, int SecondDiscrIndex = 0, bool PlotGENIE = true, bool PlotGen = false) {

	//----------------------------------------//

	Tools tools;

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		

        TString PathToFiles = "/exp/uboone/data/users/apapadop/mySTVAnalysis/myXSec/v08_00_00_52/";

	TString Extra = "";
	if (PlotGen) { Extra = "Gene"; }
	if (PlotGENIE) { Extra = "Genie"; }

	//----------------------------------------//

	vector<TString> PlotNames;
	vector< vector<double> > SliceDiscriminators;
	vector< vector< vector<double> > > SliceBinning;

	//----------------------------------------//		

	// 3D names

	PlotNames.push_back(PlotName);
	
	//----------------------------------------//	

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<TString> Runs;
//	Runs.push_back("Run1");
//	Runs.push_back("Run2");	
//	Runs.push_back("Run3");
//	Runs.push_back("Run4");
//	Runs.push_back("Run4a");
//	Runs.push_back("Run5");
	Runs.push_back("Combined");

	int NRuns = (int)(Runs.size());
	cout << "Number of Runs = " << NRuns << endl;

	//----------------------------------------//

	for (int WhichRun = 0; WhichRun < NRuns; WhichRun++) {

		vector<vector<TH1D*> > PlotsFullUncReco; PlotsFullUncReco.clear();
		vector<vector<TH1D*> > PlotsXSecReco; PlotsXSecReco.clear();
		vector<vector<TH1D*> > PlotsTotalReco; PlotsTotalReco.clear();
		vector<vector<TH1D*> > PlotsNormOnly; PlotsNormOnly.clear();		
		vector<vector<TH1D*> > PlotsReco; PlotsReco.clear();
		vector<vector<TH1D*> > PlotsCC1pReco; PlotsCC1pReco.clear();
		vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();
		vector<vector<TH1D*> > QEPlotsTrue; QEPlotsTrue.clear();
		vector<vector<TH1D*> > MECPlotsTrue; MECPlotsTrue.clear();
		vector<vector<TH1D*> > RESPlotsTrue; RESPlotsTrue.clear();
		vector<vector<TH1D*> > DISPlotsTrue; DISPlotsTrue.clear();	
		vector<vector<TH1D*> > COHPlotsTrue; COHPlotsTrue.clear();		

		gStyle->SetPalette(55); const Int_t NCont = 999; gStyle->SetNumberContours(NCont); gStyle->SetTitleSize(0.07,"t"); SetOffsetAndSize();

		vector<TString> NameOfSamples; NameOfSamples.clear();
		vector<int> Colors; Colors.clear();		
		vector<TString> Labels; Labels.clear();

		// CV

		NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("G18 ");                    

		//----------------------------------------//	

                if (PlotGENIE) {

		  NameOfSamples.push_back("GENIE_v3_4_0_AR23_20i_00_000");  Colors.push_back(GENIEv2Color); Labels.push_back("AR23 ");
		  //NameOfSamples.push_back("GENIE_v2_12_10_MEC");  Colors.push_back(GENIEv2Color); Labels.push_back("Gv2 ");
		  NameOfSamples.push_back("GENIE_v3_0_6"); Colors.push_back(Geniev3OutOfTheBoxColor); Labels.push_back("G18 NoTune ");
		  NameOfSamples.push_back("GENIE_v3_0_6_G21_11b_00_000"); Colors.push_back(SuSav2Color); Labels.push_back("G21hN ");
		  
                }

                //----------------------------------------//                                                                                                                                              

                if (PlotGen) {

		  NameOfSamples.push_back("NuWro_19_02_1"); Colors.push_back(NuWroColor); Labels.push_back("NuWro ");
		  NameOfSamples.push_back("GiBUU_2021"); Colors.push_back(GiBUUColor); Labels.push_back("GiBUU ");
		  NameOfSamples.push_back("NEUT_5_4_0_1"); Colors.push_back(kMagenta); Labels.push_back("NEUT ");
		  
                }

		//----------------------------------------//

		const int NSamples = NameOfSamples.size();
		vector<TFile*> FileSample; FileSample.clear();

		//----------------------------------------//

		TString GenXSecName = "OutputFiles/All_XSecs_3D_" + PlotNames[0] + "_" + TString(std::to_string(FirstDiscrIndex)) + "_" + Runs[WhichRun] + "_"+UBCodeVersion+".root";

		TFile* fGenXSec = TFile::Open(GenXSecName,"recreate");

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
			vector<TH1D*> QECurrentPlotsTrue; QECurrentPlotsTrue.clear();
			vector<TH1D*> MECCurrentPlotsTrue; MECCurrentPlotsTrue.clear();	
			vector<TH1D*> RESCurrentPlotsTrue; RESCurrentPlotsTrue.clear();
			vector<TH1D*> DISCurrentPlotsTrue; DISCurrentPlotsTrue.clear();	
			vector<TH1D*> COHCurrentPlotsTrue; COHCurrentPlotsTrue.clear();			

			// CV With Statistical Uncertainties

			if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

				TString FileSampleName = PathToFiles+"/WienerSVD_ExtractedXSec_"+NameOfSamples[WhichSample]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".root"; 
				FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

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

					TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QETrueSerial"+PlotNames[WhichPlot]));
					QECurrentPlotsTrue.push_back(QEhistTrue);
					TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECTrueSerial"+PlotNames[WhichPlot]));
					MECCurrentPlotsTrue.push_back(MEChistTrue);
					TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESTrueSerial"+PlotNames[WhichPlot]));
					RESCurrentPlotsTrue.push_back(REShistTrue);
					TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISTrueSerial"+PlotNames[WhichPlot]));
					DISCurrentPlotsTrue.push_back(DIShistTrue);
					TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHTrueSerial"+PlotNames[WhichPlot]));
					COHCurrentPlotsTrue.push_back(COHhistTrue);					
		
				}

			} 
			else {

                          FileSample.push_back(TFile::Open("OutputFiles/FlatTreeAnalyzerOutput_"+NameOfSamples[WhichSample]+".root"));

			  for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

			    TH1D* histTotalReco = nullptr;
			    CurrentPlotsTotalReco.push_back(histTotalReco);

			    TH1D* histXSecReco = nullptr;
			    CurrentPlotsXSecReco.push_back(histXSecReco);

			    TH1D* histFullUncReco = nullptr;
			    CurrentPlotsFullUncReco.push_back(histFullUncReco);

			    TH1D* histNormOnly = nullptr;
			    CurrentPlotsNormOnly.push_back(histNormOnly);

			    TH1D* histReco = nullptr;
			    CurrentPlotsReco.push_back(histReco);

			    TH1D* histCC1pReco = nullptr;
			    CurrentPlotsCC1pReco.push_back(histCC1pReco);

			    TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("SerialTrue"+PlotNames[WhichPlot]));
			    CurrentPlotsTrue.push_back(histTrue);

			    TH1D* QEhistTrue = (TH1D*)(FileSample[WhichSample]->Get("QESerialTrue"+PlotNames[WhichPlot]));
			    QECurrentPlotsTrue.push_back(QEhistTrue);
			    TH1D* MEChistTrue = (TH1D*)(FileSample[WhichSample]->Get("MECSerialTrue"+PlotNames[WhichPlot]));
			    MECCurrentPlotsTrue.push_back(MEChistTrue);
			    TH1D* REShistTrue = (TH1D*)(FileSample[WhichSample]->Get("RESSerialTrue"+PlotNames[WhichPlot]));
			    RESCurrentPlotsTrue.push_back(REShistTrue);
			    TH1D* DIShistTrue = (TH1D*)(FileSample[WhichSample]->Get("DISSerialTrue"+PlotNames[WhichPlot]));
			    DISCurrentPlotsTrue.push_back(DIShistTrue);
			    TH1D* COHhistTrue = (TH1D*)(FileSample[WhichSample]->Get("COHSerialTrue"+PlotNames[WhichPlot]));
			    COHCurrentPlotsTrue.push_back(COHhistTrue);

			  }

			}

			PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);
			PlotsXSecReco.push_back(CurrentPlotsXSecReco);
			PlotsTotalReco.push_back(CurrentPlotsTotalReco);
			PlotsNormOnly.push_back(CurrentPlotsNormOnly);					
			PlotsReco.push_back(CurrentPlotsReco);		
			PlotsCC1pReco.push_back(CurrentPlotsCC1pReco);
			PlotsTrue.push_back(CurrentPlotsTrue);	

			QEPlotsTrue.push_back(QECurrentPlotsTrue);			
			MECPlotsTrue.push_back(MECCurrentPlotsTrue);
			RESPlotsTrue.push_back(RESCurrentPlotsTrue);
			DISPlotsTrue.push_back(DISCurrentPlotsTrue);
			COHPlotsTrue.push_back(COHCurrentPlotsTrue);				

		}

		//----------------------------------------//

		// Loop over the plots

		vector< vector<TH1D*> > BeamOnFullUnc;
		vector< vector<TH1D*> > BeamOnXSec;
		vector< vector<TH1D*> > BeamOnStatShape;
		vector< vector<TH1D*> > BeamOnStatOnly;
		vector< vector<TH1D*> > BeamOnNormOnly;
		vector< vector< vector<TH1D*> > > MC;
		vector< vector< vector<TH1D*> > > QEMC;
		vector< vector< vector<TH1D*> > > MECMC;
		vector< vector< vector<TH1D*> > > RESMC;
		vector< vector< vector<TH1D*> > > DISMC;
		vector< vector< vector<TH1D*> > > COHMC;													

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {		

			//----------------------------------------//

			// Setting up the relevant discriminators

			SliceDiscriminators.clear();
			SliceBinning.clear();

			int StartIndex = 0;
			// Determine the global start bin number
			int GlobalIndex = 0;

			if (PlotNames[WhichPlot] == "ECal_DeltaPTDeltaAlphaTPlot") {

			  SliceDiscriminators.push_back(TwoDArrayNBinsDeltaAlphaT);
			  SliceBinning.push_back(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[FirstDiscrIndex]);

			  for (int ig = 0; ig < FirstDiscrIndex; ig++) {

			    for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[ig].size(); icolumn++) {

			      GlobalIndex += TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[ig][icolumn].size()-1;

			    }

			  }

			}

			if (PlotNames[WhichPlot] == "ECal_DeltaPtxDeltaPtyPlot") {

			  SliceDiscriminators.push_back(TwoDArrayNBinsDeltaPtx);
			  SliceBinning.push_back(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[FirstDiscrIndex]);

			  for (int ig = 0; ig < FirstDiscrIndex; ig++) {

			    for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[ig].size(); icolumn++) {

			      GlobalIndex += TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[ig][icolumn].size()-1;

			    }

			  }

			}

			if (PlotNames[WhichPlot] == "ECal_MuonCosThetaMuonMomentumPlot") {

			  SliceDiscriminators.push_back(TwoDArrayNBinsMuonMomentum);
			  SliceBinning.push_back(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[FirstDiscrIndex]);

			  for (int ig = 0; ig < FirstDiscrIndex; ig++) {

			    for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[ig].size(); icolumn++) {

			      GlobalIndex += TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[ig][icolumn].size()-1;

			    }

			  }

			}

			if (PlotNames[WhichPlot] == "ECal_ProtonCosThetaProtonMomentumPlot") {

			  SliceDiscriminators.push_back(TwoDArrayNBinsProtonMomentum);
			  SliceBinning.push_back(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[FirstDiscrIndex]);

			  for (int ig = 0; ig < FirstDiscrIndex; ig++) {

			    for (int icolumn = 0; icolumn < TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[ig].size(); icolumn++) {

			      GlobalIndex += TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[ig][icolumn].size()-1;

			    }

			  }

			}

			//----------------------------------------//

			BeamOnFullUnc.resize(N1DPlots);
			BeamOnXSec.resize(N1DPlots);
			BeamOnStatShape.resize(N1DPlots);
			BeamOnStatOnly.resize(N1DPlots);
			BeamOnNormOnly.resize(N1DPlots);							
			MC.resize(N1DPlots);
			QEMC.resize(N1DPlots);
			MECMC.resize(N1DPlots);
			RESMC.resize(N1DPlots);
			DISMC.resize(N1DPlots);
			COHMC.resize(N1DPlots);															

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("AcSerial"+PlotNames[WhichPlot]);

			TH2D* Cov = (TH2D*)FileSample[0]->Get("UnfCovSerial"+PlotNames[WhichPlot]);	
			Cov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial" + PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* NormCov = (TH2D*)FileSample[0]->Get("NormUnfCovSerial"+PlotNames[WhichPlot]);	
			NormCov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial"+PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* ShapeCov = (TH2D*)FileSample[0]->Get("ShapeUnfCovSerial"+PlotNames[WhichPlot]);	
			ShapeCov->Scale(1./TMath::Power(MultiDimScaleFactor["Serial"+PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			//----------------------------------------//

//			TH1D* UncHist = (TH1D*)(fUnc->Get("UnfUnc_"+PlotNames[WhichPlot]));

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
			      if (iBin == 0) { SerialVectorLowBin.push_back(GlobalIndex + BinCounter); }
			      // Last bin number for a given slice
			      if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(GlobalIndex + BinCounter); }

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
			QEMC[WhichPlot].resize(NSlices);
			MECMC[WhichPlot].resize(NSlices);
			RESMC[WhichPlot].resize(NSlices);
			DISMC[WhichPlot].resize(NSlices);	
			COHMC[WhichPlot].resize(NSlices);																	

			//------------------------------------//
		       
			int BinStartIndex = 0;			

			//------------------------------------//

			// MC multiplication by the additional matrix Ac
			// NOT for index 0 = overlay, that has already been done

			for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				PlotsTrue[WhichSample][WhichPlot] = Multiply(PlotsTrue[WhichSample][WhichPlot],Ac);
				QEPlotsTrue[WhichSample][WhichPlot] = Multiply(QEPlotsTrue[WhichSample][WhichPlot],Ac);
				MECPlotsTrue[WhichSample][WhichPlot] = Multiply(MECPlotsTrue[WhichSample][WhichPlot],Ac);
				RESPlotsTrue[WhichSample][WhichPlot] = Multiply(RESPlotsTrue[WhichSample][WhichPlot],Ac);
				DISPlotsTrue[WhichSample][WhichPlot] = Multiply(DISPlotsTrue[WhichSample][WhichPlot],Ac);
				COHPlotsTrue[WhichSample][WhichPlot] = Multiply(COHPlotsTrue[WhichSample][WhichPlot],Ac);																					
			}							

			//------------------------------------//

			// Loop over the N-dimensional slices
			
			for (int NDimSlice = 0; NDimSlice < NSlices; NDimSlice++) {	

				//------------------------------------//

				MC[WhichPlot][NDimSlice].resize(NSamples);
				QEMC[WhichPlot][NDimSlice].resize(NSamples);
				MECMC[WhichPlot][NDimSlice].resize(NSamples);	
				RESMC[WhichPlot][NDimSlice].resize(NSamples);
				DISMC[WhichPlot][NDimSlice].resize(NSamples);	
				COHMC[WhichPlot][NDimSlice].resize(NSamples);																		

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

				NameCopy = "Serial" + NameCopy + "_" + TString(std::to_string(SecondDiscrIndex + NDimSlice));

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

				TString CanvasName = Extra + "_" + PlotNames[WhichPlot] + "_Slice_"+TString(std::to_string(SecondDiscrIndex + NDimSlice));
				TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
				PlotCanvas->cd();
				PlotCanvas->SetBottomMargin(0.14);
				PlotCanvas->SetTopMargin(0.12);
				PlotCanvas->SetLeftMargin(0.19);
				PlotCanvas->SetRightMargin(0.03);				

				TLegend* leg = new TLegend(0.62,0.52,0.72,0.85);
				//TLegend* legChi2 = new TLegend(0.8,0.72,0.9,0.85);
				if (
				    (PlotNames[WhichPlot] == "ECal_ProtonCosThetaProtonMomentumPlot" && NDimSlice == 2) ||
				    (PlotNames[WhichPlot] == "ECal_MuonCosThetaMuonMomentumPlot" && NDimSlice == 2)
				    ) { 
				  
				  leg = new TLegend(0.22,0.52,0.32,0.85); 
				  
				}

				leg->SetBorderSize(0);
				leg->SetTextSize(0.05);
				leg->SetTextFont(FontStyle);
				leg->SetNColumns(1);
				leg->SetMargin(0.15);
				leg->SetFillStyle(0);

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

				//------------------------------//									
		
				// Alternative MC

				for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

				  MC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(PlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->SetLineColor(Colors[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->SetMarkerColor(Colors[WhichSample]);
					MC[WhichPlot][NDimSlice][WhichSample]->Draw("hist same");
				  				
	
					QEMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(QEPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					MECMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(MECPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);										
					RESMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(RESPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					DISMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(DISPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					COHMC[WhichPlot][NDimSlice][WhichSample] = tools.GetHistoBins(COHPlotsTrue[WhichSample][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning, NameOfSamples[WhichSample]);
					
					CalcChiSquared(MC[WhichPlot][NDimSlice][WhichSample],BeamOnStatShape[WhichPlot][NDimSlice],SliceCovMatrix,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
					TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[WhichSample],1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";
					TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][WhichSample],Labels[WhichSample] + Chi2NdofAlt,"l");
					lGenie->SetTextColor(Colors[WhichSample]); 										
					
				}			

				//------------------------------//

				// Overlay
				
				MC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(PlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				PrettyPlot(MC[WhichPlot][NDimSlice][0]);
				MC[WhichPlot][NDimSlice][0]->SetLineColor(Colors[0]);
				MC[WhichPlot][NDimSlice][0]->SetMarkerColor(Colors[0]);	
				MC[WhichPlot][NDimSlice][0]->Draw("hist same");	
				
				QEMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(QEPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				MECMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(MECPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");								
				RESMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(RESPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				DISMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(DISPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");
				COHMC[WhichPlot][NDimSlice][0] = tools.GetHistoBins(COHPlotsTrue[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"Overlay");

				CalcChiSquared(MC[WhichPlot][NDimSlice][0],BeamOnStatShape[WhichPlot][NDimSlice],SliceCovMatrix,Chi2[0],Ndof[0],pval[0]);
				TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";
				TLegendEntry* lGenie = leg->AddEntry(MC[WhichPlot][NDimSlice][0],Labels[0] + Chi2NdofAlt,"l");
				lGenie->SetTextColor(Colors[0]); 										

				//------------------------------//	
				
				// Stat Unc Only	
				
				BeamOnStatOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsTotalReco[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"StatOnly");											
				PrettyPlot(BeamOnStatOnly[WhichPlot][NDimSlice]);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetMarkerColor(BeamOnColor);
				BeamOnStatOnly[WhichPlot][NDimSlice]->SetLineWidth(1);			
				BeamOnStatOnly[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Stat Only

				// Plot again on top
				BeamOnStatShape[WhichPlot][NDimSlice]->Draw("e1x0 same"); // Total Unc (Shape + Stat)				

				//------------------------------//
				
				// Norm Unc Only

				BeamOnNormOnly[WhichPlot][NDimSlice] = tools.GetHistoBins(PlotsNormOnly[0][WhichPlot],SerialVectorLowBin.at(NDimSlice),SerialVectorHighBin.at(NDimSlice), MultiDimScaleFactor[ MapUncorCor[ NameCopy ] ], SerialSliceBinning,"NormOnly");
				PrettyPlot(BeamOnNormOnly[WhichPlot][NDimSlice]); // includes scaling factor for multi dimensional analysis			
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetFillColorAlpha(kGray+1, 0.45);	
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetLineColor(kGray+1);
				BeamOnNormOnly[WhichPlot][NDimSlice]->SetMarkerColor(kGray+1);
				BeamOnNormOnly[WhichPlot][NDimSlice]->Draw("e2 same");		

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
				if (Runs[WhichRun] == "Run4") { tor860_wcut = Fulltor860_wcut_Run4; }
				if (Runs[WhichRun] == "Run4a") { tor860_wcut = Fulltor860_wcut_Run4a; }				
				if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
				if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }
				TString Label = ToString(tor860_wcut)+" POT";			

				// ---------------------------------------------------------------------------------------------------------
				// ---------------------------------------------------------------------------------------------------------
				
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"MicroBooNE Data","ep");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],"(Stat #oplus Shape Unc)","");
				leg->AddEntry(BeamOnStatShape[WhichPlot][NDimSlice],Label,"");
				leg->AddEntry(BeamOnNormOnly[WhichPlot][NDimSlice],"Norm Unc","f");
				leg->Draw();			
								
				TLatex *textSlice = new TLatex();
				textSlice->SetTextFont(FontStyle);
				textSlice->SetTextSize(0.06);
				TString PlotNameDuplicate = NameCopy;
				TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
				textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ MapUncorCor[ReducedPlotName] ]);	

				//------------------------------------//

				// Storing xsecs in file

				fGenXSec->cd();

				// Unfolded covariance matrix
				SliceCovMatrix->Write("UnfCov_" + NameCopy);
				SliceAc->Write("Ac_" + NameCopy);

				// Data
				BeamOnStatShape[WhichPlot][NDimSlice]->Write("StatShape_" + NameCopy);
				BeamOnNormOnly[WhichPlot][NDimSlice]->Write("NormOnly_" + NameCopy);
				BeamOnStatOnly[WhichPlot][NDimSlice]->Write("StatOnly_" + NameCopy);
				BeamOnXSec[WhichPlot][NDimSlice]->Write("XSecOnly_" + NameCopy);
				BeamOnFullUnc[WhichPlot][NDimSlice]->Write("FullUnc_" + NameCopy);

				// Overlay GENIE
				MC[WhichPlot][NDimSlice][0]->Write("OverlayGENIE_" + NameCopy);
				QEMC[WhichPlot][NDimSlice][0]->Write("QEOverlayGENIE_" + NameCopy);
				MECMC[WhichPlot][NDimSlice][0]->Write("MECOverlayGENIE_" + NameCopy);
				RESMC[WhichPlot][NDimSlice][0]->Write("RESOverlayGENIE_" + NameCopy);
				DISMC[WhichPlot][NDimSlice][0]->Write("DISOverlayGENIE_" + NameCopy);
				COHMC[WhichPlot][NDimSlice][0]->Write("COHOverlayGENIE_" + NameCopy);

				// Store the remaining generator xsecs
				// Start from 1, the GENIE Overlay = 0 has been stored above

				for (int igen = 1; igen < NSamples; igen++ ) {

				  TString GenName = NameOfSamples[igen] + "_" + NameCopy;
				  MC[WhichPlot][NDimSlice][igen]->Write(GenName);
				  QEMC[WhichPlot][NDimSlice][igen]->Write("QE" + GenName);
				  MECMC[WhichPlot][NDimSlice][igen]->Write("MEC" + GenName);
				  RESMC[WhichPlot][NDimSlice][igen]->Write("RES" + GenName);
				  DISMC[WhichPlot][NDimSlice][igen]->Write("DIS" + GenName);
				  COHMC[WhichPlot][NDimSlice][igen]->Write("COH" + GenName);

				}

				//------------------------------------//
				
				// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

                                PlotCanvas->SaveAs("/uboone/data/users/apapadop/FlatTTreePlots/ThreeDXSec/ThreeDXSections_"+CanvasName+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");

				delete PlotCanvas;

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
