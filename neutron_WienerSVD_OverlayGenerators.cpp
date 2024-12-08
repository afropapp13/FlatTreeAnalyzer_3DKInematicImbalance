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

void neutron_WienerSVD_OverlayGenerators() {

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		

	TString PathToFiles = "/exp/uboone/data/users/"+UserID+"/Atmospherics/myXSec/v08_00_00_70/";

	//----------------------------------------//

	// CV Flux File

	TFile* FluxFile = TFile::Open("../mySTVAnalysis/MCC9_FluxHist_volTPCActive.root"); 
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));


	//----------------------------------------//

	vector<TString> PlotNames;
	PlotNames.push_back("ThetaVisPlot");
	//PlotNames.push_back("PMissPlot");

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

					TH1D* histTotalReco = (TH1D*)(FileSample[WhichSample]->Get("StatReco"+PlotNames[WhichPlot]));
					CurrentPlotsTotalReco.push_back(histTotalReco);

					TH1D* histXSecReco = (TH1D*)(FileSample[WhichSample]->Get("XSecReco"+PlotNames[WhichPlot]));
					CurrentPlotsXSecReco.push_back(histXSecReco);

					TH1D* histFullUncReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot]));
					CurrentPlotsFullUncReco.push_back(histFullUncReco);										

					TH1D* histNormOnly = (TH1D*)(FileSample[WhichSample]->Get("NormOnlyReco"+PlotNames[WhichPlot]));
					CurrentPlotsNormOnly.push_back(histNormOnly);					

					TH1D* histReco = (TH1D*)(FileSample[WhichSample]->Get("Reco"+PlotNames[WhichPlot]));
					if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { histReco = (TH1D*)(FileSample[WhichSample]->Get("RecoFullUnc"+PlotNames[WhichPlot])); }
					CurrentPlotsReco.push_back(histReco);

					TString TrueString = "True";

					TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get(TrueString+PlotNames[WhichPlot]));
					CurrentPlotsTrue.push_back(histTrue);

					TH1D* zeroneutronhistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("0n_"+TrueString+PlotNames[WhichPlot]));
					zeroneutronCurrentPlotsTrue.push_back(zeroneutronhistTrue);
					TH1D* oneneutronhistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("1n_"+TrueString+PlotNames[WhichPlot]));
					oneneutronCurrentPlotsTrue.push_back(oneneutronhistTrue);
					TH1D* twoneutronshistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("2n_"+TrueString+PlotNames[WhichPlot]));
					twoneutronsCurrentPlotsTrue.push_back(twoneutronshistTrue);
					TH1D* threeneutronshistTrue = (TH1D*)(MCFileSample[WhichSample]->Get("3n_"+TrueString+PlotNames[WhichPlot]));
					threeneutronsCurrentPlotsTrue.push_back(threeneutronshistTrue);
		
				}

			}

			PlotsXSecReco.push_back(CurrentPlotsXSecReco);
			PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);			
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

		for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

			//----------------------------------------//

			TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac"+PlotNames[WhichPlot]);

			TString CovString = "UnfCov"+PlotNames[WhichPlot];
			//cout << CovString << endl;
			TH2D* Cov = (TH2D*)FileSample[0]->Get(CovString);	
			Cov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* NormCov = (TH2D*)FileSample[0]->Get("NormUnfCov"+PlotNames[WhichPlot]);	
			NormCov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis

			TH2D* ShapeCov = (TH2D*)FileSample[0]->Get("ShapeUnfCov"+PlotNames[WhichPlot]);	
			ShapeCov->Scale(1./TMath::Power(MultiDimScaleFactor[PlotNames[WhichPlot]],2.)); // includes scaling factor for multi dimensional analysis								

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

			// -----------------------------------------------------------------------------------------------------------------------------			

			TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_"+Runs[WhichRun],PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
			PlotCanvas->cd();

			TPad *midPad = new TPad("midPad", "",0.005, 0., 0.995, 0.995);
			midPad->SetBottomMargin(0.16);
			midPad->SetTopMargin(0.12);
			midPad->SetLeftMargin(0.19);
			midPad->SetRightMargin(0.03);			
			midPad->Draw();

			TLegend* leg = new TLegend(0.26,0.68,0.59,0.85);
			TLegend* legMC = new TLegend(0.59,0.68,0.69,0.85);
			
			if (
				PlotNames[WhichPlot] == "MuonCosThetaPlot" || 
				PlotNames[WhichPlot].Contains("Serial")
			) { 
				
			  leg = new TLegend(0.24,0.68,0.57,0.85);	
			  legMC = new TLegend(0.56,0.68,0.66,0.85);

			}

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

			// ------------------------------------------------------------------------------------------------------------------

			// BeamOn Total Uncertainty

			PrettyPlot(PlotsReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis

			double MaxValue = PlotsReco[0][WhichPlot]->GetMaximum();
			int MaxValueBin = LocateBinWithValue(PlotsReco[0][WhichPlot],MaxValue);
			double MaxValueError = PlotsReco[0][WhichPlot]->GetBinError(MaxValueBin);

			double MinValue = PlotsReco[0][WhichPlot]->GetMinimum();
													
			PlotsReco[0][WhichPlot]->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);

			PlotsReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsReco[0][WhichPlot]->SetLineWidth(1);	
			PlotsReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);					
		
			midPad->cd();
			
			//------------------------------//

			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // Total Unc (Shape + Stat)

			PrettyPlot(PlotsTotalReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsTotalReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsTotalReco[0][WhichPlot]->SetLineWidth(1);			
			PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only

			PrettyPlot(PlotsXSecReco[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsXSecReco[0][WhichPlot]->SetLineColor(BeamOnColor);
			PlotsXSecReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
			PlotsXSecReco[0][WhichPlot]->SetLineWidth(1);
			PlotsXSecReco[0][WhichPlot]->SetMarkerSize(1.);
			PlotsXSecReco[0][WhichPlot]->SetMarkerStyle(20);
			PlotsXSecReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);									
			//PlotsXSecReco[0][WhichPlot]->Draw("e1x0 same"); // XSec Only			
			
			PrettyPlot(PlotsNormOnly[0][WhichPlot]); // includes scaling factor for multi dimensional analysis			
			PlotsNormOnly[0][WhichPlot]->SetFillColorAlpha(kGray+1, 0.75);	
			PlotsNormOnly[0][WhichPlot]->SetLineColor(kGray+1);
			PlotsNormOnly[0][WhichPlot]->SetMarkerColor(kGray+1);
			if (PlotNames[WhichPlot] != "MuonCosThetaSingleBinPlot") { PlotsNormOnly[0][WhichPlot]->Draw("e2 hist same"); } // Norm unc Only					

			// -----------------------------------------------------------------------------------------------------------------

			// Overlay GENIE v3 + Tune

			PrettyPlot(PlotsTrue[0][WhichPlot]); // includes scaling factor for multi dimensional analysis
			PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);
			PlotsTrue[0][WhichPlot]->SetLineStyle(kSolid);

			// -----------------------------------------------------------------------------------------------------------------

			// arrays for NSamples

			double Chi2[NSamples];
			double ShapeChi2[NSamples];						
			int Ndof[NSamples];
			double pval[NSamples];
			double sigma[NSamples];

			//----------------------------------------//

			// Legend & Run / POT

			double tor860_wcut = -99.;
			if (Runs[WhichRun] == "Run1") { tor860_wcut = Fulltor860_wcut_Run1; }
			if (Runs[WhichRun] == "Run2") { tor860_wcut = Fulltor860_wcut_Run2; }
			if (Runs[WhichRun] == "Run3") { tor860_wcut = Fulltor860_wcut_Run3; }
			if (Runs[WhichRun] == "Run5") { tor860_wcut = Fulltor860_wcut_Run5; }
			if (Runs[WhichRun] == "Combined") { tor860_wcut = Fulltor860_wcut_Combined; }

			TString Label = ToString(tor860_wcut).ReplaceAll("e"," #times 10").ReplaceAll("+","^{")+"} POT";	
			if (Runs[WhichRun] == "Combined") { Label = "1.30 #times 10^{21} POT"; }
	
			// ---------------------------------------------------------------------------------------------------------

			// Total Chi2
			CalcChiSquared(PlotsTrue[0][WhichPlot],PlotsReco[0][WhichPlot],CovClone,Chi2[0],Ndof[0],pval[0],sigma[0]);
			TString Chi2NdofNom = "(" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";
			if (PlotNames[WhichPlot] == "MuonCosThetaSingleBinPlot") { Chi2NdofNom = ""; }

			TLegendEntry* lGenie_GenieOverlay = legMC->AddEntry(PlotsTrue[0][WhichPlot],Labels[0]+Chi2NdofNom,"l");
			PlotsTrue[0][WhichPlot]->SetLineWidth(3); 
			PlotsTrue[0][WhichPlot]->Draw("hist same"); 
			lGenie_GenieOverlay->SetTextColor(Colors[0]); 

			// ---------------------------------------------------------------------------------------------------------
			// ---------------------------------------------------------------------------------------------------------

			PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only
			PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

			leg->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","");
			leg->AddEntry(PlotsReco[0][WhichPlot],Label,"");
			leg->AddEntry(PlotsReco[0][WhichPlot],"Stat #oplus Shape","ep");			
			leg->AddEntry(PlotsNormOnly[0][WhichPlot],"Norm","f"); 
			leg->Draw();			

			legMC->Draw();			

			TLatex *textSlice = new TLatex();
			textSlice->SetTextFont(FontStyle);
			textSlice->SetTextSize(0.06);
			TString PlotNameDuplicate = PlotNames[WhichPlot];
			TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
			TString ReducedLatexLabel = LatexLabel[ ReducedPlotName ];
			textSlice->DrawLatexNDC(0.2, 0.92, ReducedLatexLabel.ReplaceAll("All events","") );

			TLatex *textPanel = new TLatex();
			textPanel->SetTextFont(FontStyle);
			textPanel->SetTextSize(TextSize);

			//----------------------------------------//
			//----------------------------------------//

			// interaction breakdown canvas			
			// start of the loop over the mc samples

			for (int WhichSample = 0; WhichSample < NSamples; WhichSample++) {
			
				// Canvas with interaction breakdown for each generator
			
				TCanvas* inte_can = new TCanvas(NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun],NameOfSamples[WhichSample]+"_" + PlotNames[WhichPlot]+"_"+Runs[WhichRun],205,34,1024,768);
				inte_can->cd();
				inte_can->SetBottomMargin(0.14);
				inte_can->SetLeftMargin(0.16);

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
				
				PlotsReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.1);
				PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total

				//----------------------------------------//

				// zeroneutron

				zeroneutronPlotsTrue[WhichSample][WhichPlot] = Multiply(zeroneutronPlotsTrue[WhichSample][WhichPlot],Ac);
				zeroneutronPlotsTrue[WhichSample][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
				Reweight(zeroneutronPlotsTrue[WhichSample][WhichPlot]);
				zeroneutronPlotsTrue[WhichSample][WhichPlot]->SetLineColor(OverlayColor);
				zeroneutronPlotsTrue[WhichSample][WhichPlot]->SetFillColor(OverlayColor);
				thstack->Add(zeroneutronPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// oneneutron
				oneneutronPlotsTrue[WhichSample][WhichPlot] = Multiply(oneneutronPlotsTrue[WhichSample][WhichPlot],Ac);
				oneneutronPlotsTrue[WhichSample][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
				Reweight(oneneutronPlotsTrue[WhichSample][WhichPlot]);
				oneneutronPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kOrange-3);
				oneneutronPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kOrange-3);
				thstack->Add(oneneutronPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// twoneutrons

				twoneutronsPlotsTrue[WhichSample][WhichPlot] = Multiply(twoneutronsPlotsTrue[WhichSample][WhichPlot],Ac);
				twoneutronsPlotsTrue[WhichSample][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
				Reweight(twoneutronsPlotsTrue[WhichSample][WhichPlot]);
				twoneutronsPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kGreen+1);
				twoneutronsPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kGreen+1);
				thstack->Add(twoneutronsPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//

				// threeneutrons
				threeneutronsPlotsTrue[WhichSample][WhichPlot] = Multiply(threeneutronsPlotsTrue[WhichSample][WhichPlot],Ac);
				threeneutronsPlotsTrue[WhichSample][WhichPlot]->Scale(Units/(IntegratedFlux*NTargets));
				Reweight(threeneutronsPlotsTrue[WhichSample][WhichPlot]);
				threeneutronsPlotsTrue[WhichSample][WhichPlot]->SetLineColor(kRed+1);
				threeneutronsPlotsTrue[WhichSample][WhichPlot]->SetFillColor(kRed+1);
				thstack->Add(threeneutronsPlotsTrue[WhichSample][WhichPlot],"hist");
				thstack->Draw("same");

				//----------------------------------------//
				
				// plot the data points again so that they can be on top

				PlotsReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total
				PlotsTotalReco[0][WhichPlot]->Draw("e1x0 same"); // Stat Only
				PlotsNormOnly[0][WhichPlot]->Draw("e2 hist same"); // Norm only	
				PlotsTrue[0][WhichPlot]->Draw("same");

				//----------------------------------------//
			
				TH1D* hstack = (TH1D*)(thstack->GetStack()->Last());
	
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"MicroBooNE Data","");	
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"","");	
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"1.30 #times 10^{21} POT","");
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"","");	

				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"Stat#oplusShape","ep");
			
				ilegmc->AddEntry(PlotsNormOnly[0][WhichPlot],"Norm","f");
	
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"#chi^{2}/ndf = " + to_string_with_precision(Chi2[WhichSample],1.) + "/"+ to_string_with_precision(Ndof[WhichSample],0) ,"");
				
				ilegmc->AddEntry(PlotsReco[0][WhichPlot],"p-value = " + to_string_with_precision(pval[WhichSample],2.),"");
				
				double qe_frac = zeroneutronPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lqe = ilegmc->AddEntry(zeroneutronPlotsTrue[WhichSample][WhichPlot],"0n (" + to_string_with_precision(qe_frac,1.) + "%)","f");
				lqe->SetTextColor(OverlayColor);	

				double mec_frac = oneneutronPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lmec = ilegmc->AddEntry(oneneutronPlotsTrue[WhichSample][WhichPlot],"1n (" + to_string_with_precision(mec_frac,1.) + "%)","f");
				lmec->SetTextColor(kOrange-3);	

				double res_frac = twoneutronsPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* lres = ilegmc->AddEntry(twoneutronsPlotsTrue[WhichSample][WhichPlot],"2n (" + to_string_with_precision(res_frac,1.) + "%)","f");
				lres->SetTextColor(kGreen+1);	

				double dis_frac = threeneutronsPlotsTrue[WhichSample][WhichPlot]->Integral() / hstack->Integral() * 100.;
				TLegendEntry* ldis = ilegmc->AddEntry(threeneutronsPlotsTrue[WhichSample][WhichPlot],"3(+)n (" + to_string_with_precision(dis_frac,1.) + "%)","f");
				ldis->SetTextColor(kRed+1);	

				textSlice->DrawLatexNDC(0.17, 0.92, LatexLabel[ ReducedPlotName ] );
				ilegmc->Draw();

				//----------------------------------------//

				// Save canvas with interaction breakdown
	
				gPad->RedrawAxis();
				inte_can->SaveAs("/exp/uboone/data/users/"+UserID+"/FlatTTreePlots/Atmospherics/neutronbreak_" + NameOfSamples[WhichSample] + "_XSections_"+PlotNames[WhichPlot]+"_"+Runs[WhichRun]+"_"+UBCodeVersion+".pdf");
				delete inte_can;

			} // end of the loop over the mc samples

			//----------------------------------------//

		} // End of the loop over the plots

		//----------------------------------------//					

	} // End of the loop over the runs	

} // End of the program 
