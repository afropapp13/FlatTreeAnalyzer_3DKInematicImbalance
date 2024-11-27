#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>

#include "../myClasses/myFunctions.cpp"
#include "../myClasses/Constants.h"
#include "../myClasses/STV_Tools.h"
#include "../myClasses/Tools.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 40.; // so that we can have xsecs per nucleus

	int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

	//----------------------------------------//	

	Tools tools;

	//----------------------------------------//	

        // Output root files

	// Root

	TString FileNameAndPath = "OutputFiles/" + fweights + "FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	// Plot declaration

	// Post FSI

	// Neutron breakdown

	int nneutrons = 3; //0,1,2,3
	TH1D* TrueThetaVis_NeutronMultiPlot[nneutrons];
	
	// loop over the neutrons
	
	for (int ineutron = 0; ineutron <= nneutrons; ineutron++) {
	
		TrueThetaVis_NeutronMultiPlot[ineutron] = new TH1D("NeutronMulti_" + ToString(ineutron) + "_TrueThetaVis_NeutronMultiPlot",";#theta_{vis} [deg]",25,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);

	}

	// 1D Fine Binning

	TH1D* TrueFineBinNeutronMultiplicityPlot[NInte];
	TH1D* TrueFineBinMuonCosThetaPlot[NInte];
	TH1D* TrueFineBinProtonCosThetaPlot[NInte];
	TH1D* TrueFineBinThetaVisPlot[NInte];
	TH1D* TrueFineBinCosThetaVisPlot[NInte];
	TH1D* TrueFineBinEvPlot[NInte];
	TH1D* TrueFineBinECalPlot[NInte];
	TH1D* TrueFineBinPMissPlot[NInte];

	// 1D Nominal Binning

	TH1D* TrueMuonCosThetaSingleBinPlot[NInte];	
	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueProtonCosThetaPlot[NInte];
	TH1D* TrueThetaVisPlot[NInte];
	TH1D* TrueCosThetaVisPlot[NInte];
	TH1D* TruePMissPlot[NInte];

	// 2D plots

	TH1D* TrueThetaVis_InECalTwoDPlot[NInte][TwoDNBinsECal];
	TH1D* SerialTrueThetaVis_InECalPlot[NInte];

	TH1D* TrueThetaVis_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
	TH1D* SerialTrueThetaVis_InDeltaPnPlot[NInte];

	TH1D* TrueThetaVis_InPMissTwoDPlot[NInte][TwoDNBinsPMiss];
	TH1D* SerialTrueThetaVis_InPMissPlot[NInte];


	// Relate pL_GKI to pL_vis = pT * tan(theta_vis)
	
	TH2D* PnGKIvsPnVis[NInte]; 
	TH2D* PnGKIvsPnTrue[NInte]; 
	TH2D* ThetaVisvsNeutronTheta[NInte]; 

	//----------------------------------------//

	// Pre FSI

	// 1D Fine Binning

	TH1D* NoFSITrueFineBinMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueFineBinProtonCosThetaPlot[NInte];
	TH1D* NoFSITrueFineBinThetaVisPlot[NInte];
	TH1D* NoFSITrueFineBinCosThetaVisPlot[NInte];
	TH1D* NoFSITrueFineBinPMissPlot[NInte];

	// 1D Nominal Binning

	TH1D* NoFSITrueMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueProtonCosThetaPlot[NInte];
	TH1D* NoFSITrueThetaVisPlot[NInte];
	TH1D* NoFSITrueCosThetaVisPlot[NInte];
	TH1D* NoFSITruePMissPlot[NInte];

	// 2D plots

	TH1D* NoFSITrueThetaVis_InECalTwoDPlot[NInte][TwoDNBinsECal];
	TH1D* NoFSISerialTrueThetaVis_InECalPlot[NInte];

	TH1D* NoFSITrueThetaVis_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
	TH1D* NoFSISerialTrueThetaVis_InDeltaPnPlot[NInte];

	TH1D* NoFSITrueThetaVis_InPMissTwoDPlot[NInte][TwoDNBinsPMiss];
	TH1D* NoFSISerialTrueThetaVis_InPMissPlot[NInte];


	// Relate pL_GKI to pL_vis = pT * tan(theta_vis)
	
	TH2D* NoFSIPLGKIvsPLVis[NInte]; 

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  // Post FSI

	  // 1D Fine Binning

	  TrueFineBinNeutronMultiplicityPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinNeutronMultiplicityPlot",";# neutrons",6,-0.5,5.5);
	  TrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinMuonCosThetaPlot",";cos#theta_{#mu}",20,-1.,1.);
	  TrueFineBinProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinProtonCosThetaPlot",";cos#theta_{p}",20,-1.,1.);
	  TrueFineBinThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinThetaVisPlot",";#theta_{vis} [deg]",25,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  TrueFineBinCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinCosThetaVisPlot",";cos#theta_{vis}",25,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);
	  TrueFineBinEvPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinEvPlot",";E_{#nu} [GeV]",NBinsEv,ArrayNBinEv);
	  TrueFineBinECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinECalPlot",";E^{Cal} [GeV]",NBinsEv,ArrayNBinEv);
	  TrueFineBinPMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinPMissPlot",";p_{miss} [GeV/c]",20,ArrayNBinsPMiss[0],ArrayNBinsPMiss[NBinsPMiss]);

	  // 1D Nominal Binning

	  TrueMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaSingleBinPlot",";",1,0,1);	
	  TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos#theta_{#mu}",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	  TrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonCosThetaPlot",";cos#theta_{p}",NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
	  TrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaVisPlot","#theta_{vis} [deg]",NBinsThetaVis,ArrayNBinsThetaVis);
	  TrueCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosThetaVisPlot","cos#theta_{vis}",NBinsCosThetaVis,ArrayNBinsCosThetaVis);
	  TruePMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TruePMissPlot","p_{miss} [GeV/c]",NBinsPMiss,ArrayNBinsPMiss);

	  for (int WhichECal = 0; WhichECal < TwoDNBinsECal; WhichECal++) {

		TString ThetaVisTwoDInECalLabel = "ThetaVis_ECal_"+tools.ConvertToString(TwoDArrayNBinsECal[WhichECal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[WhichECal+1])+"Plot";			
		TrueThetaVis_InECalTwoDPlot[inte][WhichECal] = new TH1D(InteractionLabels[inte]+"True"+ThetaVisTwoDInECalLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInECalSlices[WhichECal].size()-1,&TwoDArrayNBinsThetaVisInECalSlices[WhichECal][0]);

	  }	

	 SerialTrueThetaVis_InECalPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialThetaVis_ECalPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInECalSlices)[0]);

	  //--------------------------------------------------//

	  for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

		TString ThetaVisTwoDInDeltaPnLabel = "ThetaVis_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";			
		TrueThetaVis_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"True"+ThetaVisTwoDInDeltaPnLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsThetaVisInDeltaPnSlices[WhichDeltaPn][0]);

	  }	

	 SerialTrueThetaVis_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialThetaVis_DeltaPnPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInDeltaPnSlices)[0]);

	//--------------------------------------------------//
	
	for (int WhichPMiss = 0; WhichPMiss < TwoDNBinsPMiss; WhichPMiss++) {

		TString ThetaVisTwoDInPMissLabel = "ThetaVis_PMiss_"+tools.ConvertToString(TwoDArrayNBinsPMiss[WhichPMiss])+"To"+tools.ConvertToString(TwoDArrayNBinsPMiss[WhichPMiss+1])+"Plot";			
		TrueThetaVis_InPMissTwoDPlot[inte][WhichPMiss] = new TH1D(InteractionLabels[inte]+"True"+ThetaVisTwoDInPMissLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInPMissSlices[WhichPMiss].size()-1,&TwoDArrayNBinsThetaVisInPMissSlices[WhichPMiss][0]);

	  }	

	 SerialTrueThetaVis_InPMissPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSerialThetaVis_PMissPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInPMissSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInPMissSlices)[0]);


	//--------------------------------------------------//

	PnGKIvsPnVis[inte] = new TH2D(InteractionLabels[inte]+"PnGKIvsPnVis",";p_{n}^{GKI} [GeV/c];p_{n}^{vis} [GeV/c]",20,0.,0.85,20,0.,0.85); 
	PnGKIvsPnTrue[inte] = new TH2D(InteractionLabels[inte]+"PnGKIvsPnTrue",";p_{n}^{GKI} [GeV/c];p_{n}^{true} [GeV/c]",20,0.,0.85,20,0.,0.85); 
	ThetaVisvsNeutronTheta[inte] = new TH2D(InteractionLabels[inte]+"ThetaVisvsNeutronTheta",";#theta_{vis}^{reco} [deg];#theta_{n}^{true} [deg]",20,0.,180,20,0.,180); 


	  //--------------------------------------------------//

	  // Pre FSI

	  // 1D Fine Binning

	  NoFSITrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinMuonCosThetaPlot",";cos#theta_{#mu}",20,-1.,1.);
	  NoFSITrueFineBinProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinProtonCosThetaPlot",";cos#theta_{p}",20,-1.,1.);
	  NoFSITrueFineBinThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinThetaVisPlot",";#theta_{vis} [deg]",25,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  NoFSITrueFineBinCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinCosThetaVisPlot",";cos#theta_{vis}",25,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);
	  NoFSITrueFineBinPMissPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinPMissPlot",";p_{miss} [GeV/c]",20,ArrayNBinsPMiss[0],ArrayNBinsPMiss[NBinsPMiss]);

	  // 1D Nominal Binning

	  NoFSITrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueMuonCosThetaPlot",";cos#theta_{#mu}",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]);
	  NoFSITrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueProtonCosThetaPlot",";cos#theta_{p}",NBinsProtonCosTheta,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]);
	  NoFSITrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueThetaVisPlot",";#theta_{vis} [deg]",NBinsThetaVis,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  NoFSITrueCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueCosThetaVisPlot",";cos#theta_{vis}",NBinsCosThetaVis,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);
	  NoFSITruePMissPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITruePMissPlot",";p_{miss} [GeV/c]",NBinsPMiss,ArrayNBinsPMiss[0],ArrayNBinsPMiss[NBinsPMiss]);

	  //--------------------------------------------------//

	  for (int WhichECal = 0; WhichECal < TwoDNBinsECal; WhichECal++) {

		TString ThetaVisTwoDInECalLabel = "ThetaVis_ECal_"+tools.ConvertToString(TwoDArrayNBinsECal[WhichECal])+"To"+tools.ConvertToString(TwoDArrayNBinsECal[WhichECal+1])+"Plot";			
		NoFSITrueThetaVis_InECalTwoDPlot[inte][WhichECal] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+ThetaVisTwoDInECalLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInECalSlices[WhichECal].size()-1,&TwoDArrayNBinsThetaVisInECalSlices[WhichECal][0]);

	  }	

	 NoFSISerialTrueThetaVis_InECalPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueSerialThetaVis_ECalPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInECalSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInECalSlices)[0]);

	  //--------------------------------------------------//

	  for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

		TString ThetaVisTwoDInDeltaPnLabel = "ThetaVis_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";			
		NoFSITrueThetaVis_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+ThetaVisTwoDInDeltaPnLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsThetaVisInDeltaPnSlices[WhichDeltaPn][0]);

	  }	

	 NoFSISerialTrueThetaVis_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueSerialThetaVis_DeltaPnPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInDeltaPnSlices)[0]);

	NoFSIPLGKIvsPLVis[inte] = new TH2D(InteractionLabels[inte]+"NoFSIPLGKIvsPLVis",";p_{L}^{GKI} [GeV/c];p_{L}^{vis} [GeV/c]",20,-0.5,0.5,20,-0.5,0.5); 

	  //--------------------------------------------------//

	  for (int WhichPMiss = 0; WhichPMiss < TwoDNBinsPMiss; WhichPMiss++) {

		TString ThetaVisTwoDInPMissLabel = "ThetaVis_PMiss_"+tools.ConvertToString(TwoDArrayNBinsPMiss[WhichPMiss])+"To"+tools.ConvertToString(TwoDArrayNBinsPMiss[WhichPMiss+1])+"Plot";			
		NoFSITrueThetaVis_InPMissTwoDPlot[inte][WhichPMiss] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+ThetaVisTwoDInPMissLabel,LabelXAxisThetaVis,TwoDArrayNBinsThetaVisInPMissSlices[WhichPMiss].size()-1,&TwoDArrayNBinsThetaVisInPMissSlices[WhichPMiss][0]);

	  }	

	 NoFSISerialTrueThetaVis_InPMissPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueSerialThetaVis_PMissPlot",LabelXAxisThetaVis,tools.Return2DNBins(TwoDArrayNBinsThetaVisInPMissSlices),&tools.Return2DBinIndices(TwoDArrayNBinsThetaVisInPMissSlices)[0]);

	  //--------------------------------------------------//

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;

	//----------------------------------------//

	TFile* fweights_file = nullptr;
	TTree* tweights = nullptr;
	float cv_weight = -99.;

	if (fweights == "Weights") {

	  if (fOutputFile == "GENIE_v3_0_6") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_Nominal.root"); }
	  if (fOutputFile == "GENIE_v3_4_0_G18_10a_02_11a") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_10a_02_11a.root"); }
	  if (fOutputFile == "GENIE_v3_4_0_G18_10a_02_11b") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_10a_02_11b.root"); }
	  if (fOutputFile == "GENIE_v3_4_0_G18_10a_03_320") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_10a_03_320.root"); }
	  if (fOutputFile == "GENIE_v3_4_0_G18_10a_03_330") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_10a_03_330.root"); }

	  tweights = (TTree*)fweights_file->Get("GenericVectors__VARS");
	  tweights->SetBranchAddress("Weight", &cv_weight);

	}

	//----------------------------------------//

	TFile* fhondaspline_file = nullptr;
	TH1D* hspline = nullptr;

	if (fOutputFile == "GENIE_v3_0_6_BNBToHonda") {

		TFile* fspline = new TFile("OutputFiles/spline.root","readonly");
		hspline = (TH1D*)(fspline->Get("TrueFineBinEvPlot"));

	}

	if (fOutputFile == "GENIE_v3_0_6_BNBToHondaECal") {

		TFile* fspline = new TFile("OutputFiles/spline.root","readonly");
		hspline = (TH1D*)(fspline->Get("TrueFineBinECalPlot"));

	}


	//----------------------------------------//
	
	// Loop over the events
	cout << "nentries = " << nentries << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	
	
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//	

	  double t2kweight = 1.;

	  if (fweights == "Weights") {

	    tweights->GetEntry(jentry); t2kweight = cv_weight;

	  }

	  //----------------------------------------//			

	  float honda_weight = 1.;

          if (fOutputFile == "GENIE_v3_0_6_BNBToHonda") {

		int bin_honda = LocateClosetsBinWithValue(hspline,Enu_true);
		honda_weight = hspline->GetBinContent(bin_honda);

	 }

	  //----------------------------------------//			

	  double weight = fScaleFactor*Units*A*Weight * t2kweight * honda_weight;	
	  if (fOutputFile == "ACHILLES") { weight = fScaleFactor*Units*Weight; }
	  //if (jentry%1000 == 0) { std::cout << "fScaleFactor = " << fScaleFactor << ", weight = " << weight << std::endl; }

	  //	  cout << "fScaleFactor = " << fScaleFactor << ", Weight = " << Weight << endl;

	  //double weight = 1.;	
	  if (fOutputFile == "GiBUU_2023") { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 105 samples	
	  if (fOutputFile == "GiBUU_2023_medium") { weight = weight/150.; } // To increase the stats, the GiBUU sample has been produced in 73 samples	
	  if (fOutputFile == "ACHILLES") { weight = weight*1000./(40./12.); } // ACHILLES scaling still under discussion

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  // ACHILLES doesn't know how to handle the cc branch yet
	  if (fOutputFile != "ACHILLES") {
	    if (cc != 1) { continue; } // make sure that we have only CC interactions		
	  }

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, TrueHeavierMesonCounter = 0;
	  int ElectronTagging = 0, PhotonTagging = 0, neutron_counter = 0, neutron_counter_round = 0;
	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();		

	  int NoFSIProtonTagging = 0, NoFSIChargedPionTagging = 0, NoFSINeutralPionTagging = 0, NoFSIMuonTagging = 0, NoFSITrueHeavierMesonCounter = 0;
	  int NoFSIElectronTagging = 0, NoFSIPhotonTagging = 0;
	  vector <int> NoFSIProtonID; NoFSIProtonID.clear();
	  vector <int> NoFSIMuonID; NoFSIMuonID.clear();

	  //----------------------------------------//	

	  // Loop over the particles before the interaction
	  int neutron_index = -1;
	  for (int i = 0; i < ninitp; i++) {

		if (pdg_init[i] == 2112) { neutron_index = i; }

	  }

	  TVector3 init_neutron(px_init[neutron_index], py_init[neutron_index], pz_init[neutron_index]);
	  double pn_true = init_neutron.Mag();

	  // Loop over the final state particles / post FSI

	  for (int i = 0; i < nfsp; i++) {
		
	    double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
	  
            if (pdg[i] == 13 && (pf > 0.1 && pf < 1.2) ) {

	      MuonTagging ++;
	      MuonID.push_back(i);

	    }

	    if (pdg[i] == 2212 && (pf > 0.3 && pf < 1.) ) {

	      ProtonTagging ++;
	      ProtonID.push_back(i);

	      double eff_mass = TMath::Sqrt(E[i]*E[i] - pf*pf);
	      //if ( string( fOutputFile ).find("GiBUU") != std::string::npos && eff_mass < ProtonMass_GeV) { cout << "eff_mass = " << eff_mass << endl; }

	    }

	    if (pdg[i] == 2112 ) { neutron_counter++; }

	    if (fabs(pdg[i]) == 211 && pf > 0.07)  {

	      ChargedPionTagging ++;

	    }

	    if (pdg[i] == 111)  {

	      NeutralPionTagging ++;

	    }

	    if (fabs(pdg[i]) == 11)  {

	      ElectronTagging ++;

	    }

	    if (fabs(pdg[i]) == 22)  {

	      PhotonTagging ++;

	    }


	    if ( pdg[i] != NeutralPionPdg && fabs(pdg[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdg[i]) ) { TrueHeavierMesonCounter++; }


	  } // End of the loop over the final state particles / post FSI

	neutron_counter_round = neutron_counter;
	if ( neutron_counter_round > nneutrons ) { neutron_counter_round = nneutrons; }

	  //----------------------------------------//	

	  // Loop over final state particles / pre FSI

	  for (int i = 0; i < nvertp; i++) {
		
	    double pi = TMath::Sqrt( px_vert[i]*px_vert[i] + py_vert[i]*py_vert[i] + pz_vert[i]*pz_vert[i]);

	    if (pdg_vert[i] == 13 && (pi > 0.1 && pi < 1.2) ) {

	      NoFSIMuonTagging++;
	      NoFSIMuonID.push_back(i);

	    }

	    if (pdg_vert[i] == 2212 && (pi > 0.3 && pi < 1.) ) {

	      NoFSIProtonTagging++;
	      NoFSIProtonID.push_back(i);

	    }

	    if (fabs(pdg_vert[i]) == 211 && pi > 0.07)  {

	      NoFSIChargedPionTagging ++;

	    }

	    if (pdg_vert[i] == 111)  {

	      NoFSINeutralPionTagging ++;

	    }

	    if (fabs(pdg_vert[i]) == 11)  {

	      NoFSIElectronTagging ++;

	    }

	    if (fabs(pdg_vert[i]) == 22)  {

	      NoFSIPhotonTagging ++;

	    }


	    if ( pdg_vert[i] != NeutralPionPdg && fabs(pdg_vert[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdg_vert[i]) ) { NoFSITrueHeavierMesonCounter++; }


	  } // End of the loop over the final state particles / pre FSI

	  //----------------------------------------//	

	  // Classify the events based on the interaction type

	  // https://arxiv.org/pdf/2106.15809.pdf

	  int genie_mode = -1.;

	  if (fOutputFile ==  "ACHILLES") {

	    genie_mode = 1; // ACHILLES has only QE for now

	  } else {

	    if (TMath::Abs(Mode) == 1) { genie_mode = 1; } // QE
	    else if (TMath::Abs(Mode) == 2) { genie_mode = 2; } // MEC
	    else if (
		   TMath::Abs(Mode) == 10 ||
		   TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
		   TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
		   ) { genie_mode = 3; } // RES
	    else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { genie_mode = 4; } // DIS
	    else if (TMath::Abs(Mode) == 16) { genie_mode = 5;} // COH
	    else { continue; }  

	  }

	  // Feb 8 2022: Only case that is not covered is 15 = diffractive

	  //----------------------------------------//	

	  // If the signal definition post-FSI  is satisfied
	  if ( ProtonTagging == 1 && ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging == 1 && TrueHeavierMesonCounter == 0) { 

	    CounterEventsPassedSelection++;

	    // Kinematics of muon & proton in the final state

	    TLorentzVector Muon4Vector(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]], E[MuonID[0]]);
	    TLorentzVector Proton4Vector(px[ProtonID[0]], py[ProtonID[0]], pz[ProtonID[0]], E[ProtonID[0]]);

	    //----------------------------------------//

	    // Variables of interest
	    // Assign twice to keep tracl of the old values as well

	    //	    STV_Tools reco_stv_tool(Muon4Vector.Vect(),Proton4Vector.Vect(),Muon4Vector.E(),Proton4Vector.E());
	    STV_Tools reco_stv_tool(Muon4Vector.Vect(),Proton4Vector.Vect(),Muon4Vector.E(),TMath::Sqrt( TMath::Power(Proton4Vector.Rho(),2.) + TMath::Power(ProtonMass_GeV,2.) ) );

	    double MuonMomentum = Muon4Vector.Rho();
	    double ProtonMomentum = Proton4Vector.Rho();
	    double MuonCosTheta = Muon4Vector.CosTheta();
	    double ProtonCosTheta = Proton4Vector.CosTheta();
	    double DeltaPT = reco_stv_tool.ReturnPt();
	    TVector3 DeltaPT_vec = reco_stv_tool.ReturnPt_vec();
	    int pty_sgn = (DeltaPT_vec.Y() > 0) ? 1 : ( (DeltaPT_vec.Y() < 0) ? -1 : 0);
	    double DeltaPL = reco_stv_tool.ReturnPL();
	    double DeltaPn = reco_stv_tool.ReturnPn();
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double ECal = reco_stv_tool.ReturnECalMB();
	    double ThetaVis = reco_stv_tool.ReturnThetaVis(); // deg
	    double CosThetaVis = TMath::Cos(ThetaVis * TMath::Pi() / 180.);
	    double PLVis = DeltaPT / TMath::Tan(ThetaVis * TMath::Pi() / 180.) - ECal;
	    TVector3 beam_vec = reco_stv_tool.ReturnBeamVector();
	    double PnVis = TMath::Sqrt( beam_vec.Mag2() + ECal*ECal - 2*ECal*beam_vec.Mag()*CosThetaVis );
	    double init_neutron_theta_vis = init_neutron.Theta() *180./TMath::Pi() ;	
	    double pmiss = reco_stv_tool.ReturnMTilde(); // GeV/c^2
	   
	    //----------------------------------------//	

	    // Underflow / overflow
            if (ThetaVis < ArrayNBinsThetaVis[0]) { ThetaVis = (ArrayNBinsThetaVis[0] + ArrayNBinsThetaVis[1])/2.; }
            if (ThetaVis > ArrayNBinsThetaVis[NBinsThetaVis]) { ThetaVis = (ArrayNBinsThetaVis[NBinsThetaVis] + ArrayNBinsThetaVis[NBinsThetaVis-1])/2.; }

            if (ECal < ArrayNBinsECal[0]) { ECal = (ArrayNBinsECal[0] + ArrayNBinsECal[1])/2.; }
            if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1])/2.; }

            if (DeltaPn < ArrayNBinsDeltaPn[0]) { DeltaPn = (ArrayNBinsDeltaPn[0] + ArrayNBinsDeltaPn[1])/2.; }
            if (DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { DeltaPn = (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1])/2.; }

            if (pmiss < ArrayNBinsPMiss[0]) { pmiss = (ArrayNBinsPMiss[0] + ArrayNBinsPMiss[1])/2.; }
            if (pmiss > ArrayNBinsPMiss[NBinsPMiss]) { pmiss = (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1])/2.; }

	    //----------------------------------------//	

	    // 2D indices

	    int ECalTwoDIndex = tools.ReturnIndex(ECal, TwoDArrayNBinsECal);
	    int SerialThetaVisInECalIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInECalSlices,ECalTwoDIndex,ThetaVis);

	    int DeltaPnTwoDIndex = tools.ReturnIndex(DeltaPn, TwoDArrayNBinsDeltaPn);
	    int SerialThetaVisInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInDeltaPnSlices,DeltaPnTwoDIndex,ThetaVis);

	    int PMissTwoDIndex = tools.ReturnIndex( TMath::Abs(pmiss), TwoDArrayNBinsPMiss);
	    int SerialThetaVisInPMissIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInPMissSlices,PMissTwoDIndex,ThetaVis);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    TrueFineBinNeutronMultiplicityPlot[0]->Fill(neutron_counter,weight);
	    TrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueFineBinProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    TrueFineBinThetaVisPlot[0]->Fill(ThetaVis,weight);
	    TrueFineBinCosThetaVisPlot[0]->Fill(CosThetaVis,weight);
	    TrueFineBinEvPlot[0]->Fill(Enu_true,weight);
	    TrueFineBinECalPlot[0]->Fill(ECal,weight);
	    PnGKIvsPnVis[0]->Fill(DeltaPn,PnVis,weight);
	    PnGKIvsPnTrue[0]->Fill(DeltaPn,pn_true,weight);
	    ThetaVisvsNeutronTheta[0]->Fill(ThetaVis,init_neutron_theta_vis,weight);
	    TrueFineBinPMissPlot[0]->Fill(pmiss,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    TrueFineBinNeutronMultiplicityPlot[genie_mode]->Fill(neutron_counter,weight);	
	    TrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueFineBinProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    TrueFineBinThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    TrueFineBinCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);
	    TrueFineBinEvPlot[genie_mode]->Fill(Enu_true,weight);
	    TrueFineBinECalPlot[genie_mode]->Fill(ECal,weight);
	    PnGKIvsPnVis[genie_mode]->Fill(DeltaPn,PnVis,weight);
	    PnGKIvsPnTrue[genie_mode]->Fill(DeltaPn,pn_true,weight);
	    ThetaVisvsNeutronTheta[genie_mode]->Fill(ThetaVis,init_neutron_theta_vis,weight);
	    TrueFineBinPMissPlot[genie_mode]->Fill(pmiss,weight);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Nominal Binning

	    TrueThetaVis_NeutronMultiPlot[neutron_counter_round]->Fill(ThetaVis,weight);

	    TrueMuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);	
	    TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    TrueThetaVisPlot[0]->Fill(ThetaVis,weight);
	    TrueCosThetaVisPlot[0]->Fill(CosThetaVis,weight);
	    TruePMissPlot[0]->Fill(pmiss,weight);

	    // 2D

	    TrueThetaVis_InECalTwoDPlot[0][ECalTwoDIndex]->Fill(ThetaVis,weight);
	    SerialTrueThetaVis_InECalPlot[0]->Fill(SerialThetaVisInECalIndex,weight);								
	    TrueThetaVis_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(ThetaVis,weight);
	    SerialTrueThetaVis_InDeltaPnPlot[0]->Fill(SerialThetaVisInDeltaPnIndex,weight);							
	    SerialTrueThetaVis_InPMissPlot[0]->Fill(SerialThetaVisInPMissIndex,weight);	

	    // filling in the histo based on the interaction mode

	    // 1D Nominal Binning

	    TrueMuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);	
	    TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    TrueThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    TrueCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);
	    TruePMissPlot[genie_mode]->Fill(pmiss,weight);

	    // 2D

	    TrueThetaVis_InECalTwoDPlot[genie_mode][ECalTwoDIndex]->Fill(ThetaVis,weight);
	    SerialTrueThetaVis_InECalPlot[genie_mode]->Fill(SerialThetaVisInECalIndex,weight);								
	    TrueThetaVis_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(ThetaVis,weight);
	    SerialTrueThetaVis_InDeltaPnPlot[genie_mode]->Fill(SerialThetaVisInDeltaPnIndex,weight);							
	    SerialTrueThetaVis_InPMissPlot[genie_mode]->Fill(SerialThetaVisInPMissIndex,weight);

	    //----------------------------------------//

	  } // End of the post-FSI selection

	  //----------------------------------------//

	  // If the signal definition pre-FSI is satisfied
	  if ( NoFSIProtonTagging == 1 && NoFSIChargedPionTagging == 0 && NoFSINeutralPionTagging == 0 && NoFSIMuonTagging == 1 && NoFSITrueHeavierMesonCounter == 0) { 

	    // Kinematics of muon & proton in the final state pre FSI

	    TLorentzVector Muon4Vector(px_vert[NoFSIMuonID[0]], py_vert[NoFSIMuonID[0]], pz_vert[NoFSIMuonID[0]], E_vert[NoFSIMuonID[0]]);
	    TLorentzVector Proton4Vector(px_vert[NoFSIProtonID[0]], py_vert[NoFSIProtonID[0]], pz_vert[NoFSIProtonID[0]], E_vert[NoFSIProtonID[0]]);

	    //----------------------------------------//

	    // Variables of interest
	    // Assign twice so that you can keep track of the default values

	    STV_Tools reco_stv_tool(Muon4Vector.Vect(),Proton4Vector.Vect(),Muon4Vector.E(),Proton4Vector.E());

	    double MuonCosTheta = Muon4Vector.CosTheta();
	    double ProtonCosTheta = Proton4Vector.CosTheta();
	    double DeltaPT = reco_stv_tool.ReturnPt();
	    double DeltaPn = reco_stv_tool.ReturnPn();
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
            double ECal = reco_stv_tool.ReturnECalMB();
            double ThetaVis = reco_stv_tool.ReturnThetaVis(); // deg
	    double CosThetaVis = TMath::Cos(ThetaVis * TMath::Pi() / 180.);
	    TVector3 DeltaPT_vec = reco_stv_tool.ReturnPt_vec();
	    int pty_sgn = (DeltaPT_vec.Y() > 0) ? 1 : ( (DeltaPT_vec.Y() < 0) ? -1 : 0);
	    double DeltaPL = reco_stv_tool.ReturnPL();
	    double PLVis = DeltaPT / TMath::Tan(ThetaVis * TMath::Pi() / 180.) - ECal;
	    double pmiss = reco_stv_tool.ReturnMTilde();
	
	    //----------------------------------------//	

	    // Underflow / overflow
            if (ThetaVis < ArrayNBinsThetaVis[0]) { ThetaVis = (ArrayNBinsThetaVis[0] + ArrayNBinsThetaVis[1])/2.; }
            if (ThetaVis > ArrayNBinsThetaVis[NBinsThetaVis]) { ThetaVis = (ArrayNBinsThetaVis[NBinsThetaVis] + ArrayNBinsThetaVis[NBinsThetaVis-1])/2.; }

            if (ECal < ArrayNBinsECal[0]) { ECal = (ArrayNBinsECal[0] + ArrayNBinsECal[1])/2.; }
            if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1])/2.; }

            if (DeltaPn < ArrayNBinsDeltaPn[0]) { DeltaPn = (ArrayNBinsDeltaPn[0] + ArrayNBinsDeltaPn[1])/2.; }
            if (DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { DeltaPn = (ArrayNBinsDeltaPn[NBinsDeltaPn] + ArrayNBinsDeltaPn[NBinsDeltaPn-1])/2.; }

            if (pmiss < ArrayNBinsPMiss[0]) { pmiss = (ArrayNBinsPMiss[0] + ArrayNBinsPMiss[1])/2.; }
            if (pmiss > ArrayNBinsPMiss[NBinsPMiss]) { pmiss = (ArrayNBinsPMiss[NBinsPMiss] + ArrayNBinsPMiss[NBinsPMiss-1])/2.; }

	    //----------------------------------------//	

	    // 2D indices

	    int ECalTwoDIndex = tools.ReturnIndex(ECal, TwoDArrayNBinsECal);
	    int SerialThetaVisInECalIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInECalSlices,ECalTwoDIndex,ThetaVis);

	    int DeltaPnTwoDIndex = tools.ReturnIndex(DeltaPn, TwoDArrayNBinsDeltaPn);
	    int SerialThetaVisInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInDeltaPnSlices,DeltaPnTwoDIndex,ThetaVis);

	    int PMissTwoDIndex = tools.ReturnIndex( TMath::Abs(pmiss), TwoDArrayNBinsPMiss);
	    int SerialThetaVisInPMissIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsThetaVisInPMissSlices,PMissTwoDIndex,ThetaVis);


	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    NoFSITrueFineBinThetaVisPlot[0]->Fill(ThetaVis,weight);
	    NoFSITrueFineBinCosThetaVisPlot[0]->Fill(CosThetaVis,weight);
	    NoFSITrueFineBinPMissPlot[0]->Fill(pmiss,weight);
	    NoFSIPLGKIvsPLVis[0]->Fill(DeltaPL,PLVis,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    NoFSITrueFineBinThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    NoFSITrueFineBinCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);
	    NoFSITrueFineBinPMissPlot[genie_mode]->Fill(pmiss,weight);

	    //----------------------------------------//

	    // 1D Nominal Binning

	    NoFSITrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    NoFSITrueThetaVisPlot[0]->Fill(ThetaVis,weight);
	    NoFSITrueCosThetaVisPlot[0]->Fill(CosThetaVis,weight);
	    NoFSITruePMissPlot[0]->Fill(pmiss,weight);

	    NoFSITrueThetaVis_InECalTwoDPlot[0][ECalTwoDIndex]->Fill(ThetaVis,weight);
	    NoFSISerialTrueThetaVis_InECalPlot[0]->Fill(SerialThetaVisInECalIndex,weight);								
	    NoFSITrueThetaVis_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(ThetaVis,weight);
	    NoFSISerialTrueThetaVis_InDeltaPnPlot[0]->Fill(SerialThetaVisInDeltaPnIndex,weight);							
	    NoFSISerialTrueThetaVis_InPMissPlot[0]->Fill(SerialThetaVisInPMissIndex,weight);							

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    NoFSITrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    NoFSITrueProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    NoFSITrueThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    NoFSITrueCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);
	    NoFSITruePMissPlot[genie_mode]->Fill(pmiss,weight);
	    NoFSIPLGKIvsPLVis[genie_mode]->Fill(DeltaPL,PLVis,weight);

	    // 2D

	    NoFSITrueThetaVis_InECalTwoDPlot[genie_mode][ECalTwoDIndex]->Fill(ThetaVis,weight);
	    NoFSISerialTrueThetaVis_InECalPlot[genie_mode]->Fill(SerialThetaVisInECalIndex,weight);								
	    NoFSITrueThetaVis_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(ThetaVis,weight);
	    NoFSISerialTrueThetaVis_InDeltaPnPlot[genie_mode]->Fill(SerialThetaVisInDeltaPnIndex,weight);							
	    NoFSISerialTrueThetaVis_InPMissPlot[genie_mode]->Fill(SerialThetaVisInPMissIndex,weight);							

	  } // End of the post-FSI selection

	  //----------------------------------------//
	
	} // End of the loop over the events

	//----------------------------------------//	

	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	//----------------------------------------//	
	//----------------------------------------//	

	// The serial plots don't need to be rebinned because the bin width is 1

	// The nominal bin plots also shouldn't be divided by the bin width because we first need to multiply by Ac and then divide by the bin width

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		//----------------------------------------//
	
	        // 1D Fine Binning Post FSI

		Reweight(TrueFineBinMuonCosThetaPlot[inte]);
		Reweight(TrueFineBinProtonCosThetaPlot[inte]);
		Reweight(TrueFineBinThetaVisPlot[inte]);
		Reweight(TrueFineBinCosThetaVisPlot[inte]);
		Reweight(TrueFineBinEvPlot[inte]);
		Reweight(TrueFineBinECalPlot[inte]);
		Reweight(TrueFineBinPMissPlot[inte]);

	        // 1D Fine Binning Pre FSI

		Reweight(NoFSITrueFineBinMuonCosThetaPlot[inte]);
		Reweight(NoFSITrueFineBinProtonCosThetaPlot[inte]);
		Reweight(NoFSITrueFineBinThetaVisPlot[inte]);
		Reweight(NoFSITrueFineBinCosThetaVisPlot[inte]);
		Reweight(NoFSITrueFineBinPMissPlot[inte]);

		// 2D analysis

		for (int iecal = 0; iecal < TwoDNBinsECal; iecal++) {

			Reweight(TrueThetaVis_InECalTwoDPlot[inte][iecal]);
			Reweight(NoFSITrueThetaVis_InECalTwoDPlot[inte][iecal]);

		}

		for (int ideltapn = 0; ideltapn < TwoDNBinsECal; ideltapn++) {

			Reweight(TrueThetaVis_InDeltaPnTwoDPlot[inte][ideltapn]);
			Reweight(NoFSITrueThetaVis_InDeltaPnTwoDPlot[inte][ideltapn]);

		}

		for (int ipmiss = 0; ipmiss < TwoDNBinsPMiss; ipmiss++) {

			Reweight(TrueThetaVis_InPMissTwoDPlot[inte][ipmiss]);
			Reweight(NoFSITrueThetaVis_InPMissTwoDPlot[inte][ipmiss]);

		}

		//----------------------------------------//

	} // End of the loop over the interaction processes		

	for (int ineutron = 0; ineutron <= nneutrons; ineutron++) {
	
		Reweight(TrueThetaVis_NeutronMultiPlot[ineutron]);

	}

	//----------------------------------------//		

	file->cd();
	file->Write();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created" << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	//----------------------------------------//		

} // End of the program

//----------------------------------------//		
