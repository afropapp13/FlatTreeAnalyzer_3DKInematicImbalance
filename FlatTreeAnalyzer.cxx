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

	// 1D Fine Binning

	TH1D* TrueFineBinMuonCosThetaPlot[NInte];
	TH1D* TrueFineBinProtonCosThetaPlot[NInte];
	TH1D* TrueFineBinThetaVisPlot[NInte];
	TH1D* TrueFineBinCosThetaVisPlot[NInte];

	// 1D Nominal Binning

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueProtonCosThetaPlot[NInte];
	TH1D* TrueThetaVisPlot[NInte];
	TH1D* TrueCosThetaVisPlot[NInte];

	//----------------------------------------//

	// Pre FSI

	// 1D Fine Binning

	TH1D* NoFSITrueFineBinMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueFineBinProtonCosThetaPlot[NInte];
	TH1D* NoFSITrueFineBinThetaVisPlot[NInte];
	TH1D* NoFSITrueFineBinCosThetaVisPlot[NInte];

	// 1D Nominal Binning

	TH1D* NoFSITrueMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueProtonCosThetaPlot[NInte];
	TH1D* NoFSITrueThetaVisPlot[NInte];
	TH1D* NoFSITrueCosThetaVisPlot[NInte];

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  // Post FSI

	  // 1D Fine Binning

	  TrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinMuonCosThetaPlot",";cos#theta_{#mu}",20,-1.,1.);
	  TrueFineBinProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinProtonCosThetaPlot",";cos#theta_{p}",20,-1.,1.);
	  TrueFineBinThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinThetaVisPlot",";#theta_{vis} [deg]",25,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  TrueFineBinCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinCosThetaVisPlot",";cos#theta_{vis}",25,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);

	  // 1D Nominal Binning

	  TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos#theta_{#mu}",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	  TrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueProtonCosThetaPlot",";cos#theta_{p}",NBinsProtonCosTheta,ArrayNBinsProtonCosTheta);
	  TrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaVisPlot","#theta_{vis} [deg]",NBinsThetaVis,ArrayNBinsThetaVis);
	  TrueCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosThetaVisPlot","cos#theta_{vis}",NBinsCosThetaVis,ArrayNBinsCosThetaVis);

	  //--------------------------------------------------//

	  // Pre FSI

	  // 1D Fine Binning

	  NoFSITrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinMuonCosThetaPlot",";cos#theta_{#mu}",20,-1.,1.);
	  NoFSITrueFineBinProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinProtonCosThetaPlot",";cos#theta_{p}",20,-1.,1.);
	  NoFSITrueFineBinThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinThetaVisPlot",";#theta_{vis} [deg]",25,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  NoFSITrueFineBinCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinCosThetaVisPlot",";cos#theta_{vis}",25,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);

	  // 1D Nominal Binning

	  NoFSITrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueMuonCosThetaPlot",";cos#theta_{#mu}",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]);
	  NoFSITrueProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueProtonCosThetaPlot",";cos#theta_{p}",NBinsProtonCosTheta,ArrayNBinsProtonCosTheta[0],ArrayNBinsProtonCosTheta[NBinsProtonCosTheta]);
	  NoFSITrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueThetaVisPlot",";#theta_{vis} [deg]",NBinsThetaVis,ArrayNBinsThetaVis[0],ArrayNBinsThetaVis[NBinsThetaVis]);
	  NoFSITrueCosThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueCosThetaVisPlot",";cos#theta_{vis}",NBinsCosThetaVis,ArrayNBinsCosThetaVis[0],ArrayNBinsCosThetaVis[NBinsCosThetaVis]);

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

	  double weight = fScaleFactor*Units*A*Weight * t2kweight;	
	  if (fOutputFile == "ACHILLES") { weight = fScaleFactor*Units*Weight; }
	  //if (jentry%1000 == 0) { std::cout << "fScaleFactor = " << fScaleFactor << ", weight = " << weight << std::endl; }

	  //	  cout << "fScaleFactor = " << fScaleFactor << ", Weight = " << Weight << endl;

	  //double weight = 1.;	
	  if (fOutputFile == "GiBUU_2023") { weight = weight/500.; } // To increase the stats, the GiBUU sample has been produced in 50 samples	
	  if (fOutputFile == "ACHILLES") { weight = weight*1000./(40./12.); } // ACHILLES scaling still under discussion

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  // ACHILLES doesn't know how to handle the cc branch yet
	  if (fOutputFile != "ACHILLES") {
	    if (cc != 1) { continue; } // make sure that we have only CC interactions		
	  }

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, TrueHeavierMesonCounter = 0;
	  int ElectronTagging = 0, PhotonTagging = 0;
	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();		

	  int NoFSIProtonTagging = 0, NoFSIChargedPionTagging = 0, NoFSINeutralPionTagging = 0, NoFSIMuonTagging = 0, NoFSITrueHeavierMesonCounter = 0;
	  int NoFSIElectronTagging = 0, NoFSIPhotonTagging = 0;
	  vector <int> NoFSIProtonID; NoFSIProtonID.clear();
	  vector <int> NoFSIMuonID; NoFSIMuonID.clear();

	  //----------------------------------------//	

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

	    }

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
	  if ( ProtonTagging == 1 && ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging == 1 && TrueHeavierMesonCounter == 0 && ElectronTagging == 0 && PhotonTagging == 0) { 

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
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double ECal = reco_stv_tool.ReturnECalMB();
	    double ThetaVis = reco_stv_tool.ReturnThetaVis(); // deg
	    double CosThetaVis = TMath::Cos(ThetaVis * TMath::Pi() / 180.);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    TrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueFineBinProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    TrueFineBinThetaVisPlot[0]->Fill(ThetaVis,weight);
	    TrueFineBinCosThetaVisPlot[0]->Fill(CosThetaVis,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    TrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueFineBinProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    TrueFineBinThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    TrueFineBinCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    TrueThetaVisPlot[0]->Fill(ThetaVis,weight);
	    TrueCosThetaVisPlot[0]->Fill(CosThetaVis,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    TrueThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    TrueCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);

	    //----------------------------------------//

	  } // End of the post-FSI selection

	  //----------------------------------------//

	  // If the signal definition pre-FSI is satisfied
	  if ( NoFSIProtonTagging == 1 && NoFSIChargedPionTagging == 0 && NoFSINeutralPionTagging == 0 && NoFSIMuonTagging == 1 && NoFSITrueHeavierMesonCounter == 0 && NoFSIElectronTagging == 0 && NoFSIPhotonTagging == 0) { 

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
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
            double ECal = reco_stv_tool.ReturnECalMB();
            double ThetaVis = reco_stv_tool.ReturnThetaVis(); // deg
	    double CosThetaVis = TMath::Cos(ThetaVis * TMath::Pi() / 180.);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    NoFSITrueFineBinThetaVisPlot[0]->Fill(ThetaVis,weight);
	    NoFSITrueFineBinCosThetaVisPlot[0]->Fill(CosThetaVis,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinProtonCosThetaPlot[genie_mode]->Fill(ProtonCosTheta,weight);
	    NoFSITrueFineBinThetaVisPlot[genie_mode]->Fill(ThetaVis,weight);
	    NoFSITrueFineBinCosThetaVisPlot[genie_mode]->Fill(CosThetaVis,weight);

	    // 1D Nominal Binning

	    NoFSITrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueProtonCosThetaPlot[0]->Fill(ProtonCosTheta,weight);
	    NoFSITrueThetaVisPlot[0]->Fill(ThetaVis,weight);
	    NoFSITrueCosThetaVisPlot[0]->Fill(CosThetaVis,weight);

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

	        // Pre FSI

	        // 1D Fine Binning

		Reweight(NoFSITrueFineBinMuonCosThetaPlot[inte]);
		Reweight(NoFSITrueFineBinProtonCosThetaPlot[inte]);
		Reweight(NoFSITrueFineBinThetaVisPlot[inte]);
		Reweight(NoFSITrueFineBinCosThetaVisPlot[inte]);

		//----------------------------------------//

	} // End of the loop over the interaction processes		

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
