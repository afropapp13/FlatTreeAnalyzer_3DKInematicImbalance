#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>

#include "constants.h"
#include "STV_Tools.h"
#include "Tools.h"

using namespace std;
using namespace constants;

//Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h, double SF = 1.);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//----------------------------------------//	

	Tools tools;

	//----------------------------------------//	

        // Output root & txt files

	// Root

	TString FileNameAndPath = "OutputFiles/" + fweights + "FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;

	//----------------------------------------//

	// Plot declaration

	// Post FSI

	// 1D Nominal Binning

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueThetaVisPlot[NInte];
	TH1D* TrueEnuPlot[NInte];

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  // Post FSI

	  // 1D Nominal Binning

	  TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",20,-1.,1.);
	  TrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaVisPlot",";#theta_{vis} [deg]",20,0.,180.);
	  TrueEnuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueEnuPlot",";E_{#nu} [GeV]",20,0.,1.5);
	  
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
	  if (fOutputFile == "GENIE_v3_4_0_G18_10a_02_11a_RS") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_10a_02_11a_RS.root"); }
	  if (fOutputFile == "GENIE_v3_4_0_G18_02a_02_11a") { fweights_file = TFile::Open("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/myWeights_uB_Tune_v3_4_0_G18_02a_02_11a.root"); }

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
	  if ( fOutputFile.Contains("ACHILLES") || fOutputFile.Contains("achilles")) { weight = fScaleFactor*Units*Weight * A; }
	  //if (jentry%1000 == 0) { std::cout << "fScaleFactor = " << fScaleFactor << ", weight = " << weight << std::endl; }

	  //	  cout << "fScaleFactor = " << fScaleFactor << ", Weight = " << Weight << endl;

	  //double weight = 1.;	
	  if (fOutputFile == "GiBUU_2021") { weight = weight/100.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2021_Inclusive") { weight = weight/10.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2023_ME_DOW") { weight = weight/5.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2023") { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 50 samples	
	  if (fOutputFile == "GiBUU_2023_Patch1") { weight = weight/500.; } // To increase the stats, the GiBUU sample has been produced in 500 samples	
	  if (fOutputFile == "GiBUU_2021_NoFSI") { weight = weight/100.; } // To increase the stats, Ben Bogart produced the no fsi GiBUU sample in 100 samples

	  // GiBUU Ben-Ulrich (BU) in medium investigation
	  if (fOutputFile == "GiBUU_2023_BU") { weight = weight/105.; }
	  if (fOutputFile == "GiBUU_2023_BU_NoFSI") { weight = weight/1119.; }
	  if (fOutputFile == "GiBUU_2023_BU_flagScreen") { weight = weight/76.; }
	  if (fOutputFile == "GiBUU_2023_BU_flagInMedium") { weight = weight/73.; }

	  //if (fOutputFile == "ACHILLES") { weight = weight*1000./(40./12.); } // ACHILLES scaling still under discussion

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state

	  // ACHILLES doesn't know how to handle the cc branch yet
	  if ( !( fOutputFile.Contains("ACHILLES") || fOutputFile.Contains("achilles") ) ) {
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

	  // Classify the events based on the interaction type

	  // https://arxiv.org/pdf/2106.15809.pdf

	  int genie_mode = -1.;

	  if (fOutputFile.Contains("ACHILLES") || fOutputFile.Contains("achilles")) {

	    if (TMath::Abs(Mode) == 200) { genie_mode = 1; } // QE
	    else if (TMath::Abs(Mode) == 700) { genie_mode = 2; } // MEC
	    else if (TMath::Abs(Mode) == 400) { genie_mode = 3; } // RES
	    else { continue; }

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

	    STV_Tools reco_stv_tool(Muon4Vector.Vect(),Proton4Vector.Vect(),Muon4Vector.E(),TMath::Sqrt( TMath::Power(Proton4Vector.Rho(),2.) + TMath::Power(ProtonMass_GeV,2.) ) );

	    double MuonCosTheta = Muon4Vector.CosTheta();

		TVector3 b = Muon4Vector.Vect()+ Proton4Vector.Vect();
		double thetavis = b.Theta() * 180./TMath::Pi();

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueEnuPlot[0]->Fill(Enu_true,weight);
		TrueThetaVisPlot[0]->Fill(thetavis,weight);

	    //----------------------------------------//

	    // filling in the histo based on the interaction mode

	    TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
		TrueEnuPlot[genie_mode]->Fill(Enu_true,weight);
		TrueThetaVisPlot[genie_mode]->Fill(thetavis,weight);

	    //----------------------------------------//

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

		Reweight(TrueMuonCosThetaPlot[inte]);	 
		Reweight(TrueEnuPlot[inte]);	 
		Reweight(TrueThetaVisPlot[inte]);	       

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

void Reweight(TH1D* h, double SF) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = SF * CurrentEntry / h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = SF * CurrentError / h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}

//----------------------------------------//		
