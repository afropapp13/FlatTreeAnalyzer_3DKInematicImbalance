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

#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"
#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.h"
#include "/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.h"

using namespace std;
using namespace Constants;

//Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h, double SF = 1.);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

	//----------------------------------------//	

	Tools tools;

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 40.; // so that we can have xsecs per nucleus

	int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

	//----------------------------------------//	

        // Output root & txt files

	// Root

	TString FileNameAndPath = "OutputFiles/" + fweights + "FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	// Plot declaration

	TH1D* InclusiveMuonCosThetaSingleBinPlot[NInte];	
	TH1D* MuonCosThetaSingleBinPlot[NInte];
	TH1D* InclusiveMuonCosThetaPlot[NInte];	
	TH1D* InclusiveMuonEnergyPlot[NInte];	
	TH1D* InclusiveNeutronMultiPlot[NInte];	
	
	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  InclusiveMuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"InclusiveMuonCosThetaSingleBinPlot",";",1,0.,1.);	
	  MuonCosThetaSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"MuonCosThetaSingleBinPlot",";",1,0.,1.);
	
	  InclusiveMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"InclusiveMuonCosThetaPlot",";cos(#theta_{#mu})",20,-1.,1.);	
	  InclusiveMuonEnergyPlot[inte] = new TH1D(InteractionLabels[inte]+"InclusiveMuonEnergyPlot",";E_{#mu} [GeV]",20,0.1,1.3);	
	  InclusiveNeutronMultiPlot[inte] = new TH1D(InteractionLabels[inte]+"InclusiveNeutronMultiPlot",";#neutrons",6,-0.5,5.5);	
	  
	  //--------------------------------------------------//

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;

	//----------------------------------------//
	
	// Loop over the events
	cout << "nentries = " << nentries << endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	
	
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//			

	  double weight = fScaleFactor*Units*A*Weight;	
	  if (fOutputFile == "GiBUU_2021") { weight = weight/100.; } // To increase the stats, the GiBUU sample has been produced in 50 samples

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  // ACHILLES doesn't know how to handle the cc branch yet
	  if (fOutputFile != "ACHILLES") {
	    if (cc != 1) { continue; } // make sure that we have only CC interactions		
	  }

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, TrueHeavierMesonCounter = 0, NeutronCounter = 0;
	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();		
	  vector <int> NeutronID; NeutronID.clear();		

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

	    if (pdg[i] == 2112 && (pf > 0.5 && pf < 5.) ) {

	      NeutronCounter ++;
	      NeutronID.push_back(i);

	    }


	    if (fabs(pdg[i]) == 211 && pf > 0.07)  {

	      ChargedPionTagging ++;

	    }

	    if (pdg[i] == 111)  {

	      NeutralPionTagging ++;

	    }

	    if ( pdg[i] != NeutralPionPdg && fabs(pdg[i]) != AbsChargedPionPdg && tools.is_meson_or_antimeson(pdg[i]) ) { TrueHeavierMesonCounter++; }


	  } // End of the loop over the final state particles / post FSI

	  //----------------------------------------//	

	  // Classify the events based on the interaction type

	  // https://arxiv.org/pdf/2106.15809.pdf

	  int genie_mode = -1.;

	  if (TMath::Abs(Mode) == 1) { genie_mode = 1; } // QE
	  else if (TMath::Abs(Mode) == 2) { genie_mode = 2; } // MEC
	  else if (
		   TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
		   TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
		   ) { genie_mode = 3; } // RES
	  else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { genie_mode = 4; } // DIS
	  else if (TMath::Abs(Mode) == 16) { genie_mode = 5;} // COH
	  else { continue; }  

	  //----------------------------------------//	

	  // If the numu CC inclusive signal definition
	  // is not satisfied, continue

	  if (MuonTagging != 1 || ChargedPionTagging != 0 || NeutralPionTagging != 0) { continue; }

          TLorentzVector Muon4Vector(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]], E[MuonID[0]]);
	  double Emu = E[MuonID[0]];
	  //if (E[MuonID[0]] < 0.4) { continue; }
	  double costheta_mu = Muon4Vector.CosTheta();
	  //if (costheta_mu > 0.4) { continue; }

	  //filling in the histo regardless of interaction mode

	  InclusiveMuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);
	  InclusiveMuonCosThetaPlot[0]->Fill(costheta_mu,weight);
	  InclusiveMuonEnergyPlot[0]->Fill(Emu,weight);
	  InclusiveNeutronMultiPlot[0]->Fill(NeutronID.size(),weight);
	
	  //filling in the histo based on the interaction mode

	  InclusiveMuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);
	  InclusiveMuonCosThetaPlot[genie_mode]->Fill(costheta_mu,weight);
	  InclusiveMuonEnergyPlot[genie_mode]->Fill(Emu,weight);
	  InclusiveNeutronMultiPlot[genie_mode]->Fill(NeutronID.size(),weight);

	  //----------------------------------------//

	  // If the CC1p0pi signal definition post-FSI  is satisfied

	  if ( ProtonTagging == 1 && ChargedPionTagging == 0 && NeutralPionTagging == 0 && TrueHeavierMesonCounter == 0) { 

	    CounterEventsPassedSelection++;

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    MuonCosThetaSingleBinPlot[0]->Fill(0.5,weight);

	    // filling in the histo based on the interaction mode

	    MuonCosThetaSingleBinPlot[genie_mode]->Fill(0.5,weight);

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

		//----------------------------------------//
	
	        // 1D Fine Binning Post FSI

		Reweight(InclusiveMuonCosThetaSingleBinPlot[inte]);
		Reweight(MuonCosThetaSingleBinPlot[inte]);
		Reweight(InclusiveMuonCosThetaPlot[inte]);
		Reweight(InclusiveMuonEnergyPlot[inte]);

		//----------------------------------------//

	} // End of the loop over the interaction processes		

	file->cd();
	file->Write();
	fFile->Close();

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
