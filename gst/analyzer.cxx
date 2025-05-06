#define analyzer_cxx
#include "analyzer.h"

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
#include "helper_functions.cxx"
#include "Tools.cxx"

using namespace std;
using namespace constants;

//----------------------------------------//

void analyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//----------------------------------------//	

	Tools tools;

	//----------------------------------------//	

	TString FileNameAndPath = "output_files/" + fweights + "analyzer_output_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;

	//----------------------------------------//

	// Plot declaration

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueThetaVisPlot[NInte];
	TH1D* TrueThetaNuPlot[NInte];
	TH1D* TrueEnuPlot[NInte];
	TH2D* TrueTheta_vsThetaVisPlot[NInte];

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

		TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",20,-1.,1.);
		TrueThetaVisPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaVisPlot",";#theta_{vis} [deg]",20,0.,180.);
		TrueThetaNuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueThetaNuPlot",";#theta_{#nu_{#mu}} [deg]",20,0.,180.);
		TrueEnuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueEnuPlot",";E_{#nu} [GeV]",20,0.,1.5);
		TrueTheta_vsThetaVisPlot[inte] = new TH2D(InteractionLabels[inte]+"TrueTheta_vsThetaVisPlot",";#theta_{#nu_{#mu}} [deg];#theta_{vis} [deg]",20,0.,180.,20,0.,180.);

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

		//weights
		/*double weight = fScaleFactor*Units*A*Weight;	

		//----------------------------------------//	

		// Signal definition

		if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
		if (cc != 1) { continue; } // make sure that we have only CC interactions		

		int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, TrueHeavierMesonCounter = 0;
		int ElectronTagging = 0, PhotonTagging = 0;
		vector <int> ProtonID; ProtonID.clear();
		vector <int> MuonID; MuonID.clear();		

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

		//----------------------------------------//

		// before the neutrino-nucleus interaction
*/
		double numu_theta = -999.;
		/*for (int i = 0; i < ninitp; i++) {

			if (pdg_init[i] == 14) {*/

					TVector3 numu(pxv,pyv,pzv);
					numu_theta = TMath::ACos(numu.CosTheta()) * 180./TMath::Pi();
					cout << numu_theta << endl;

			/*} 
	
		}*/

	/*	//----------------------------------------//	

		// If the signal definition post-FSI  is satisfied
		if ( ProtonTagging == 1 && ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging == 1 && TrueHeavierMesonCounter == 0) { 

			// Kinematics of muon & proton in the final state

			TLorentzVector Muon4Vector(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]], E[MuonID[0]]);
			TLorentzVector Proton4Vector(px[ProtonID[0]], py[ProtonID[0]], pz[ProtonID[0]], E[ProtonID[0]]);

			//----------------------------------------//

			double MuonCosTheta = Muon4Vector.CosTheta();

			TVector3 b = Muon4Vector.Vect()+ Proton4Vector.Vect();
			double thetavis = b.Theta() * 180./TMath::Pi();

			//----------------------------------------//

			// filling in the histo regardless of interaction mode

			TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
			TrueEnuPlot[0]->Fill(Enu_true,weight);
			TrueThetaVisPlot[0]->Fill(thetavis,weight);
			TrueThetaNuPlot[0]->Fill(numu_theta,weight);

			TrueTheta_vsThetaVisPlot[0]->Fill(numu_theta,thetavis,weight);

			//----------------------------------------//

			// filling in the histo based on the interaction mode

			TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
			TrueEnuPlot[genie_mode]->Fill(Enu_true,weight);
			TrueThetaVisPlot[genie_mode]->Fill(thetavis,weight);
			TrueThetaNuPlot[genie_mode]->Fill(numu_theta,weight);

			//cout << "numu_theta = " << numu_theta << " thetavis = " << thetavis << endl;
			TrueTheta_vsThetaVisPlot[genie_mode]->Fill(thetavis,numu_theta,weight);

			//----------------------------------------//

		} // End of the post-FSI selection

		//----------------------------------------//
	*/
	} // End of the loop over the events

	//----------------------------------------//	

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes

	/*for (int inte = 0; inte < NInte; inte++) {

		divide_bin_width(TrueMuonCosThetaPlot[inte]);	 
		divide_bin_width(TrueEnuPlot[inte]);	 
		divide_bin_width(TrueThetaVisPlot[inte]);

	} // End of the loop over the interaction processes		

	//----------------------------------------//		
		*/
	file->cd();
	file->Write();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created" << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	//----------------------------------------//		

} // End of the program