#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

namespace constants {

	//----------------------------------------//

    // Kerberos user name
  
	TString UserID = getenv("USER");

	//----------------------------------------//

	// Argon 

	static const double A = 40.;
	static const double Z = 18.;

	const int FontStyle = 132;
	const double TextSize = 0.07;
	const int NCont = 999; 

	//----------------------------------------//

	// Labels / Ranges & Label  map
	// max values

	static std::map<TString,std::pair<double,double> > XSecRange =
	{
		{ "MuonCosThetaPlot",  std::make_pair(0, 24) },
		{ "ThetaVisPlot",  std::make_pair(0, 14) },
		{ "EnuPlot",  std::make_pair(0, 14) },						
								
	};	
	
	//----------------------------------------//

	static std::map<TString,TString> VarLabel =
	{
		{ "MuonCosThetaPlot",  "#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "ThetaVisPlot",  "#frac{d#sigma}{d#theta_{vis}} #left[10^{-38} #frac{cm^{2}}{Ar deg}#right]" },
		{ "EnuPlot",  "#frac{d#sigma}{dE_{#nu}} #left[10^{-38} #frac{cm^{2}}{Ar GeV}#right]" },
					
	};	
	
	static std::map<TString,TString> LatexLabel =
	{

		{ "MuonCosThetaPlot",  "all events" },
		{ "ThetaVisPlot",  "all events" },
		{ "EnuPlot",  "all events" },

	};	

	// -------------------------------------------------------------------------------------------------------------------------			

	// Global Constants

	static const double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
//	static const double NTargets = 1.203E30; // Argon nuclei, not nucleons
	static const double NTargets = 1.05E30; // Argon nuclei, not nucleons
	
	static const int NuMuPdg = 14, MuonPdg = 13, ProtonPdg = 2212, AbsChargedPionPdg = 211, NeutralPionPdg = 111;
	static const int ElectronPdg = 11, PhotonPdg = 22, NeutronPdg = 2112, KaonPdg = 321;
	static const int DeuteriumPdg = 1000010020, HeliumPdg = 1000020040, ArgonPdg = 1000180400;

	static const double MuonMass = 106, ProtonMass = 938.272, NeutronMass = 939.565; // MeV
	static const double MuonMass_GeV = 0.106, ProtonMass_GeV = 0.938272, NeutronMass_GeV = 0.939565; // GeV
	static const double DeltaM2 = TMath::Power(NeutronMass_GeV,2.) - TMath::Power(ProtonMass_GeV,2.); // GeV^2	
	
	//----------------------------------------//
	
	// Interaction labels
	
	//const std::vector<int> InteBreakColors{kBlack,kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};
	const std::vector<int> InteBreakColors{kBlack,kAzure-4,kOrange-3,kGreen+1,kRed+1,kBlue};		
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();
	
	//----------------------------------------//		

}
#endif