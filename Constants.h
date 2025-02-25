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

namespace Constants {

	//----------------------------------------//

        // Kerberos user name
  
	TString UserID = "apapadop";

	//----------------------------------------//

	// Argon 

	static const double A = 40.;
	static const double Z = 18.;

	const int FontStyle = 132;
	const double TextSize = 0.07;
	const int NCont = 999; 
	const TString CutExtension = "_NoCuts_PID_NuScore";

	// UBCodeVersion

	const TString UBCodeVersion = "v08_00_00_52";

	//----------------------------------------//

	// Exporting paths

	const TString PathToFiles = "/uboone/data/users/"+UserID+"/myEvents/OutputFiles/"+UBCodeVersion+"/";
	const TString PathToExtractedXSec = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myXSec/"+UBCodeVersion+"/"; 
	const TString FileEfficienciesPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myEfficiencies/"+UBCodeVersion+"/";
	const TString PlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myPlots/"+UBCodeVersion+"/"; 
	const TString CCQEPlotPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myCCQEPlots/"+UBCodeVersion+"/"; 
	const TString PathToSystematics = "/uboone/data/users/"+UserID+"/mySTVAnalysis/mySystematics/"+UBCodeVersion+"/";
	const TString TrueSTVPath = PathToFiles;
	const TString MigrationMatrixPath = "/uboone/data/users/"+UserID+"/mySTVAnalysis/myMigrationMatrices/"+UBCodeVersion+"/"; 

	//----------------------------------------//

	// Zarko Pavlovic, Jun 22 2020
	// /pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist/readme.txt

	double Nominal_UB_XY_Surface = 256.35*233.; // cm2
	double SoftFidSurface = 236. * 210.;  // cm2
	double POTPerSpill = 4997.*5e8;

	//----------------------------------------//

	std::vector<double> range{0.,0.1,0.2,0.3,10.};
	int NRanges = range.size() - 1; 
	int NBins = 80;
	double LowEdge = 0.;
	double HighEdge = 0.8;
	
	std::vector<TString> RangeLabel{"","0 < p_{n} < 0.1 GeV/c","0.1 < p_{n} < 0.2 GeV/c","0.2 < p_{n} < 0.3 GeV/c","p_{n} > 0.3 GeV/c"};

	//----------------------------------------//

	// Labels / Ranges & Label  map
	// max values

	static std::map<TString,std::pair<double,double> > XSecRange =
	{
		{ "DeltaPTPlot",  std::make_pair(0, 59) }, // 57
		{ "DeltaAlphaTPlot",  std::make_pair(0, 0.12) },
		{ "DeltaAlpha3DqPlot",  std::make_pair(0, 0.18) },
		{ "DeltaAlpha3DMuPlot",  std::make_pair(0, 0.16) },				
		{ "DeltaPhiTPlot",  std::make_pair(0., 0.39) },
		{ "DeltaPhi3DPlot",  std::make_pair(0., 0.39) },		
		{ "MuonMomentumPlot",  std::make_pair(0, 29) },
		{ "MuonCosThetaPlot",  std::make_pair(0, 24) },
		{ "MuonCosThetaSingleBinPlot",  std::make_pair(0, 14) },
		{ "MuonPhiPlot",  std::make_pair(0, 0.049) },
		{ "ProtonMomentumPlot",  std::make_pair(0, 34) },
		{ "ProtonCosThetaPlot",  std::make_pair(0, 24) },
		{ "ProtonPhiPlot",  std::make_pair(0, 0.045) },
		{ "ECalPlot",  std::make_pair(0, 29.) },			
		{ "EQEPlot",  std::make_pair(0, 29.) },
		{ "Q2Plot",  std::make_pair(0, 29) },	
		{ "CCQEMuonMomentumPlot",  std::make_pair(0, 24) },
		{ "CCQEMuonCosThetaPlot",  std::make_pair(0, 9) },
		{ "CCQEProtonMomentumPlot",  std::make_pair(0, 37) },
		{ "CCQEProtonCosThetaPlot",  std::make_pair(0, 11) },
		{ "CCQEECalPlot",  std::make_pair(0, 24) },
		{ "CCQEQ2Plot",  std::make_pair(0, 24) },	
		{ "DeltaPLPlot",  std::make_pair(0, 47) },
		{ "DeltaPnPlot",  std::make_pair(0, 49) },
		{ "DeltaPnPerpPlot",  std::make_pair(0, 64) },
		{ "DeltaPnPerpxPlot",  std::make_pair(0, 49) },
		{ "DeltaPnPerpyPlot",  std::make_pair(0, 54) },				
		{ "DeltaPnParPlot",  std::make_pair(0, 39) },				
		{ "DeltaPtxPlot",  std::make_pair(0, 54) },
		{ "DeltaPtyPlot",  std::make_pair(0, 39) },	
		{ "APlot",  std::make_pair(0, 34) },	
		{ "kMissPlot",  std::make_pair(0, 49) },
		{ "PMissPlot",  std::make_pair(0, 49) },	
		{ "PMissMinusPlot",  std::make_pair(0, 39) },
		{ "ECal_DeltaPT_0_00To0_20Plot",  std::make_pair(0, 63.9) },	
		{ "ECal_DeltaPT_0_20To0_40Plot",  std::make_pair(0, 36.9) },	
		{ "ECal_DeltaPT_0_40To1_00Plot",  std::make_pair(0, 7.9) },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot",  std::make_pair(0, 10.9) },
		{ "ECal_DeltaPtx_Minus0_15To0_15Plot",  std::make_pair(0, 46.9) },
		{ "ECal_DeltaPtx_0_15To0_55Plot",  std::make_pair(0, 11.9) },
		{ "ECal_DeltaPty_Minus0_75ToMinus0_15Plot",  std::make_pair(0, 12.9) },
		{ "ECal_DeltaPty_Minus0_15To0_15Plot",  std::make_pair(0, 46.9) },
		{ "ECal_DeltaPty_0_15To0_45Plot",  std::make_pair(0, 6.9) },											
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot",  std::make_pair(0, 0.37) },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot",  std::make_pair(0, 0.19) },				
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot",  std::make_pair(0, 0.059) },
		{ "DeltaAlphaT_MuonMomentum_0_10To0_40Plot",  std::make_pair(0, 0.19) },
		{ "DeltaAlphaT_MuonMomentum_0_40To0_60Plot",  std::make_pair(0, 0.21) },				
		{ "DeltaAlphaT_MuonMomentum_0_60To1_20Plot",  std::make_pair(0, 0.089) },
		{ "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot",  std::make_pair(0, 0.39) },
		{ "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot",  std::make_pair(0, 0.23) },				
		{ "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot",  std::make_pair(0, 0.079) },				
		{ "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot",  std::make_pair(0, 0.32) },
		{ "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot",  std::make_pair(0, 0.44) },				
		{ "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot",  std::make_pair(0, 0.149) },
		{ "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot",  std::make_pair(0, 0.34) },
		{ "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot",  std::make_pair(0, 0.39) },				
		{ "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot",  std::make_pair(0, 0.099) },				
		{ "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 0.02) },
		{ "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot",  std::make_pair(0, 0.06) },
		{ "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot",  std::make_pair(0, 0.12) },
		{ "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot",  std::make_pair(0, 0.21) },
		{ "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 0.012) },
		{ "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot",  std::make_pair(0, 0.07) },
		{ "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot",  std::make_pair(0, 0.149) },
		{ "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot",  std::make_pair(0, 0.24) },		
		{ "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 9.) },
		{ "DeltaPT_MuonCosTheta_0_00To0_50Plot",  std::make_pair(0, 22.) },
		{ "DeltaPT_MuonCosTheta_0_50To0_75Plot",  std::make_pair(0, 48.) },
		{ "DeltaPT_MuonCosTheta_0_75To1_00Plot",  std::make_pair(0, 64.) },
		{ "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 3.9) },
		{ "DeltaPT_ProtonCosTheta_0_00To0_50Plot",  std::make_pair(0, 22.) },
		{ "DeltaPT_ProtonCosTheta_0_50To0_75Plot",  std::make_pair(0, 64.) },
		{ "DeltaPT_ProtonCosTheta_0_75To1_00Plot",  std::make_pair(0, 74.) },
		{ "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.19) },
		{ "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.19) },
		{ "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.21) },
		{ "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.27) },
		{ "DeltaPT_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.34) },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.27) },
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.28) },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.34) },
		{ "DeltaPn_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.24) },
		{ "DeltaPn_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.24) },
		{ "DeltaPn_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.28) },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.26) },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.34) },
		{ "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot",  std::make_pair(0, 0.1499) },
		{ "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot",  std::make_pair(0, 0.34) },
		{ "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot",  std::make_pair(0, 0.39) },
		{ "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot",  std::make_pair(0, 0.229) },
		{ "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot",  std::make_pair(0, 0.159) },
		{ "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot",  std::make_pair(0, 0.34) },
		{ "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot",  std::make_pair(0, 0.42) },
		{ "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot",  std::make_pair(0, 0.199) },										
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 11.) },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot",  std::make_pair(0, 19.) },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot",  std::make_pair(0, 46.) },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot",  std::make_pair(0, 36.) },
		{ "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 9.9) },
		{ "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot",  std::make_pair(0, 13.) },
		{ "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot",  std::make_pair(0, 27.) },
		{ "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot",  std::make_pair(0, 29.) },		
		{ "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot",  std::make_pair(0, 8.) },
		{ "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot",  std::make_pair(0, 21.) },
		{ "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot",  std::make_pair(0, 52.) },
		{ "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot",  std::make_pair(0, 59.) },	
		{ "DeltaPhiT_DeltaPT_0_00To0_20Plot",  std::make_pair(0, 1.8) },
		{ "DeltaPhiT_DeltaPT_0_20To0_40Plot",  std::make_pair(0, 0.35) },				
		{ "DeltaPhiT_DeltaPT_0_40To1_00Plot",  std::make_pair(0, 0.049) },	
//		{ "DeltaPn_DeltaPT_0_00To0_20Plot",  std::make_pair(0, 219.) },
//		{ "DeltaPn_DeltaPT_0_20To0_40Plot",  std::make_pair(0, 95.) },				
//		{ "DeltaPn_DeltaPT_0_40To1_00Plot",  std::make_pair(0, 19.) },
		{ "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot",  std::make_pair(0, 14.) },
		{ "DeltaPty_DeltaPtx_Minus0_15To0_15Plot",  std::make_pair(0, 89.) },				
		{ "DeltaPty_DeltaPtx_0_15To0_55Plot",  std::make_pair(0, 13.) },
		{ "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot",  std::make_pair(0, 15.4) },
		{ "DeltaPtx_DeltaPty_Minus0_15To0_15Plot",  std::make_pair(0, 99.) },
		{ "DeltaPtx_DeltaPty_0_15To0_45Plot",  std::make_pair(0, 14.) },		
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.42) },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.44) },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.47) },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.43) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.14) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.16) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.29) },	
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.29) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.006) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.0159) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.036) },	
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.11) },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot",  std::make_pair(0, 10.4) },		
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot",  std::make_pair(0, 16.4) },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot",  std::make_pair(0, 4.2) },	
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot",  std::make_pair(0, 24.) },		
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot",  std::make_pair(0, 129.) },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot",  std::make_pair(0, 16.) },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot",  std::make_pair(0, 9.) },		
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot",  std::make_pair(0, 16.) },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot",  std::make_pair(0, 2.9) },							
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot",  std::make_pair(0, 14.9) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot",  std::make_pair(0, 24.) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot",  std::make_pair(0, 19.) },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot",  std::make_pair(0, 7.) },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot",  std::make_pair(0, 2.2) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot",  std::make_pair(0, 34.) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot",  std::make_pair(0, 69.) },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot",  std::make_pair(0, 21.) },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot",  std::make_pair(0, 0.03) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot",  std::make_pair(0, 1.6) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot",  std::make_pair(0, 24.) },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot",  std::make_pair(0, 39.) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot",  std::make_pair(0, 9.9) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot",  std::make_pair(0, 40.) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot",  std::make_pair(0, 54.) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot",  std::make_pair(0, 39.) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot",  std::make_pair(0, 2.4) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot",  std::make_pair(0, 24.) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot",  std::make_pair(0, 79.) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot",  std::make_pair(0, 60.) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot",  std::make_pair(0, 0.22) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot",  std::make_pair(0, 3.5) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot",  std::make_pair(0, 29.) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot",  std::make_pair(0, 64.) },
		{ "ECal_DeltaAlphaT_0_00To45_00Plot",  std::make_pair(0, 0.099) },
		{ "ECal_DeltaAlphaT_45_00To90_00Plot",  std::make_pair(0, 0.11) },
		{ "ECal_DeltaAlphaT_90_00To135_00Plot",  std::make_pair(0, 0.149) },	
		{ "ECal_DeltaAlphaT_135_00To180_00Plot",  std::make_pair(0, 0.18) },		
		{ "SerialDeltaPT_MuonCosThetaPlot",  std::make_pair(0, 0.8) },
		{ "SerialDeltaPT_ProtonCosThetaPlot",  std::make_pair(0, 1.1) },
		{ "SerialMuonMomentum_MuonCosThetaPlot",  std::make_pair(0, 1.59) },
		{ "SerialProtonMomentum_ProtonCosThetaPlot",  std::make_pair(0, 1.34) },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot",  std::make_pair(0, 1.29) },							
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot",  std::make_pair(0, 1.73) },
		{ "SerialDeltaAlphaT_DeltaPTPlot",  std::make_pair(0, 1.83) },
		{ "SerialDeltaAlphaT_MuonMomentumPlot",  std::make_pair(0, 1.73) },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot",  std::make_pair(0, 1.73) },				
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot",  std::make_pair(0, 1.79) },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot",  std::make_pair(0, 1.79) },		
		{ "SerialProtonMomentum_DeltaAlphaTPlot",  std::make_pair(0, 1.2) },
		{ "SerialDeltaPT_DeltaAlphaTPlot",  std::make_pair(0, 1.2) },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot",  std::make_pair(0, 1.9) },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot",  std::make_pair(0, 1.3) },		
		{ "SerialDeltaPn_DeltaAlphaTPlot",  std::make_pair(0, 1.39) },		
		{ "SerialDeltaPhiT_DeltaPTPlot",  std::make_pair(0, 4.09) },
//		{ "SerialDeltaPn_DeltaPTPlot",  std::make_pair(0, 2.5) },
		{ "SerialECal_DeltaPTPlot",  std::make_pair(0, 1.6) },
		{ "SerialECal_DeltaPtxPlot",  std::make_pair(0, 1.6) },
		{ "SerialECal_DeltaPtyPlot",  std::make_pair(0, 1.6) },				
		{ "SerialECal_DeltaAlphaTPlot",  std::make_pair(0, 1.4) },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot",  std::make_pair(0, 0.9) },
		{ "SerialECal_DeltaPtxDeltaPtyPlot",  std::make_pair(0, 1.9) },					
		{ "SerialECal_MuonCosThetaMuonMomentumPlot",  std::make_pair(0, 2.4) },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot",  std::make_pair(0, 1.1) },		
		{ "SerialProtonCosTheta_MuonCosThetaPlot",  std::make_pair(0, 1.39) },
		{ "SerialDeltaPty_DeltaPtxPlot",  std::make_pair(0, 2.9) },
		{ "SerialDeltaPtx_DeltaPtyPlot",  std::make_pair(0, 2.4) },	
		{ "VertexXPlot",  std::make_pair(0, 0.079) },
		{ "VertexYPlot",  std::make_pair(0, 0.079) },
		{ "VertexZPlot",  std::make_pair(0, 0.014) },							
								
	};	
	
	//----------------------------------------//

	static std::map<TString,TString> VarLabel =
	{
		{ "DeltaPTPlot",  "#frac{d#sigma}{d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaAlphaTPlot",  "#frac{d#sigma}{d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlpha3DqPlot",  "#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlpha3DMuPlot",  "#frac{d#sigma}{d#alpha_{3D,#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },		
		{ "DeltaPhiTPlot",  "#frac{d#sigma}{d#delta#phi_{T}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaPhi3DPlot",  "#frac{d#sigma}{d#phi_{3D}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },		
		{ "MuonMomentumPlot",  "#frac{d#sigma}{dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "MuonCosThetaPlot",  "#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "MuonCosThetaSingleBinPlot",  "#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "MuonPhiPlot",  "#frac{d#sigma}{d#phi_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "ProtonMomentumPlot",  "#frac{d#sigma}{dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "ProtonCosThetaPlot",  "#frac{d#sigma}{dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "ProtonPhiPlot",  "#frac{d#sigma}{d#phi_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "ECalPlot",  "#frac{d#sigma}{dE^{Cal}} #left[10^{-38} #frac{cm^{2}}{GeV Ar}#right]" },
		{ "ECal_DeltaPT_0_00To0_20Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },	
		{ "ECal_DeltaPT_0_20To0_40Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_40To1_00Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_15To0_15Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPtx_0_15To0_55Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },	
		{ "ECal_DeltaPty_Minus0_75ToMinus0_15Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPty_Minus0_15To0_15Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPty_0_15To0_45Plot",  "#frac{d#sigma}{dE^{Cal}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },							
		{ "EQEPlot",  "#frac{d#sigma}{dE^{QE}} #left[10^{-38} #frac{cm^{2}}{GeV Ar}#right]" },
		{ "Q2Plot",  "#frac{d#sigma}{dQ^{2}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },	
		{ "CCQEMuonMomentumPlot",  "#frac{d#sigma}{dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "CCQEMuonCosThetaPlot",  "#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "CCQEProtonMomentumPlot",  "#frac{d#sigma}{dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "CCQEProtonCosThetaPlot",  "#frac{d#sigma}{dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "CCQEECalPlot",  "#frac{d#sigma}{dE^{Cal}} #left[10^{-38} #frac{cm^{2}}{GeV Ar}#right]" },
		{ "CCQEQ2Plot",  "#frac{d#sigma}{dQ^{2}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },	
		{ "DeltaPLPlot",  "#frac{d#sigma}{d#deltap_{L}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPnPlot",  "#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPnPerpPlot",  "#frac{d#sigma}{dp_{n#perp} }   #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPnPerpxPlot",  "#frac{d#sigma}{dp_{n#perp ,x} }   #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPnPerpyPlot",  "#frac{d#sigma}{dp_{n#perp ,y} }   #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },				
		{ "DeltaPnParPlot",  "#frac{d#sigma}{dp_{n#parallel}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },				
		{ "DeltaPtxPlot",  "#frac{d#sigma}{d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPtyPlot",  "#frac{d#sigma}{d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "APlot",  "#frac{d#sigma}{d#alpha} #left[10^{-38} #frac{cm^{2}}{c Ar}#right]" },
		{ "kMissPlot",  "#frac{d#sigma}{dk_{Miss}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "PMissPlot",  "#frac{d#sigma}{dP_{Miss}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "PMissMinusPlot",  "#frac{d#sigma}{dp_{Miss}^{-}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlphaT_MuonMomentum_0_10To0_40Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlphaT_MuonMomentum_0_40To0_60Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlphaT_MuonMomentum_0_60To1_20Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },				
		{ "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "ProtonMomentum_DeltaAlphaT_0_00To45Plot",  "#frac{d^{2}#sigma}{d#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "ProtonMomentum_DeltaAlphaT_0_00To45Plot",  "#frac{d^{2}#sigma}{d#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPT_DeltaAlphaT_0_00To45Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaPn_DeltaAlphaT_0_00To45Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlphaT_45_00To90_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaPn_DeltaAlphaT_90_00To135_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },	
		{ "DeltaPn_DeltaAlpha3Dq_0_00To45Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },	
		{ "DeltaPn_DeltaAlpha3DMu_0_00To45Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot",  "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot",  "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },		
		{ "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_MuonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_MuonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_MuonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_ProtonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_ProtonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPT_ProtonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },		
		{ "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot",  "#frac{d^{2}#sigma}{dcos#theta_{p}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{dcos#theta_{p}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{dcos#theta_{p}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{dcos#theta_{p}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot", "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot", "#frac{d^{2}#sigma}{dp_{p}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot", "#frac{d^{2}#sigma}{dp_{p}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot", "#frac{d^{2}#sigma}{dp_{p}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot", "#frac{d^{2}#sigma}{dp_{p}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "DeltaPhiT_DeltaPT_0_00To0_20Plot",  "#frac{d^{2}#sigma}{d#delta#phi_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "DeltaPhiT_DeltaPT_0_20To0_40Plot",  "#frac{d^{2}#sigma}{d#delta#phi_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/cAr}#right]" },		
		{ "DeltaPhiT_DeltaPT_0_40To1_00Plot",  "#frac{d^{2}#sigma}{d#delta#phi_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },	
//		{ "DeltaPn_DeltaPT_0_00To0_20Plot",  "#frac{d^{2}#sigma}{dp_{n}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
//		{ "DeltaPn_DeltaPT_0_20To0_40Plot",  "#frac{d^{2}#sigma}{dp_{n}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },		
//		{ "DeltaPn_DeltaPT_0_40To1_00Plot",  "#frac{d^{2}#sigma}{dp_{n}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,y}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPty_DeltaPtx_Minus0_15To0_15Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,y}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPty_DeltaPtx_0_15To0_55Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,y}d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPtx_DeltaPty_Minus0_15To0_15Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "DeltaPtx_DeltaPty_0_15To0_45Plot",  "#frac{d^{2}#sigma}{d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot", "#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },	
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },	
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot","#frac{d^{3}#sigma}{dE^{Cal}d#deltap_{T,x}d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },	
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{#mu}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot","#frac{d^{3}#sigma}{dE^{Cal}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot","#frac{d^{3}#sigma}{dE^{Cal}dcos#theta_{p}dp_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]"},	
		{ "ECal_DeltaAlphaT_0_00To45_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "ECal_DeltaAlphaT_45_00To90_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV Ar}#right]" },
		{ "ECal_DeltaAlphaT_90_00To135_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV Ar}#right]" },	
		{ "ECal_DeltaAlphaT_135_00To180_00Plot","#frac{d^{3}#sigma}{dE^{Cal}d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV Ar}#right]" },
		{ "SerialDeltaPT_MuonCosThetaPlot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "SerialDeltaPT_ProtonCosThetaPlot", "#frac{d^{2}#sigma}{d#deltap_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "SerialMuonMomentum_MuonCosThetaPlot", "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot", "#frac{d^{2}#sigma}{dp_{p}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV/c Ar}#right]" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },				
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{deg Ar}#right]" },
		{ "SerialDeltaAlphaT_DeltaPTPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot", "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot", "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "SerialProtonMomentum_DeltaAlphaTPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaPT_DeltaAlphaTPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaPT_MuonMomentumPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{#mu}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaPT_ProtonMomentumPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{p}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },				
		{ "SerialDeltaPn_DeltaAlphaTPlot", "#frac{d^{2}#sigma}{d#delta#alpha_{T}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "SerialDeltaPn_DeltaAlpha3DqPlot", "#frac{d^{2}#sigma}{d#alpha_{3D}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot", "#frac{d^{2}#sigma}{d#alpha_{3D,#mu}dp_{n}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },		
		{ "SerialDeltaPhiT_DeltaPTPlot", "#frac{d^{2}#sigma}{d#delta#phi_{T}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV/c Ar}#right]" },
//		{ "SerialDeltaPn_DeltaPTPlot",  "#frac{d^{2}#sigma}{dp_{n}d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "SerialECal_DeltaPTPlot", "#frac{d^{2}#sigma}{dE^{Cal} d#deltap_{T}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "SerialECal_DeltaPtxPlot", "#frac{d^{2}#sigma}{dE^{Cal} d#deltap_{T,x}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "SerialECal_DeltaPtyPlot", "#frac{d^{2}#sigma}{dE^{Cal} d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },				
		{ "SerialECal_DeltaAlphaTPlot", "#frac{d^{2}#sigma}{dE^{Cal} d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV Ar}#right]" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot", "#frac{d^{3}#sigma}{dE^{Cal} d#deltap_{T} d#delta#alpha_{T}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{2}/c Ar}#right]" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot", "#frac{d^{3}#sigma}{dE^{Cal} d#deltap_{T,x} d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{deg GeV^{3}/c^{2} Ar}#right]" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot", "#frac{d^{3}#sigma}{dE^{Cal} dp_{#mu} dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot", "#frac{d^{3}#sigma}{dE^{Cal} dp_{p} dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c Ar}#right]" },		
		{ "SerialProtonCosTheta_MuonCosThetaPlot",  "#frac{d^{2}#sigma}{dcos#theta_{#mu} dcos#theta_{p}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]" },
		{ "SerialDeltaPty_DeltaPtxPlot", "#frac{d^{2}#sigma}{d#deltap_{T,x} d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },
		{ "SerialDeltaPtx_DeltaPtyPlot", "#frac{d^{2}#sigma}{d#deltap_{T,x} d#deltap_{T,y}} #left[10^{-38} #frac{cm^{2}}{GeV^{2}/c^{2} Ar}#right]" },	
		{ "VertexXPlot",  "#frac{d#sigma}{dx} #left[10^{-38} #frac{cm^{2}}{cm Ar}#right]" },
		{ "VertexYPlot",  "#frac{d#sigma}{dy} #left[10^{-38} #frac{cm^{2}}{cm Ar}#right]" },
		{ "VertexZPlot",  "#frac{d#sigma}{dz} #left[10^{-38} #frac{cm^{2}}{cm Ar}#right]" },			
					
	};	
	
	static std::map<TString,TString> LatexLabel =
	{
		{ "DeltaPTPlot",  "All events" },
		{ "DeltaAlphaTPlot",  "All events" },
		{ "DeltaAlpha3DqPlot",  "All events" },
		{ "DeltaAlpha3DMuPlot",  "All events" },		
		{ "DeltaPhiTPlot",  "All events" },
		{ "DeltaPhi3DPlot",  "All events" },		
		{ "MuonMomentumPlot",  "All events" },
		{ "MuonCosThetaPlot",  "All events" },
		{ "MuonCosThetaSingleBinPlot",  "All events" },
		{ "MuonPhiPlot",  "All events" },
		{ "ProtonMomentumPlot",  "All events" },
		{ "ProtonCosThetaPlot",  "All events" },
		{ "ProtonPhiPlot",  "All events" },
		{ "ECalPlot",  "All events" },
		{ "ECal_DeltaPT_0_00To0_20Plot",  "#deltap_{T} < 0.2 GeV/c" },	
		{ "ECal_DeltaPT_0_20To0_40Plot",  "0.2 < #deltap_{T} < 0.4 GeV/c" },
		{ "ECal_DeltaPT_0_40To1_00Plot",  "#deltap_{T} > 0.4 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot",  "#deltap_{T,x} < -0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_15To0_15Plot",  "|#deltap_{T,x}| < 0.15 GeV/c" },				
		{ "ECal_DeltaPtx_0_15To0_55Plot",  "#deltap_{T,x} > 0.15 GeV/c" },
		{ "ECal_DeltaPty_Minus0_75ToMinus0_15Plot","#deltap_{T,y} < -0.15 GeV/c" },
		{ "ECal_DeltaPty_Minus0_15To0_15Plot","|#deltap_{T,y}| < 0.15 GeV/c" },
		{ "ECal_DeltaPty_0_15To0_45Plot","#deltap_{T,y} > 0.15 GeV/c" },								
		{ "EQEPlot",  "All events" },
		{ "Q2Plot",  "All events" },	
		{ "CCQEMuonMomentumPlot",  "" },
		{ "CCQEMuonCosThetaPlot",  "" },
		{ "CCQEProtonMomentumPlot",  "" },
		{ "CCQEProtonCosThetaPlot",  "" },
		{ "CCQEECalPlot",  "" },
		{ "CCQEQ2Plot",  "" },	
		{ "DeltaPLPlot",  "All events" },
		{ "DeltaPnPlot",  "All events" },
		{ "DeltaPnPerpPlot",  "All events" },
		{ "DeltaPnPerpxPlot",  "All events" },
		{ "DeltaPnPerpyPlot",  "All events" },				
		{ "DeltaPnParPlot",  "All events" },				
		{ "DeltaPtxPlot",  "All events" },
		{ "DeltaPtyPlot",  "All events" },
		{ "APlot",  "All events" },
		{ "kMissPlot",  "All events" },
		{ "PMissPlot",  "All events" },
		{ "PMissMinusPlot",  "All events" },
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot",  "#deltap_{T} < 0.2 GeV/c" },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot",  "0.2 < #deltap_{T} < 0.4 GeV/c" },				
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot",  "#deltap_{T} > 0.4 GeV/c" },
		{ "DeltaAlphaT_MuonMomentum_0_10To0_40Plot",  "0.1 < p_{#mu} < 0.4 GeV/c" },
		{ "DeltaAlphaT_MuonMomentum_0_40To0_60Plot",  "0.4 < p_{#mu} < 0.6 GeV/c" },				
		{ "DeltaAlphaT_MuonMomentum_0_60To1_20Plot",  "0.6 < p_{#mu} < 1.2 GeV/c" },
		{ "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot",  "0.3 < p_{p} < 0.5 GeV/c" },
		{ "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot",  "0.5 < p_{p} < 0.7 GeV/c" },				
		{ "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot",  "0.7 < p_{p} < 1 GeV/c" },				
		{ "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot",  "p_{n} < 0.2 GeV/c" },
		{ "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot",  "0.2 < p_{n} < 0.4 GeV/c" },				
		{ "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot",  "p_{n} > 0.4 GeV/c" },
		{ "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot",  "p_{n} < 0.2 GeV/c" },
		{ "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot",  "0.2 < p_{n} < 0.4 GeV/c" },				
		{ "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot",  "p_{n} > 0.4 GeV/c" },			
		{ "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot",  "-1 < cos#theta_{#mu} < 0" },
		{ "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot", "0 < cos#theta_{#mu} < 0.5" },
		{ "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{#mu} < 0.75" },
		{ "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{#mu} < 1" },
		{ "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot",  "-1 < cos#theta_{p} < 0" },
		{ "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot", "0 < cos#theta_{p} < 0.5" },
		{ "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{p} < 0.75" },
		{ "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{p} < 1" },
		{ "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot",  "#delta#alpha_{T} < 45^{o}" },
		{ "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot",  "45^{o} < #delta#alpha_{T} < 90^{o}" },				
		{ "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot",  "90^{o} < #delta#alpha_{T} < 135^{o}" },
		{ "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot",  "135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "DeltaPT_DeltaAlphaT_0_00To45_00Plot",  "#delta#alpha_{T} < 45^{o}" },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot",  "45^{o} < #delta#alpha_{T} < 90^{o}" },				
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot",  "90^{o} < #delta#alpha_{T} < 135^{o}" },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot",  "135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "DeltaPn_DeltaAlphaT_0_00To45_00Plot",  "#delta#alpha_{T} < 45^{o}" },
		{ "DeltaPn_DeltaAlphaT_45_00To90_00Plot",  "45^{o} < #delta#alpha_{T} < 90^{o}" },				
		{ "DeltaPn_DeltaAlphaT_90_00To135_00Plot",  "90^{o} < #delta#alpha_{T} < 135^{o}" },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot",  "135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot",  "#alpha_{3D} < 45^{o}" },
		{ "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot",  "45^{o} < #alpha_{3D} < 90^{o}" },				
		{ "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot",  "90^{o} < #alpha_{3D} < 135^{o}" },
		{ "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot",  "135^{o} < #alpha_{3D} < 180^{o}" },
		{ "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot",  "#alpha_{3D,#mu} < 45^{o}" },
		{ "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot",  "45^{o} < #alpha_{3D,#mu} < 90^{o}" },				
		{ "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot",  "90^{o} < #alpha_{3D,#mu} < 135^{o}" },
		{ "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot",  "135^{o} < #alpha_{3D,#mu} < 180^{o}" },								
		{ "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot",  "-1 < cos#theta_{#mu} < 0" },
		{ "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot", "0 < cos#theta_{#mu} < 0.5" },
		{ "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{#mu} < 0.75" },
		{ "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{#mu} < 1" },				
		{ "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot", "-1 < cos#theta_{#mu} < 0" },
		{ "DeltaPT_MuonCosTheta_0_00To0_50Plot", "0 < cos#theta_{#mu} < 0.5" },
		{ "DeltaPT_MuonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{#mu} < 0.75" },
		{ "DeltaPT_MuonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{#mu} < 1" },
		{ "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot", "-1 < cos#theta_{p} < 0" },
		{ "DeltaPT_ProtonCosTheta_0_00To0_50Plot", "0 < cos#theta_{p} < 0.5" },
		{ "DeltaPT_ProtonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{p} < 0.75" },
		{ "DeltaPT_ProtonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{p} < 1" },		
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot", "-1 < cos#theta_{#mu} < 0" },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot", "0 < cos#theta_{#mu} < 0.5" },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{#mu} < 0.75" },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{#mu} < 1" },
		{ "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot", "-1 < cos#theta_{p} < 0" },
		{ "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot", "0 < cos#theta_{p} < 0.5" },
		{ "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot", "0.5 < cos#theta_{p} < 0.75" },
		{ "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot", "0.75 < cos#theta_{p} < 1" },
		{ "DeltaPhiT_DeltaPT_0_00To0_20Plot",  "#deltap_{T} < 0.2 GeV/c" },
		{ "DeltaPhiT_DeltaPT_0_20To0_40Plot",  "0.2 < #deltap_{T} < 0.4 GeV/c" },				
		{ "DeltaPhiT_DeltaPT_0_40To1_00Plot",  "#deltap_{T} > 0.4 GeV/c" },	
//		{ "DeltaPn_DeltaPT_0_00To0_20Plot",  "#deltap_{T} < 0.2 GeV/c" },
//		{ "DeltaPn_DeltaPT_0_20To0_40Plot",  "0.2 < #deltap_{T} < 0.4 GeV/c" },				
//		{ "DeltaPn_DeltaPT_0_40To1_00Plot",  "#deltap_{T} > 0.4 GeV/c" },	
		{ "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot",  "#deltap_{T,x} < -0.15 GeV/c" },
		{ "DeltaPty_DeltaPtx_Minus0_15To0_15Plot",  "|#deltap_{T,x}| < 0.15 GeV/c" },				
		{ "DeltaPty_DeltaPtx_0_15To0_55Plot",  "#deltap_{T,x} > 0.15 GeV/c" },
		{ "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot","#deltap_{T,y} < -0.15 GeV/c" },
		{ "DeltaPtx_DeltaPty_Minus0_15To0_15Plot","|#deltap_{T,y}| < 0.15 GeV/c" },
		{ "DeltaPtx_DeltaPty_0_15To0_45Plot","#deltap_{T,y} > 0.15 GeV/c" },		
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot","#deltap_{T} < 0.2 GeV/c, #delta#alpha_{T} < 45^{o}" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot","#deltap_{T} < 0.2 GeV/c, 45^{o} < #delta#alpha_{T} < 90^{o}" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot","#deltap_{T} < 0.2 GeV/c, 90^{o} < #delta#alpha_{T} < 135^{o}" },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot","#deltap_{T} < 0.2 GeV/c, 135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot","0.2 < #deltap_{T} < 0.4 GeV/c, #delta#alpha_{T} < 45^{o}" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot","0.2 < #deltap_{T} < 0.4 GeV/c, 45^{o} < #delta#alpha_{T} < 90^{o}" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot","0.2 < #deltap_{T} < 0.4 GeV/c, 90^{o} < #delta#alpha_{T} < 135^{o}" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot","0.2 < #deltap_{T} < 0.4 GeV/c, 135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot","#deltap_{T} > 0.4 GeV/c, #delta#alpha_{T} < 45^{o}" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot","#deltap_{T} > 0.4 GeV/c, 45^{o} < #delta#alpha_{T} < 90^{o}" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot","#deltap_{T} > 0.4 GeV/c, 90^{o} < #delta#alpha_{T} < 135^{o}" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot","#deltap_{T} > 0.4 GeV/c, 135^{o} < #delta#alpha_{T} < 180^{o}" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot","#deltap_{T,x} < -0.15 GeV/c, #deltap_{T,y} < -0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot","#deltap_{T,x} < -0.15 GeV/c, |#deltap_{T,y}| < 0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot","#deltap_{T,x} < -0.15 GeV/c, #deltap_{T,y} > 0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot","|#deltap_{T,x}| < 0.15 GeV/c, #deltap_{T,y} < -0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot","|#deltap_{T,x}| < 0.15 GeV/c, |#deltap_{T,y}| < 0.15 GeV/c" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot","|#deltap_{T,x}| < 0.15 GeV/c, #deltap_{T,y} > 0.15 GeV/c" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot","#deltap_{T,x} > 0.15 GeV/c, #deltap_{T,y} < -0.15 GeV/c" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot","#deltap_{T,x} > 0.15 GeV/c, |#deltap_{T,y}| < 0.15 GeV/c" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot","#deltap_{T,x} > 0.15 GeV/c, #deltap_{T,y} > 0.15 GeV/c" },			
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot","-1 < cos#theta_{#mu} < 0, 0.1 < P_{#mu} < 0.4 GeV/c" },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot","0 < cos#theta_{#mu} < 0.5, 0.1 < P_{#mu} < 0.4 GeV/c"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot","0.5 < cos#theta_{#mu} < 0.75, 0.1 < P_{#mu} < 0.4 GeV/c" },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot","0.75 < cos#theta_{#mu} < 1, 0.1 < P_{#mu} < 0.4 GeV/c"},
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot","-1 < cos#theta_{#mu} < 0, 0.4 < P_{#mu} < 0.6 GeV/c"},
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot","0 < cos#theta_{#mu} < 0.5, 0.4 < P_{#mu} < 0.6 GeV/c"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot","0.5 < cos#theta_{#mu} < 0.75, 0.4 < P_{#mu} < 0.6 GeV/c"},
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot","0.75 < cos#theta_{#mu} < 1, 0.4 < P_{#mu} < 0.6 GeV/c"},
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot","-1 < cos#theta_{#mu} < 0, 0.6 < P_{#mu} < 1.2 GeV/c"},
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot","0 < cos#theta_{#mu} < 0.5, 0.6 < P_{#mu} < 1.2 GeV/c"},
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot","0.5 < cos#theta_{#mu} < 0.75, 0.6 < P_{#mu} < 1.2 GeV/c"},
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot","0.75 < cos#theta_{#mu} < 1, 0.6 < P_{#mu} < 1.2 GeV/c"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot","-1 < cos#theta_{p} < 0, 0.3 < P_{p} < 0.5 GeV/c"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot","0 < cos#theta_{p} < 0.5, 0.3 < P_{p} < 0.5 GeV/c"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot","0.5 < cos#theta_{p} < 0.75, 0.3 < P_{p} < 0.5 GeV/c"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot","0.75 < cos#theta_{p} < 1, 0.3 < P_{p} < 0.5 GeV/c"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot","-1 < cos#theta_{p} < 0, 0.5 < P_{p} < 0.7 GeV/c"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot","0 < cos#theta_{p} < 0.5, 0.5 < P_{p} < 0.7 GeV/c"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot","0.5 < cos#theta_{p} < 0.75, 0.5 < P_{p} < 0.7 GeV/c"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot","0.75 < cos#theta_{p} < 1, 0.5 < P_{p} < 0.7 GeV/c"},
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot","-1 < cos#theta_{p} < 0, 0.7 < P_{p} < 1 GeV/c"},
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot","0 < cos#theta_{p} < 0.5, 0.7 < P_{p} < 1 GeV/c"},
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot","0.5 < cos#theta_{p} < 0.75, 0.7 < P_{p} < 1 GeV/c"},
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot","0.75 < cos#theta_{p} < 1, 0.7 < P_{p} < 1 GeV/c"},
		{ "ECal_DeltaAlphaT_0_00To45_00Plot","#delta#alpha_{T} < 45^{o}" },
		{ "ECal_DeltaAlphaT_45_00To90_00Plot","45^{o} < #delta#alpha_{T} < 90^{o}" },
		{ "ECal_DeltaAlphaT_90_00To135_00Plot","90^{o} < #delta#alpha_{T} < 135^{o}" },	
		{ "ECal_DeltaAlphaT_135_00To180_00Plot","135^{o} < #delta#alpha_{T} < 180^{o}" },		
		{ "SerialDeltaPT_MuonCosThetaPlot", "cos#theta_{#mu} slices" },
		{ "SerialDeltaPT_ProtonCosThetaPlot", "cos#theta_{p} slices" },
		{ "SerialMuonMomentum_MuonCosThetaPlot", "cos#theta_{#mu} slices" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot", "cos#theta_{p} slices" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot", "cos#theta_{#mu} slices" },					
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot", "cos#theta_{p} slices" },
		{ "SerialDeltaAlphaT_DeltaPTPlot", "#deltap_{T} slices" },
		{ "SerialDeltaAlphaT_MuonMomentumPlot", "p_{#mu} slices" },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot", "p_{p} slices" },				
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot", "p_{n} slices" },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot", "p_{n} slices" },		
		{ "SerialDeltaPT_DeltaAlphaTPlot", "#delta#alpha_{T} slices" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot", "#alpha_{3D} slices" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot", "#alpha_{3D,#mu} slices" },				
		{ "SerialProtonMomentum_DeltaAlphaTPlot", "#delta#alpha_{T} slices" },		
		{ "SerialDeltaPn_DeltaAlphaTPlot", "#delta#alpha_{T} slices" },		
		{ "SerialDeltaPhiT_DeltaPTPlot", "#deltap_{T} slices" },
//		{ "SerialDeltaPn_DeltaPTPlot",  "#deltap_{T} slices" },
		{ "SerialECal_DeltaPTPlot", "#deltap_{T} slices" },
		{ "SerialECal_DeltaPtxPlot", "#deltap_{T,x} slices" },
		{ "SerialECal_DeltaPtyPlot", "#deltap_{T,y} slices" },				
		{ "SerialECal_DeltaAlphaTPlot", "#delta#alpha_{T} slices" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot", "#deltap_{T} & #delta#alpha_{T} slices" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot", "#deltap_{T,x} & #deltap_{T,y} slices" },		
		{ "SerialECal_MuonCosThetaMuonMomentumPlot", "cos#theta_{#mu} & p_{#mu} slices" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot", "cos#theta_{p} & p_{p} slices" },								
		{ "SerialProtonCosTheta_MuonCosThetaPlot", "cos#theta_{#mu} slices" },
		{ "SerialDeltaPty_DeltaPtxPlot", "#deltap_{T,x} slices" },
		{ "SerialDeltaPtx_DeltaPtyPlot", "#deltap_{T,y} slices" },	
		{ "VertexXPlot",  "All events" },
		{ "VertexYPlot",  "All events" },
		{ "VertexZPlot",  "All events" },			
	};	
	
	//----------------------------------------//
	
	static std::map<TString,TString> MapUncorCor =
	{

		{ "VertexXPlot", "VertexXPlot" },
		{ "VertexYPlot", "VertexYPlot" },
		{ "VertexZPlot", "VertexZPlot" },				
		{ "MuonCosThetaPlot", "MuonCosThetaPlot" },
		{ "MuonPhiPlot", "MuonPhiPlot" },		
		{ "MuonMomentumPlot", "MuonMomentumPlot" },
		{ "ProtonCosThetaPlot", "ProtonCosThetaPlot" },
		{ "ProtonPhiPlot", "ProtonPhiPlot" },		
		{ "ProtonMomentumPlot", "ProtonMomentumPlot" },				
		{ "DeltaPTPlot", "DeltaPTPlot" },
		{ "DeltaPnPlot", "DeltaPnPlot" },
		{ "DeltaPnPerpPlot", "DeltaPnPerpPlot" },
		{ "DeltaPnPerpxPlot", "DeltaPnPerpxPlot" },
		{ "DeltaPnPerpyPlot", "DeltaPnPerpyPlot" },				
		{ "DeltaPnParPlot", "DeltaPnParPlot" },						
		{ "DeltaAlphaTPlot", "DeltaAlphaTPlot" },
		{ "DeltaAlpha3DqPlot", "DeltaAlpha3DqPlot" },
		{ "DeltaAlpha3DMuPlot", "DeltaAlpha3DMuPlot" },				
		{ "DeltaPhiTPlot", "DeltaPhiTPlot" },
		{ "DeltaPhi3DPlot", "DeltaPhi3DPlot" },		
		{ "DeltaPtxPlot", "DeltaPtxPlot" },
		{ "DeltaPtyPlot", "DeltaPtyPlot" },
		{ "ECalPlot", "ECalPlot" },
		{ "EQEPlot", "EQEPlot" },		
		{ "Q2Plot", "Q2Plot" },	
		
		{ "SerialDeltaPT_MuonCosThetaPlot", "SerialDeltaPT_MuonCosThetaPlot" },
		{ "SerialDeltaPT_ProtonCosThetaPlot", "SerialDeltaPT_ProtonCosThetaPlot" },
		{ "SerialMuonMomentum_MuonCosThetaPlot", "SerialMuonMomentum_MuonCosThetaPlot" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot", "SerialProtonMomentum_ProtonCosThetaPlot" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot", "SerialDeltaAlphaT_MuonCosThetaPlot" },					
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot", "SerialDeltaAlphaT_ProtonCosThetaPlot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot", "SerialDeltaAlphaT_DeltaPTPlot" },
		{ "SerialDeltaAlphaT_MuonMomentumPlot", "SerialDeltaAlphaT_MuonMomentumPlot" },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot", "SerialDeltaAlphaT_ProtonMomentumPlot" },				
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot", "SerialDeltaAlpha3Dq_DeltaPnPlot" },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot", "SerialDeltaAlpha3DMu_DeltaPnPlot" },		
		{ "SerialDeltaPT_DeltaAlphaTPlot", "SerialDeltaPT_DeltaAlphaTPlot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot", "SerialDeltaPn_DeltaAlpha3DqPlot" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot", "SerialDeltaPn_DeltaAlpha3DMuPlot" },				
		{ "SerialProtonMomentum_DeltaAlphaTPlot", "SerialProtonMomentum_DeltaAlphaTPlot" },		
		{ "SerialDeltaPn_DeltaAlphaTPlot", "SerialDeltaPn_DeltaAlphaTPlot" },		
		{ "SerialDeltaPhiT_DeltaPTPlot", "SerialDeltaPhiT_DeltaPTPlot" },
//		{ "SerialDeltaPn_DeltaPTPlot",  "SerialDeltaPn_DeltaPTPlot" },
		{ "SerialECal_DeltaPTPlot", "SerialECal_DeltaPTPlot" },
		{ "SerialECal_DeltaPtxPlot", "SerialECal_DeltaPtxPlot" },
		{ "SerialECal_DeltaPtyPlot", "SerialECal_DeltaPtyPlot" },				
		{ "SerialECal_DeltaAlphaTPlot", "SerialECal_DeltaAlphaTPlot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot", "SerialECal_DeltaPTDeltaAlphaTPlot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot", "SerialECal_DeltaPtxDeltaPtyPlot" },		
		{ "SerialECal_MuonCosThetaMuonMomentumPlot", "SerialECal_MuonCosThetaMuonMomentumPlot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot", "SerialECal_ProtonCosThetaProtonMomentumPlot" },								
		{ "SerialProtonCosTheta_MuonCosThetaPlot", "SerialProtonCosTheta_MuonCosThetaPlot" },									

		{ "SerialDeltaPT_MuonCosThetaPlot_0", "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialDeltaPT_MuonCosThetaPlot_1", "DeltaPT_MuonCosTheta_0_00To0_50Plot" },
		{ "SerialDeltaPT_MuonCosThetaPlot_2", "DeltaPT_MuonCosTheta_0_50To0_75Plot" },
		{ "SerialDeltaPT_MuonCosThetaPlot_3", "DeltaPT_MuonCosTheta_0_75To1_00Plot" },
								
		{ "SerialDeltaPT_ProtonCosThetaPlot_0", "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialDeltaPT_ProtonCosThetaPlot_1", "DeltaPT_ProtonCosTheta_0_00To0_50Plot" },
		{ "SerialDeltaPT_ProtonCosThetaPlot_2", "DeltaPT_ProtonCosTheta_0_50To0_75Plot" },	
		{ "SerialDeltaPT_ProtonCosThetaPlot_3", "DeltaPT_ProtonCosTheta_0_75To1_00Plot" },
						
		{ "SerialMuonMomentum_MuonCosThetaPlot_0", "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialMuonMomentum_MuonCosThetaPlot_1", "MuonMomentum_MuonCosTheta_0_00To0_50Plot" },
		{ "SerialMuonMomentum_MuonCosThetaPlot_2", "MuonMomentum_MuonCosTheta_0_50To0_75Plot" },
		{ "SerialMuonMomentum_MuonCosThetaPlot_3", "MuonMomentum_MuonCosTheta_0_75To1_00Plot" },
								
		{ "SerialProtonMomentum_ProtonCosThetaPlot_0", "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot_1", "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot_2", "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot" },
		{ "SerialProtonMomentum_ProtonCosThetaPlot_3", "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot" },						
		
		{ "SerialDeltaAlphaT_MuonCosThetaPlot_0", "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot_1", "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot_2", "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot" },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot_3", "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot" },
		
		{ "SerialDeltaAlphaT_MuonMomentumPlot_0", "DeltaAlphaT_MuonMomentum_0_10To0_40Plot" },
		{ "SerialDeltaAlphaT_MuonMomentumPlot_1", "DeltaAlphaT_MuonMomentum_0_40To0_60Plot" },
		{ "SerialDeltaAlphaT_MuonMomentumPlot_2", "DeltaAlphaT_MuonMomentum_0_60To1_20Plot" },						
							
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_0", "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_1", "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_2", "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot" },
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot_3", "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot" },	
		
		{ "SerialDeltaAlphaT_ProtonMomentumPlot_0", "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot" },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot_1", "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot" },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot_2", "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot" },					
		
		{ "SerialDeltaAlphaT_DeltaPTPlot_0", "DeltaAlphaT_DeltaPT_0_00To0_20Plot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot_1", "DeltaAlphaT_DeltaPT_0_20To0_40Plot" },
		{ "SerialDeltaAlphaT_DeltaPTPlot_2", "DeltaAlphaT_DeltaPT_0_40To1_00Plot" },			
		
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_0", "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot" },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_1", "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot" },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot_2", "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot" },
		
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot_0", "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot" },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot_1", "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot" },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot_2", "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot" },				

		{ "SerialProtonMomentum_DeltaAlphaTPlot_0", "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialProtonMomentum_DeltaAlphaTPlot_1", "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialProtonMomentum_DeltaAlphaTPlot_2", "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot" },
		{ "SerialProtonMomentum_DeltaAlphaTPlot_3", "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot" },
						
		{ "SerialDeltaPT_DeltaAlphaTPlot_0", "DeltaPT_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialDeltaPT_DeltaAlphaTPlot_1", "DeltaPT_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialDeltaPT_DeltaAlphaTPlot_2", "DeltaPT_DeltaAlphaT_90_00To135_00Plot" },
		{ "SerialDeltaPT_DeltaAlphaTPlot_3", "DeltaPT_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "SerialDeltaPn_DeltaAlphaTPlot_0", "DeltaPn_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialDeltaPn_DeltaAlphaTPlot_1", "DeltaPn_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialDeltaPn_DeltaAlphaTPlot_2", "DeltaPn_DeltaAlphaT_90_00To135_00Plot" },
		{ "SerialDeltaPn_DeltaAlphaTPlot_3", "DeltaPn_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_0", "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_1", "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_2", "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot_3", "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" },
		
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot_0", "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot_1", "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot_2", "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot" },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot_3", "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot" },		
										
		{ "SerialDeltaPhiT_DeltaPTPlot_0", "DeltaPhiT_DeltaPT_0_00To0_20Plot" },
		{ "SerialDeltaPhiT_DeltaPTPlot_1", "DeltaPhiT_DeltaPT_0_20To0_40Plot" },
		{ "SerialDeltaPhiT_DeltaPTPlot_2", "DeltaPhiT_DeltaPT_0_40To1_00Plot" },
								
//		{ "SerialDeltaPn_DeltaPTPlot_0",  "DeltaPn_DeltaPT_0_00To0_20Plot" },
//		{ "SerialDeltaPn_DeltaPTPlot_1",  "DeltaPn_DeltaPT_0_20To0_40Plot" },
//		{ "SerialDeltaPn_DeltaPTPlot_2",  "DeltaPn_DeltaPT_0_40To1_00Plot" },
								
		{ "SerialECal_DeltaPTPlot_0", "ECal_DeltaPT_0_00To0_20Plot" },
		{ "SerialECal_DeltaPTPlot_1", "ECal_DeltaPT_0_20To0_40Plot" },
		{ "SerialECal_DeltaPTPlot_2", "ECal_DeltaPT_0_40To1_00Plot" },
		
		{ "SerialECal_DeltaPtxPlot_0", "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot" },
		{ "SerialECal_DeltaPtxPlot_1", "ECal_DeltaPtx_Minus0_15To0_15Plot" },
		{ "SerialECal_DeltaPtxPlot_2", "ECal_DeltaPtx_0_15To0_55Plot" },
		
		{ "SerialECal_DeltaPtyPlot_0", "ECal_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialECal_DeltaPtyPlot_1", "ECal_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialECal_DeltaPtyPlot_2", "ECal_DeltaPty_0_15To0_45Plot" },										
		
		{ "SerialECal_DeltaAlphaTPlot_0", "ECal_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialECal_DeltaAlphaTPlot_1", "ECal_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialECal_DeltaAlphaTPlot_2", "ECal_DeltaAlphaT_90_00To135_00Plot" },	
		{ "SerialECal_DeltaAlphaTPlot_3", "ECal_DeltaAlphaT_135_00To180_00Plot" },

		{ "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot", "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "DeltaPT_MuonCosTheta_0_00To0_50Plot", "DeltaPT_MuonCosTheta_0_00To0_50Plot" },
		{ "DeltaPT_MuonCosTheta_0_50To0_75Plot", "DeltaPT_MuonCosTheta_0_50To0_75Plot" },
		{ "DeltaPT_MuonCosTheta_0_75To1_00Plot", "DeltaPT_MuonCosTheta_0_75To1_00Plot" },
								
		{ "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot", "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "DeltaPT_ProtonCosTheta_0_00To0_50Plot", "DeltaPT_ProtonCosTheta_0_00To0_50Plot" },
		{ "DeltaPT_ProtonCosTheta_0_50To0_75Plot", "DeltaPT_ProtonCosTheta_0_50To0_75Plot" },	
		{ "DeltaPT_ProtonCosTheta_0_75To1_00Plot", "DeltaPT_ProtonCosTheta_0_75To1_00Plot" },
						
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot", "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot", "MuonMomentum_MuonCosTheta_0_00To0_50Plot" },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot", "MuonMomentum_MuonCosTheta_0_50To0_75Plot" },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot", "MuonMomentum_MuonCosTheta_0_75To1_00Plot" },
								
		{ "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot", "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot", "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot" },
		{ "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot", "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot" },
		{ "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot", "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot" },						
		
		{ "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot", "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot", "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot" },
		{ "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot", "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot" },
		{ "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot", "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot" },
		
		{ "DeltaAlphaT_MuonMomentum_0_10To0_40Plot", "DeltaAlphaT_MuonMomentum_0_10To0_40Plot" },
		{ "DeltaAlphaT_MuonMomentum_0_40To0_60Plot", "DeltaAlphaT_MuonMomentum_0_40To0_60Plot" },
		{ "DeltaAlphaT_MuonMomentum_0_60To1_20Plot", "DeltaAlphaT_MuonMomentum_0_60To1_20Plot" },						
							
		{ "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot", "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot" },
		{ "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot", "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot" },
		{ "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot", "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot" },
		{ "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot", "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot" },	
		
		{ "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot", "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot" },
		{ "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot", "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot" },
		{ "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot", "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot" },					
		
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot", "DeltaAlphaT_DeltaPT_0_00To0_20Plot" },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot", "DeltaAlphaT_DeltaPT_0_20To0_40Plot" },
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot", "DeltaAlphaT_DeltaPT_0_40To1_00Plot" },			
		
		{ "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot", "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot" },
		{ "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot", "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot" },
		{ "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot", "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot" },
		
		{ "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot", "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot" },
		{ "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot", "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot" },
		{ "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot", "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot" },		

		{ "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot", "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot" },
		{ "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot", "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot" },
		{ "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot", "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot" },
		{ "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot", "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot" },
						
		{ "DeltaPT_DeltaAlphaT_0_00To45_00Plot", "DeltaPT_DeltaAlphaT_0_00To45_00Plot" },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot", "DeltaPT_DeltaAlphaT_45_00To90_00Plot" },
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot", "DeltaPT_DeltaAlphaT_90_00To135_00Plot" },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot", "DeltaPT_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "DeltaPn_DeltaAlphaT_0_00To45_00Plot", "DeltaPn_DeltaAlphaT_0_00To45_00Plot" },
		{ "DeltaPn_DeltaAlphaT_45_00To90_00Plot", "DeltaPn_DeltaAlphaT_45_00To90_00Plot" },
		{ "DeltaPn_DeltaAlphaT_90_00To135_00Plot", "DeltaPn_DeltaAlphaT_90_00To135_00Plot" },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot", "DeltaPn_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot", "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot" },
		{ "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot", "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot" },
		{ "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot", "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot" },
		{ "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot", "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot" },
		
		{ "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot", "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot" },
		{ "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot", "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot" },
		{ "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot", "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot" },
		{ "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot", "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot" },						
										
		{ "DeltaPhiT_DeltaPT_0_00To0_20Plot", "DeltaPhiT_DeltaPT_0_00To0_20Plot" },
		{ "DeltaPhiT_DeltaPT_0_20To0_40Plot", "DeltaPhiT_DeltaPT_0_20To0_40Plot" },
		{ "DeltaPhiT_DeltaPT_0_40To1_00Plot", "DeltaPhiT_DeltaPT_0_40To1_00Plot" },
								
//		{ "DeltaPn_DeltaPT_0_00To0_20Plot",  "DeltaPn_DeltaPT_0_00To0_20Plot" },
//		{ "DeltaPn_DeltaPT_0_20To0_40Plot",  "DeltaPn_DeltaPT_0_20To0_40Plot" },
//		{ "DeltaPn_DeltaPT_0_40To1_00Plot",  "DeltaPn_DeltaPT_0_40To1_00Plot" },
								
		{ "ECal_DeltaPT_0_00To0_20Plot", "ECal_DeltaPT_0_00To0_20Plot" },
		{ "ECal_DeltaPT_0_20To0_40Plot", "ECal_DeltaPT_0_20To0_40Plot" },
		{ "ECal_DeltaPT_0_40To1_00Plot", "ECal_DeltaPT_0_40To1_00Plot" },
		
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot", "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot" },
		{ "ECal_DeltaPtx_Minus0_15To0_15Plot", "ECal_DeltaPtx_Minus0_15To0_15Plot" },
		{ "ECal_DeltaPtx_0_15To0_55Plot", "ECal_DeltaPtx_0_15To0_55Plot" },
		
		{ "ECal_DeltaPty_Minus0_75ToMinus0_15Plot", "ECal_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "ECal_DeltaPty_Minus0_15To0_15Plot", "ECal_DeltaPty_Minus0_15To0_15Plot" },
		{ "ECal_DeltaPty_0_15To0_45Plot", "ECal_DeltaPty_0_15To0_45Plot" },										
		
		{ "ECal_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaAlphaT_0_00To45_00Plot" },
		{ "ECal_DeltaAlphaT_45_00To90_00Plot", "ECal_DeltaAlphaT_45_00To90_00Plot" },
		{ "ECal_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaAlphaT_90_00To135_00Plot" },	
		{ "ECal_DeltaAlphaT_135_00To180_00Plot", "ECal_DeltaAlphaT_135_00To180_00Plot" },		
		
		
		
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_0", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_1", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_2", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot" },	
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_3", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_4", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_5", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_6", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot" },	
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_7", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_8", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_9", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot" },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_10", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot" },	
		{ "SerialECal_DeltaPTDeltaAlphaTPlot_11", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "SerialECal_DeltaPtxDeltaPtyPlot_0", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_1", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_2", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot" },	
		{ "SerialECal_DeltaPtxDeltaPtyPlot_3", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_4", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_5", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_6", "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot" },	
		{ "SerialECal_DeltaPtxDeltaPtyPlot_7", "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialECal_DeltaPtxDeltaPtyPlot_8", "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot" },
					
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_0", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_1", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_2", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot" },	
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_3", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_4", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_5", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_6", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_7", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_8", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_9", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_10", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot" },
		{ "SerialECal_MuonCosThetaMuonMomentumPlot_11", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot" },
		
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_0", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_1", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_2", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot" },	
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_3", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_4", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_5", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_6", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_7", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_8", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_9", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_10", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot" },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot_11", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot" },							
						
		{ "SerialProtonCosTheta_MuonCosThetaPlot_0", "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "SerialProtonCosTheta_MuonCosThetaPlot_1", "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot" },
		{ "SerialProtonCosTheta_MuonCosThetaPlot_2", "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot" },
		{ "SerialProtonCosTheta_MuonCosThetaPlot_3", "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot" },
								
		{ "SerialDeltaPty_DeltaPtxPlot_0", "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot" },
		{ "SerialDeltaPty_DeltaPtxPlot_1", "DeltaPty_DeltaPtx_Minus0_15To0_15Plot" },
		{ "SerialDeltaPty_DeltaPtxPlot_2", "DeltaPty_DeltaPtx_0_15To0_55Plot" },
						
		{ "SerialDeltaPtx_DeltaPtyPlot_0", "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_1", "DeltaPtx_DeltaPty_Minus0_15To0_15Plot" },
		{ "SerialDeltaPtx_DeltaPtyPlot_2", "DeltaPtx_DeltaPty_0_15To0_45Plot" },
		
		
		
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot" },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot" },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot" },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot" },	
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot" },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot" },	
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot" },
		
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot" },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot", "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot" },	
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot" },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot", "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot", "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot" },	
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot", "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot" },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot", "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot" },
					
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot" },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot" },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot", "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot" },	
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot" },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot" },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot" },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot" },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot" },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot", "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot" },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot" },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot" },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot", "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot" },
		
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot" },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot" },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot", "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot" },	
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot" },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot" },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot", "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot" },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot" },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot" },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot", "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot" },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot" },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot" },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot", "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot" },					
						
		{ "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot", "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot" },
		{ "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot", "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot" },
		{ "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot", "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot" },
		{ "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot", "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot" },
								
		{ "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot", "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot" },
		{ "DeltaPty_DeltaPtx_Minus0_15To0_15Plot", "DeltaPty_DeltaPtx_Minus0_15To0_15Plot" },
		{ "DeltaPty_DeltaPtx_0_15To0_55Plot", "DeltaPty_DeltaPtx_0_15To0_55Plot" },
						
		{ "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot", "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot" },
		{ "DeltaPtx_DeltaPty_Minus0_15To0_15Plot", "DeltaPtx_DeltaPty_Minus0_15To0_15Plot" },
		{ "DeltaPtx_DeltaPty_0_15To0_45Plot", "DeltaPtx_DeltaPty_0_15To0_45Plot" },		
						
	};	

	// -------------------------------------------------------------------------------------------------------------------------

	// Genie Constants

	static const double SuSav2FluxIntegratedXSection = 37.2171; // e-38 cm^2
//	static const double FluxIntegratedXSection = 26.5732; // e-38 cm^2 MCC8 number
	static const double G18_10a_02_11a_FluxIntegratedXSection = 27.1682; // e-38 cm^2 MCC9 number GENIE v3.0.6
	static const double R_2_12_10_FluxIntegratedXSection = 35.116; // e-38 cm^2 MCC9 number GENIE v2.12.10
	static const double R_2_12_10_LFG_FluxIntegratedXSection = 36.5471; // e-38 cm^2 MCC9 number GENIE v2.12.10 w/ LFG
	static const double R_2_12_10_EffSF_FluxIntegratedXSection = 35.0266; // e-38 cm^2 MCC9 number GENIE v2.12.10 w/ EffSF		
	static const double R_3_0_4_FluxIntegratedXSection = 26.7883; // e-38 cm^2 MCC9 number GENIE v3.0.4
	static const double R_3_0_6_G18_10a_02_11a_NoCoulomb_FluxIntegratedXSection = 26.134; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a w/o Coulomb corrections
	static const double R_3_0_6_G18_10a_02_11a_NoRPA_FluxIntegratedXSection = 27.7777; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a w/o RPA corrections
	static const double R_3_0_6_G18_10a_02_11a_RFG_FluxIntegratedXSection = 25.76; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a and RFG instead of LFG
	static const double R_3_0_6_G18_10a_02_11a_RFG_NoRPA_FluxIntegratedXSection = 27.2908; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a and RFG instead of LFG w/o RPA	
	static const double R_3_0_6_G18_10a_02_11a_EffSF_FluxIntegratedXSection = 24.5238; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a and EffSF instead of LFG
	static const double R_3_0_6_G18_10a_02_11a_EffSF_NoRPA_FluxIntegratedXSection = 25.3184; // e-38 cm^2 MCC9 number GENIE v3.0.6 / G18_10a_02_11a and EffSF instead of LFG	w/o RPA
	static const double v3_2_0_G21_11b_00_000_FluxIntegratedXSection = 31.9663; // e-38 cm^2 MCC9
	static const double v3_4_0_G18_02a_02_11a_FluxIntegratedXSection = 27.2797; // e-38 cm^2 MCC9
	static const double v3_4_0_G21_11b_00_000_FluxIntegratedXSection = 31.9544; // e-38 cm^2 MCC9
	static const double v3_0_6_G18_10d_02_11a_FluxIntegratedXSection = 27.1433; // e-38 cm^2 MCC9
	static const double v3_2_0_G18_10d_02_11a_FluxIntegratedXSection = 27.144; // e-38 cm^2 MCC9				

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

	// -------------------------------------------------------------------------------------------------------------------------

	// Plots to be included for xsec extraction purposes

	vector<TString> PlotNames{
				  "DeltaPTPlot"
				 ,"DeltaAlphaTPlot"
				 ,"DeltaAlpha3DqPlot"
				 ,"DeltaAlpha3DMuPlot"
				 ,"DeltaPhiTPlot"
				 ,"DeltaPhi3DPlot"
				 ,"MuonMomentumPlot"
				 ,"MuonCosThetaPlot"
				 ,"MuonCosThetaSingleBinPlot"
				 ,"MuonPhiPlot"
				 ,"ProtonMomentumPlot"
				 ,"ProtonCosThetaPlot","ProtonPhiPlot","Q2Plot","ECalPlot","EQEPlot"
				 ,"ECal_DeltaPT_0_00To0_20Plot","ECal_DeltaPT_0_20To0_40Plot","ECal_DeltaPT_0_40To1_00Plot"
				 ,"ECal_DeltaPtx_Minus0_55ToMinus0_15Plot","ECal_DeltaPtx_Minus0_15To0_15Plot","ECal_DeltaPtx_0_15To0_55Plot"		 
				 ,"ECal_DeltaPty_Minus0_75ToMinus0_15Plot","ECal_DeltaPty_Minus0_15To0_15Plot","ECal_DeltaPty_0_15To0_45Plot"								 
				 ,"CCQEMuonMomentumPlot"
				 ,"CCQEMuonCosThetaPlot"
				 ,"CCQEProtonMomentumPlot"
				 ,"CCQEProtonCosThetaPlot"
				 ,"CCQEQ2Plot"
				 ,"CCQEECalPlot"
				 //,"DeltaPLPlot"
				 ,"DeltaPnPlot"
				 ,"DeltaPnPerpPlot"
				 ,"DeltaPnPerpxPlot"
				 ,"DeltaPnPerpyPlot"				 				 
				 ,"DeltaPnParPlot"
				 ,"DeltaPtxPlot"
				 ,"DeltaPtyPlot"
				 //,"APlot","kMissPlot","PMissPlot","PMissMinusPlot" 	 
				 ,"DeltaAlphaT_DeltaPT_0_00To0_20Plot","DeltaAlphaT_DeltaPT_0_20To0_40Plot","DeltaAlphaT_DeltaPT_0_40To1_00Plot"
				 ,"DeltaAlphaT_MuonMomentum_0_10To0_40Plot","DeltaAlphaT_MuonMomentum_0_40To0_60Plot","DeltaAlphaT_MuonMomentum_0_60To1_20Plot"
				 ,"DeltaAlphaT_ProtonMomentum_0_30To0_50Plot","DeltaAlphaT_ProtonMomentum_0_50To0_70Plot","DeltaAlphaT_ProtonMomentum_0_70To1_00Plot"	 		
				 ,"DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot","DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot","DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
				 ,"DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot","DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot","DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot"				 
				 ,"DeltaPT_DeltaAlphaT_0_00To45_00Plot","DeltaPT_DeltaAlphaT_45_00To90_00Plot","DeltaPT_DeltaAlphaT_90_00To135_00Plot","DeltaPT_DeltaAlphaT_135_00To180_00Plot"
				 ,"DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot","DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot"
				 ,"DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot","DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot"
				 ,"DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot","DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot"
				 ,"DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot","DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot"				 				 
				 ,"ProtonMomentum_DeltaAlphaT_0_00To45_00Plot","ProtonMomentum_DeltaAlphaT_45_00To90_00Plot","ProtonMomentum_DeltaAlphaT_90_00To135_00Plot","ProtonMomentum_DeltaAlphaT_135_00To180_00Plot"				 
				 //,"DeltaPn_DeltaAlphaT_0_00To45_00Plot","DeltaPn_DeltaAlphaT_45_00To90_00Plot","DeltaPn_DeltaAlphaT_90_00To135_00Plot","DeltaPn_DeltaAlphaT_135_00To180_00Plot"				 
				 ,"DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_MuonCosTheta_0_00To0_50Plot","DeltaAlphaT_MuonCosTheta_0_50To0_75Plot","DeltaAlphaT_MuonCosTheta_0_75To1_00Plot"
				 ,"DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot","DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot","DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot"				 
				 ,"DeltaPT_MuonCosTheta_Minus1_00To0_00Plot","DeltaPT_MuonCosTheta_0_00To0_50Plot","DeltaPT_MuonCosTheta_0_50To0_75Plot","DeltaPT_MuonCosTheta_0_75To1_00Plot"
				 ,"DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaPT_ProtonCosTheta_0_00To0_50Plot","DeltaPT_ProtonCosTheta_0_50To0_75Plot","DeltaPT_ProtonCosTheta_0_75To1_00Plot"				 
				 ,"ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot","ProtonCosTheta_MuonCosTheta_0_00To0_50Plot","ProtonCosTheta_MuonCosTheta_0_50To0_75Plot","ProtonCosTheta_MuonCosTheta_0_75To1_00Plot"			
				 ,"MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot","MuonMomentum_MuonCosTheta_0_00To0_50Plot","MuonMomentum_MuonCosTheta_0_50To0_75Plot","MuonMomentum_MuonCosTheta_0_75To1_00Plot"				 
				 ,"ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot","ProtonMomentum_ProtonCosTheta_0_00To0_50Plot","ProtonMomentum_ProtonCosTheta_0_50To0_75Plot","ProtonMomentum_ProtonCosTheta_0_75To1_00Plot"				 
				 ,"DeltaPhiT_DeltaPT_0_00To0_20Plot","DeltaPhiT_DeltaPT_0_20To0_40Plot","DeltaPhiT_DeltaPT_0_40To1_00Plot"
//				 ,"DeltaPn_DeltaPT_0_00To0_20Plot","DeltaPn_DeltaPT_0_20To0_40Plot","DeltaPn_DeltaPT_0_40To1_00Plot"				 
				 ,"DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot","DeltaPty_DeltaPtx_Minus0_15To0_15Plot","DeltaPty_DeltaPtx_0_15To0_55Plot"
				 ,"DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot","DeltaPtx_DeltaPty_Minus0_15To0_15Plot","DeltaPtx_DeltaPty_0_15To0_45Plot"
				 ,"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot"
				 ,"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot"
				 ,"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot"
				 ,"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot"
				 ,"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot"
				 ,"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot"	
				 ,"ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot","ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot"
				 ,"ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot","ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot"	
				 ,"ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot","ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot"	
				 ,"ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot","ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot"
				 ,"ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot"	 				 			 
				 ,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot"
				 ,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot"
				 ,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot","ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot"
				 ,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot"
				 ,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot","ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot"
				 ,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot"
				 ,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot"
				 ,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot"
				 ,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot"
				 ,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot"
				 ,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot"
				 ,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot" 
				 ,"ECal_DeltaAlphaT_0_00To45_00Plot"
				 , "ECal_DeltaAlphaT_45_00To90_00Plot"	
				 ,"ECal_DeltaAlphaT_90_00To135_00Plot"
				 , "ECal_DeltaAlphaT_135_00To180_00Plot"				 			 
				 ,"SerialDeltaPT_MuonCosThetaPlot"
				 ,"SerialDeltaPT_ProtonCosThetaPlot"
				 ,"SerialMuonMomentum_MuonCosThetaPlot"
				 ,"SerialProtonMomentum_ProtonCosThetaPlot"
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot"
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot"
				 ,"SerialDeltaAlphaT_DeltaPTPlot"
				 ,"SerialDeltaAlphaT_MuonMomentumPlot"
				 ,"SerialDeltaAlphaT_ProtonMomentumPlot"				 				 
				 ,"SerialDeltaAlpha3Dq_DeltaPnPlot"
				 ,"SerialDeltaAlpha3DMu_DeltaPnPlot"				 				 
				 ,"SerialDeltaPT_DeltaAlphaTPlot"
				 ,"SerialDeltaPn_DeltaAlpha3DqPlot"
				 ,"SerialDeltaPn_DeltaAlpha3DMuPlot"				 				 
				 ,"SerialProtonMomentum_DeltaAlphaTPlot"				 
				 //,"SerialDeltaPn_DeltaAlphaTPlot"				 
				 ,"SerialDeltaPhiT_DeltaPTPlot"
//				 ,"SerialDeltaPn_DeltaPTPlot"
				 ,"SerialECal_DeltaPTPlot"
				 ,"SerialECal_DeltaPtxPlot"
				 ,"SerialECal_DeltaPtyPlot"				 				 				 
				 ,"SerialECal_DeltaAlphaTPlot"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot"			 
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot"				 
				 ,"SerialProtonCosTheta_MuonCosThetaPlot"
				 ,"SerialDeltaPty_DeltaPtxPlot"
				 ,"SerialDeltaPtx_DeltaPtyPlot"
				 ,"VertexXPlot"
				 ,"VertexYPlot"
				 ,"VertexZPlot"
				 				 				 				 
				 };

	// "ECalPlot" "EQEPlot" "Q2Plot" "CCQEECalPlot" "CCQEEQEPlot" "CCQEQ2Plot"
	
	//----------------------------------------//
	
	vector<TString> OneDimXSec = {
				  // 1D
				  "DeltaPTPlot"
				 ,"DeltaAlphaTPlot"
				 ,"DeltaAlpha3DqPlot"
				 ,"DeltaAlpha3DMuPlot"				 				 
				 ,"DeltaPhiTPlot"
				 ,"DeltaPhi3DPlot"
				 ,"MuonMomentumPlot"
				 ,"MuonCosThetaPlot"
				 ,"MuonCosThetaSingleBinPlot"
				 ,"ProtonMomentumPlot"
				 ,"ProtonCosThetaPlot"
				 ,"ECalPlot"
				 ,"EQEPlot"				 
				 ,"DeltaPnPlot"
				 ,"DeltaPnPerpPlot"
				 ,"DeltaPnPerpxPlot"
				 ,"DeltaPnPerpyPlot"				 				 
				 ,"DeltaPnParPlot"				 				 
				 ,"DeltaPtxPlot"
				 ,"DeltaPtyPlot"
				 ,"VertexXPlot"
				 ,"VertexYPlot"
				 ,"VertexZPlot"
				 // 2D Uncorrelated results
				 ,"DeltaPT_DeltaAlphaT_0_00To45_00Plot","DeltaPT_DeltaAlphaT_45_00To90_00Plot"
				 ,"DeltaPT_DeltaAlphaT_90_00To135_00Plot","DeltaPT_DeltaAlphaT_135_00To180_00Plot"
				 ,"DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot","DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot"
				 ,"DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot","DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot"
				 ,"DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot","DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot"
				 ,"DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot","DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot"				 				 	
				 ,"DeltaPT_MuonCosTheta_Minus1_00To0_00Plot","DeltaPT_MuonCosTheta_0_00To0_50Plot"
				 ,"DeltaPT_MuonCosTheta_0_50To0_75Plot","DeltaPT_MuonCosTheta_0_75To1_00Plot"
				 ,"DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaPT_ProtonCosTheta_0_00To0_50Plot"
				 ,"DeltaPT_ProtonCosTheta_0_50To0_75Plot","DeltaPT_ProtonCosTheta_0_75To1_00Plot"
				 ,"DeltaAlphaT_DeltaPT_0_00To0_20Plot","DeltaAlphaT_DeltaPT_0_20To0_40Plot","DeltaAlphaT_DeltaPT_0_40To1_00Plot"
				 ,"DeltaAlphaT_MuonMomentum_0_10To0_40Plot","DeltaAlphaT_MuonMomentum_0_40To0_60Plot","DeltaAlphaT_MuonMomentum_0_60To1_20Plot"
				 ,"DeltaAlphaT_ProtonMomentum_0_30To0_50Plot","DeltaAlphaT_ProtonMomentum_0_50To0_70Plot","DeltaAlphaT_ProtonMomentum_0_70To1_00Plot"
				 ,"DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot","DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot","DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot"
				 ,"DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot","DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot","DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot"				 
				 ,"DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_MuonCosTheta_0_00To0_50Plot"
				 ,"DeltaAlphaT_MuonCosTheta_0_50To0_75Plot","DeltaAlphaT_MuonCosTheta_0_75To1_00Plot"
				 ,"DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot"
				 ,"DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot","DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot"	
				 ,"DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot","DeltaPty_DeltaPtx_Minus0_15To0_15Plot","DeltaPty_DeltaPtx_0_15To0_55Plot"
				 ,"DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot","DeltaPtx_DeltaPty_Minus0_15To0_15Plot","DeltaPtx_DeltaPty_0_15To0_45Plot"
				 ,"DeltaPhiT_DeltaPT_0_00To0_20Plot","DeltaPhiT_DeltaPT_0_20To0_40Plot","DeltaPhiT_DeltaPT_0_40To1_00Plot"		
				 ,"ECal_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaAlphaT_45_00To90_00Plot"	
				 ,"ECal_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaAlphaT_135_00To180_00Plot"	
				 ,"ECal_DeltaPT_0_00To0_20Plot","ECal_DeltaPT_0_20To0_40Plot","ECal_DeltaPT_0_40To1_00Plot"
				 ,"ECal_DeltaPtx_Minus0_55ToMinus0_15Plot","ECal_DeltaPtx_Minus0_15To0_15Plot","ECal_DeltaPtx_0_15To0_55Plot"
				 ,"ECal_DeltaPty_Minus0_75ToMinus0_15Plot","ECal_DeltaPty_Minus0_15To0_15Plot","ECal_DeltaPty_0_15To0_45Plot"
				 // 2D Correlated results
				 ,"SerialDeltaPT_MuonCosThetaPlot"
				 ,"SerialDeltaPT_ProtonCosThetaPlot"
				 ,"SerialProtonMomentum_DeltaAlphaTPlot"
				 ,"SerialDeltaPT_DeltaAlphaTPlot"
				 ,"SerialDeltaPn_DeltaAlpha3DqPlot"
				 ,"SerialDeltaPn_DeltaAlpha3DMuPlot"				 				 
				 ,"SerialMuonMomentum_MuonCosThetaPlot"
				 ,"SerialProtonMomentum_ProtonCosThetaPlot"
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot"
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot"
				 ,"SerialDeltaAlphaT_DeltaPTPlot"
				 ,"SerialDeltaAlphaT_MuonMomentumPlot"
				 ,"SerialDeltaAlphaT_ProtonMomentumPlot"				 
				 ,"SerialDeltaAlpha3Dq_DeltaPnPlot"
				 ,"SerialDeltaAlpha3DMu_DeltaPnPlot"				 				 
				 ,"SerialDeltaPhiT_DeltaPTPlot"
//				 ,"SerialDeltaPn_DeltaPTPlot"
				 //,"SerialDeltaPn_DeltaAlphaTPlot"
				 ,"SerialECal_DeltaPTPlot"
				 ,"SerialECal_DeltaPtxPlot"
				 ,"SerialECal_DeltaPtyPlot"				 				 						 
				 ,"SerialECal_DeltaAlphaTPlot"				 
				 ,"SerialProtonCosTheta_MuonCosThetaPlot"
				 ,"SerialDeltaPty_DeltaPtxPlot"
				 ,"SerialDeltaPtx_DeltaPtyPlot"
				 //3D Correlated results
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot"		 				 
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot"	 				 
				 				 
			};	
	
	//----------------------------------------//
	
	vector<TString> MultiDimXSec = {
				  "DeltaPTPlot"
				 ,"DeltaAlphaTPlot"
				 ,"DeltaAlpha3DqPlot"
				 ,"DeltaAlpha3DMuPlot"				 				 
				 ,"DeltaPhiTPlot"
				 ,"DeltaPhi3DPlot"
				 ,"MuonMomentumPlot"
				 ,"MuonCosThetaPlot"
				 ,"MuonCosThetaSingleBinPlot"
				 //,"MuonPhiPlot"
				 ,"ProtonMomentumPlot"
				 ,"ProtonCosThetaPlot"
				 //,"ProtonPhiPlot"
				 //,"Q2Plot"
				 ,"ECalPlot"
				 ,"EQEPlot"				 
				 //,"ECal_DeltaPT_0_00To0_20Plot","ECal_DeltaPT_0_20To0_40Plot","ECal_DeltaPT_0_40To1_00Plot"
				 //,"ECal_DeltaPtx_Minus0_55ToMinus0_15Plot","ECal_DeltaPtx_Minus0_15To0_15Plot","ECal_DeltaPtx_0_15To0_55Plot"
				 //,"ECal_DeltaPty_Minus0_75ToMinus0_15Plot","ECal_DeltaPty_Minus0_15To0_15Plot","ECal_DeltaPty_0_15To0_45Plot"										 
				 //,"CCQEMuonMomentumPlot","CCQEMuonCosThetaPlot","CCQEProtonMomentumPlot","CCQEProtonCosThetaPlot"
				 //,"CCQEQ2Plot","CCQEECalPlot"
				 //,"DeltaPLPlot"
				 ,"DeltaPnPlot"
				 ,"DeltaPnPerpPlot"
				 ,"DeltaPnPerpxPlot"
				 ,"DeltaPnPerpyPlot"				 				 
				 ,"DeltaPnParPlot"				 				 
				 ,"DeltaPtxPlot"
				 ,"DeltaPtyPlot"
				 //,"APlot"
				 //,"kMissPlot"
				 //,"PMissPlot"
				 //,"PMissMinusPlot" 	 
				 //,"DeltaAlphaT_DeltaPT_0_00To0_20Plot","DeltaAlphaT_DeltaPT_0_20To0_40Plot","DeltaAlphaT_DeltaPT_0_40To1_00Plot"
				 //,"DeltaPT_DeltaAlphaT_0_00To45_00Plot","DeltaPT_DeltaAlphaT_45_00To90_00Plot"
				 //,"DeltaPT_DeltaAlphaT_90_00To135_00Plot","DeltaPT_DeltaAlphaT_135_00To180_00Plot"				 
				 //,"DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_MuonCosTheta_0_00To0_50Plot"
				 //,"DeltaAlphaT_MuonCosTheta_0_50To0_75Plot","DeltaAlphaT_MuonCosTheta_0_75To1_00Plot"
				 //,"DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot"
				 //,"DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot","DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot"				 
				 //,"DeltaPT_MuonCosTheta_Minus1_00To0_00Plot","DeltaPT_MuonCosTheta_0_00To0_50Plot"
				 //,"DeltaPT_MuonCosTheta_0_50To0_75Plot","DeltaPT_MuonCosTheta_0_75To1_00Plot"
				 //,"DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot","DeltaPT_ProtonCosTheta_0_00To0_50Plot"
				 //,"DeltaPT_ProtonCosTheta_0_50To0_75Plot","DeltaPT_ProtonCosTheta_0_75To1_00Plot"				 
				 //,"ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot","ProtonCosTheta_MuonCosTheta_0_00To0_50Plot"
				 //,"ProtonCosTheta_MuonCosTheta_0_50To0_75Plot","ProtonCosTheta_MuonCosTheta_0_75To1_00Plot"			
				 //,"MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot","MuonMomentum_MuonCosTheta_0_00To0_50Plot"
				 //,"MuonMomentum_MuonCosTheta_0_50To0_75Plot","MuonMomentum_MuonCosTheta_0_75To1_00Plot"				 
				 //,"ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot","ProtonMomentum_ProtonCosTheta_0_00To0_50Plot"
				 //,"ProtonMomentum_ProtonCosTheta_0_50To0_75Plot","ProtonMomentum_ProtonCosTheta_0_75To1_00Plot"				 
				 //,"DeltaPhiT_DeltaPT_0_00To0_20Plot","DeltaPhiT_DeltaPT_0_20To0_40Plot","DeltaPhiT_DeltaPT_0_40To1_00Plot"
				 //,"DeltaPn_DeltaPT_0_00To0_20Plot","DeltaPn_DeltaPT_0_20To0_40Plot","DeltaPn_DeltaPT_0_40To1_00Plot"				 
				 //,"DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot","DeltaPty_DeltaPtx_Minus0_15To0_15Plot","DeltaPty_DeltaPtx_0_15To0_55Plot"
				 //,"DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot","DeltaPtx_DeltaPty_Minus0_15To0_15Plot","DeltaPtx_DeltaPty_0_15To0_45Plot"
				 //,"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot"
				 //,"ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot"
				 //,"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot"
				 //,"ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot"
				 //,"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot"
				 //,"ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot","ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot"				 
				 //,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot", "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot"
				 //,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot"
				 //,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot","ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot"
				 //,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot"
				 //,"ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot","ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot"
				 //,"ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot","ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot"
				 //,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot"
				 //,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot"
				 //,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot"
				 //,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot"
				 //,"ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot","ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot"
				 //,"ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot","ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot" 
				 //,"ECal_DeltaAlphaT_0_00To45_00Plot", "ECal_DeltaAlphaT_45_00To90_00Plot"	
				 //,"ECal_DeltaAlphaT_90_00To135_00Plot", "ECal_DeltaAlphaT_135_00To180_00Plot"				 			 
				 ,"SerialDeltaPT_MuonCosThetaPlot_0"
				 ,"SerialDeltaPT_MuonCosThetaPlot_1"
				 ,"SerialDeltaPT_MuonCosThetaPlot_2"
				 ,"SerialDeltaPT_MuonCosThetaPlot_3"	
				 ,"SerialDeltaPT_ProtonCosThetaPlot_0"
				 ,"SerialDeltaPT_ProtonCosThetaPlot_1"
				 ,"SerialDeltaPT_ProtonCosThetaPlot_2"	
				 ,"SerialDeltaPT_ProtonCosThetaPlot_3"
				 ,"SerialMuonMomentum_MuonCosThetaPlot_0"
				 ,"SerialMuonMomentum_MuonCosThetaPlot_1"
				 ,"SerialMuonMomentum_MuonCosThetaPlot_2"
				 ,"SerialMuonMomentum_MuonCosThetaPlot_3"	
				 ,"SerialProtonMomentum_ProtonCosThetaPlot_0"
				 ,"SerialProtonMomentum_ProtonCosThetaPlot_1"
				 ,"SerialProtonMomentum_ProtonCosThetaPlot_2"
				 ,"SerialProtonMomentum_ProtonCosThetaPlot_3"						
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot_0"
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot_1"
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot_2"
				 ,"SerialDeltaAlphaT_MuonCosThetaPlot_3"						
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot_0"
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot_1"
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot_2"
				 ,"SerialDeltaAlphaT_ProtonCosThetaPlot_3"				
				 ,"SerialDeltaAlphaT_DeltaPTPlot_0"
				 ,"SerialDeltaAlphaT_DeltaPTPlot_1"
				 ,"SerialDeltaAlphaT_DeltaPTPlot_2"
				 ,"SerialDeltaAlphaT_MuonMomentumPlot_0"
				 ,"SerialDeltaAlphaT_MuonMomentumPlot_1"
				 ,"SerialDeltaAlphaT_MuonMomentumPlot_2"
				 ,"SerialDeltaAlphaT_ProtonMomentumPlot_0"
				 ,"SerialDeltaAlphaT_ProtonMomentumPlot_1"
				 ,"SerialDeltaAlphaT_ProtonMomentumPlot_2"				 				 
				 ,"SerialDeltaAlpha3Dq_DeltaPnPlot_0"
				 ,"SerialDeltaAlpha3Dq_DeltaPnPlot_1"
				 ,"SerialDeltaAlpha3Dq_DeltaPnPlot_2"
				 ,"SerialDeltaAlpha3DMu_DeltaPnPlot_0"
				 ,"SerialDeltaAlpha3DMu_DeltaPnPlot_1"
				 ,"SerialDeltaAlpha3DMu_DeltaPnPlot_2"				 				 	
				 ,"SerialProtonMomentum_DeltaAlphaTPlot_0"
				 ,"SerialProtonMomentum_DeltaAlphaTPlot_1"
				 ,"SerialProtonMomentum_DeltaAlphaTPlot_2"
				 ,"SerialProtonMomentum_DeltaAlphaTPlot_3"
				 ,"SerialDeltaPT_DeltaAlphaTPlot_0"
				 ,"SerialDeltaPT_DeltaAlphaTPlot_1"
				 ,"SerialDeltaPT_DeltaAlphaTPlot_2"
				 ,"SerialDeltaPT_DeltaAlphaTPlot_3"
				 //,"SerialDeltaPn_DeltaAlphaTPlot_0"
				 //,"SerialDeltaPn_DeltaAlphaTPlot_1"
				 //,"SerialDeltaPn_DeltaAlphaTPlot_2"
				 //,"SerialDeltaPn_DeltaAlphaTPlot_3"
				 //,"SerialDeltaPn_DeltaAlpha3DqPlot_0"
				 //,"SerialDeltaPn_DeltaAlpha3DqPlot_1"
				 //,"SerialDeltaPn_DeltaAlpha3DqPlot_2"
				 //,"SerialDeltaPn_DeltaAlpha3DqPlot_3"
				 //,"SerialDeltaPn_DeltaAlpha3DMuPlot_0"
				 //,"SerialDeltaPn_DeltaAlpha3DMuPlot_1"
				 //,"SerialDeltaPn_DeltaAlpha3DMuPlot_2"
				 //,"SerialDeltaPn_DeltaAlpha3DMuPlot_3"				 				 				 				
				 ,"SerialDeltaPhiT_DeltaPTPlot_0"
				 ,"SerialDeltaPhiT_DeltaPTPlot_1"
				 ,"SerialDeltaPhiT_DeltaPTPlot_2"		
//				 ,"SerialDeltaPn_DeltaPTPlot_0"
//				 ,"SerialDeltaPn_DeltaPTPlot_1"
//				 ,"SerialDeltaPn_DeltaPTPlot_2"		
				 ,"SerialECal_DeltaPTPlot_0"
				 ,"SerialECal_DeltaPTPlot_1"
				 ,"SerialECal_DeltaPTPlot_2"
				 ,"SerialECal_DeltaPtxPlot_0"
				 ,"SerialECal_DeltaPtxPlot_1"
				 ,"SerialECal_DeltaPtxPlot_2"
				 ,"SerialECal_DeltaPtyPlot_0"
				 ,"SerialECal_DeltaPtyPlot_1"
				 ,"SerialECal_DeltaPtyPlot_2"				 				 		 
				 ,"SerialECal_DeltaAlphaTPlot_0"
				 ,"SerialECal_DeltaAlphaTPlot_1"
				 ,"SerialECal_DeltaAlphaTPlot_2"	
				 ,"SerialECal_DeltaAlphaTPlot_3"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_0"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_1"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_2"	
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_3"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_4"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_5"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_6"	
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_7"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_8"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_9"
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_10"	
				 ,"SerialECal_DeltaPTDeltaAlphaTPlot_11"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_0"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_1"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_2"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_3"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_4"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_5"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_6"	
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_7"
				 ,"SerialECal_DeltaPtxDeltaPtyPlot_8"			 				 				 				 	
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_0"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_1"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_2"	
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_3"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_4"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_5"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_6"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_7"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_8"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_9"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_10"
				 ,"SerialECal_MuonCosThetaMuonMomentumPlot_11"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_0"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_1"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_2"	
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_3"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_4"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_5"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_6"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_7"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_8"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_9"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_10"
				 ,"SerialECal_ProtonCosThetaProtonMomentumPlot_11"							
				 ,"SerialProtonCosTheta_MuonCosThetaPlot_0"
				 ,"SerialProtonCosTheta_MuonCosThetaPlot_1"
				 ,"SerialProtonCosTheta_MuonCosThetaPlot_2"
				 ,"SerialProtonCosTheta_MuonCosThetaPlot_3"
				 ,"SerialDeltaPty_DeltaPtxPlot_0"
				 ,"SerialDeltaPty_DeltaPtxPlot_1"
				 ,"SerialDeltaPty_DeltaPtxPlot_2"						
				 ,"SerialDeltaPtx_DeltaPtyPlot_0"
				 ,"SerialDeltaPtx_DeltaPtyPlot_1"
				 ,"SerialDeltaPtx_DeltaPtyPlot_2"
						
	};	

	// --------------------------------------------------------------------------------------------------------------------------

	// POT Normalization & Runs


//	vector<TString> Runs{"Run1"};
/*	vector<TString> Runs{"Run1","Run3"};*/
/*	vector<TString> Runs{"Run1","Run2","Run3"};*/
//	vector<TString> Runs{"Run1","Run2","Run3","Run4a","Run4aRutgers","Combined"};
	vector<TString> Runs{"Run1","Run2","Run3","Run4a","Run4b","Run5","Combined"};
//	vector<TString> Runs{"Run1","Run2","Run3","Run4","Run5"};

	// ------------------------------------ //

	// v43 Run 1 

	static const double tor860_wcut_Run1 = 4.566e+19;
	static const double E1DCNT_wcut_Run1 = 10127594.0;
	static const double EXT_Run1 = 32514217.0;

	static const double Fulltor860_wcut_Run1 = 1.62e20;
	static const double FullE1DCNT_wcut_Run1 = 36139233.0;
	static const double FullEXT_Run1 = 65498807.0;
	
	// ------------------------------------ //	
	
	// v43 Run 2 

	static const double tor860_wcut_Run2 = 0.;
	static const double E1DCNT_wcut_Run2 = 1.;
	static const double EXT_Run2 = 1.;

	static const double Fulltor860_wcut_Run2 = 2.62e20;
	static const double FullE1DCNT_wcut_Run2 = 62045760.0;
	static const double FullEXT_Run2 = 152404890.0;
	
	// ------------------------------------ //	
	
	// v43 Run 3 

	static const double tor860_wcut_Run3 = 9.513e+18;
	static const double E1DCNT_wcut_Run3 = 2299517.0;
	static const double EXT_Run3 = 131356320.0;

	static const double Fulltor860_wcut_Run3 = 2.55e20;
	static const double FullE1DCNT_wcut_Run3 = 61012955.0;
	static const double FullEXT_Run3 = 214555174.0;

        // ------------------------------------ //      
        
        // Run 4a 
        
        static const double tor860_wcut_Run4a = 0.;
        static const double E1DCNT_wcut_Run4a = 0.;
        static const double EXT_Run4a = 1.;
        
        static const double Fulltor860_wcut_Run4a = 5.537e+19;
        static const double FullE1DCNT_wcut_Run4a = 12329277.0;
        static const double FullEXT_Run4a = 35295091.0;
        
        // ------------------------------------ //      
        
        // Run 4b 
        
        static const double tor860_wcut_Run4b = 0.;
        static const double E1DCNT_wcut_Run4b = 0.;
        static const double EXT_Run4b = 1.;
        
        static const double Fulltor860_wcut_Run4b = 1.47e20;
        static const double FullE1DCNT_wcut_Run4b = 34737600.0;
        static const double FullEXT_Run4b = 91167643.0;
        
        // ------------------------------------ //
        
        // Run 4c
        
        static const double tor860_wcut_Run4c = 0.;
        static const double E1DCNT_wcut_Run4c = 0.;
        static const double EXT_Run4c = 1.;
        
        static const double Fulltor860_wcut_Run4c = 9.11e19;
        static const double FullE1DCNT_wcut_Run4c = 20639283.0;
        static const double FullEXT_Run4c = 50978641.0;
        
        // ------------------------------------ //
       
        // Run 4d 
        
        static const double tor860_wcut_Run4d = 0.;
        static const double E1DCNT_wcut_Run4d = 0.;
        static const double EXT_Run4d = 1.;
        
        static const double Fulltor860_wcut_Run4d = 5.11e19;
        static const double FullE1DCNT_wcut_Run4d = 11610448.0;
        static const double FullEXT_Run4d = 75637751.0;
        
        // ------------------------------------ //
        
        // Run 5 
        
        static const double tor860_wcut_Run5 = 0.;
        static const double E1DCNT_wcut_Run5 = 0.;
        static const double EXT_Run5 = 1.;
        
        static const double Fulltor860_wcut_Run5 = 1.552e20;
        static const double FullE1DCNT_wcut_Run5 = 37159286.0;
        static const double FullEXT_Run5 = 115079907.0;
        
	// Combined POT
	// For now just Runs 1 + 2 +3

	static const double Fulltor860_wcut_Combined = Fulltor860_wcut_Run1 + Fulltor860_wcut_Run2 + Fulltor860_wcut_Run3;	
	static const double FullEXT_Combined = FullEXT_Run1 + FullEXT_Run2 + FullEXT_Run3;	
	static const double FullE1DCNT_wcut_Combined = FullE1DCNT_wcut_Run1 + FullE1DCNT_wcut_Run2 + FullE1DCNT_wcut_Run3;	

	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Binning

//	static const int NBinsA = 10; static const double ArrayNBinsA[NBinsA+1] = {0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5}; 
	static const int NBinsA = 11; static const double ArrayNBinsA[NBinsA+1] = {0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7}; 
	
//	static const int NBinskMiss = 13; static const double ArrayNBinskMiss[NBinskMiss+1] = {0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9,1.1}; 
	static const int NBinskMiss = 10; static const double ArrayNBinskMiss[NBinskMiss+1] = {0,0.07,0.14,0.2,0.3,0.4,0.47,0.55,0.65,0.75,0.9};
	
	static const int NBinsPMiss = 10; static const double ArrayNBinsPMiss[NBinsPMiss+1] = {0,0.07,0.14,0.2,0.3,0.4,0.47,0.55,0.65,0.75,0.9};
	
//	static const int NBinsPMissMinus = 14; 
//	static const double ArrayNBinsPMissMinus[NBinsPMissMinus+1] = {0.525,0.6,0.675,0.75,0.825,0.9,0.975,1.05,1.125,1.2,1.275,1.35,1.425,1.5,1.6};	
	static const int NBinsPMissMinus = 11; static const double ArrayNBinsPMissMinus[NBinsPMissMinus+1] = {0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6};
	
	static const int NBinsDeltaPT = 13; static const double ArrayNBinsDeltaPT[NBinsDeltaPT+1] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9};

//	static const int NBinsDeltaPL = 16; static const double ArrayNBinsDeltaPL[NBinsDeltaPL+1] = {-0.47,-0.39,-0.32,-0.25,-0.17,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47};
	static const int NBinsDeltaPL = 10; static const double ArrayNBinsDeltaPL[NBinsDeltaPL+1] = {-0.47,-0.37,-0.27,-0.17,-0.08,0,0.08,0.17,0.27,0.37,0.47};
		
//	static const int NBinsDeltaPn = 10; static const double ArrayNBinsDeltaPn[NBinsDeltaPn+1] = {0,0.07,0.14,0.2,0.27,0.34,0.4,0.47,0.55,0.65,0.75};
	static const int NBinsDeltaPn = 10; static const double ArrayNBinsDeltaPn[NBinsDeltaPn+1] = {0,0.07,0.14,0.2,0.3,0.4,0.47,0.55,0.65,0.75,0.85};

	static const int NBinsDeltaPnPerp = 8; static const double ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp+1] = {0,0.07,0.14,0.2,0.3,0.4,0.47,0.55,0.85};
	
	static const int NBinsDeltaPnPerpx = 11; static const double ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx+1] = {-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55};
	static const int NBinsDeltaPnPerpy = 11; static const double ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy+1] = {-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55};		

	static const int NBinsDeltaPnPar = 10; static const double ArrayNBinsDeltaPnPar[NBinsDeltaPnPar+1] = {-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45};		

//	static const int NBinsDeltaPtx = 20; static const double ArrayNBinsDeltaPtx[NBinsDeltaPtx+1] = {-0.55,-0.47,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55};
	static const int NBinsDeltaPtx = 11; static const double ArrayNBinsDeltaPtx[NBinsDeltaPtx+1] = {-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55};
	
//	static const int NBinsDeltaPty = 14; static const double ArrayNBinsDeltaPty[NBinsDeltaPty+1] = {-0.47,-0.36,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.36,0.47};
	static const int NBinsDeltaPty = 12; static const double ArrayNBinsDeltaPty[NBinsDeltaPty+1] = {-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45};
	
//	static const int NBinsDeltaAlphaT = 10; static const double ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT+1] = { 0.,15.,30.,45.,62.,77.,92.,115.,135.,160.,180. };
	static const int NBinsDeltaAlphaT = 7; static const double ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT+1] = { 0.,22.,44.,66.,88.,110.,145.,180. };
//	static const int NBinsDeltaAlpha3Dq = 8; static const double ArrayNBinsDeltaAlpha3Dq[NBinsDeltaAlpha3Dq+1] = { 0.,30.,50.,70.,90.,110.,130.,155.,180. };
	static const int NBinsDeltaAlpha3Dq = 8; static const double ArrayNBinsDeltaAlpha3Dq[NBinsDeltaAlpha3Dq+1] = { 0., 25., 50., 70., 90., 110., 130., 150.,180. };
	static const int NBinsDeltaAlpha3DMu = 8; static const double ArrayNBinsDeltaAlpha3DMu[NBinsDeltaAlpha3DMu+1] = { 0.,30.,50.,70.,90.,110.,130.,155.,180. };	

//	static const int NBinsDeltaPhiT = 13; static const double ArrayNBinsDeltaPhiT[NBinsDeltaPhiT+1] = {0.,10.,20.,30.,40.,50.,60.,75.,90.,106.,126.,145.,162.,180.};
	static const int NBinsDeltaPhiT = 12; static const double ArrayNBinsDeltaPhiT[NBinsDeltaPhiT+1] = {0.,12.5,25.,37.5,50.,60.,75.,90.,106.,126.,145.,162.,180.};
	
	static const int NBinsDeltaPhi3D = 10; static const double ArrayNBinsDeltaPhi3D[NBinsDeltaPhi3D+1] = {0.,12.5,25.,37.5,50.,60.,75.,90.,106.,126.,180.};	
	
	static const int NBinsEQE = 10; static const double ArrayNBinsEQE[NBinsEQE+1] = {0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.46,1.6}; 	
	
	// STV analysis

	static const int NBinsECal = 9; static const double ArrayNBinsECal[NBinsECal+1] = {0.2,0.35,0.5,0.65,0.8,0.95,1.1,1.25,1.4,1.6}; 

//	static const int NBinsQ2 = 8; static const double ArrayNBinsQ2[NBinsQ2+1] = { 0.0,0.1,0.2,0.3,0.4,0.5,0.65,0.8,1.}; 	
	static const int NBinsQ2 = 8; static const double ArrayNBinsQ2[NBinsQ2+1] = { 0.,0.08,0.18,0.28,0.39,0.5,0.65,0.8,1.}; 	
	
//	static const int NBinsMuonMomentum = 11; static const double ArrayNBinsMuonMomentum[NBinsMuonMomentum+1] = { 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2};
	static const int NBinsMuonMomentum = 10; static const double ArrayNBinsMuonMomentum[NBinsMuonMomentum+1] = { 0.1,0.2,0.3,0.4,0.5,0.64,0.77,0.9,1.,1.1,1.2};


//	static const int NBinsMuonPhi = 17; static const double ArrayNBinsMuonPhi[NBinsMuonPhi+1] = { -180.,-160.,-140.,-120.,-100.,-80,-60,-40,-20.,0.,20.,40.,60.,90.,120.,140.,160.,180.};
	static const int NBinsMuonPhi = 10; static const double ArrayNBinsMuonPhi[NBinsMuonPhi+1] = { -180.,-144.,-108.,-72.,-36.,0.,36.,72.,108.,144.,180.};

//	static const int NBinsMuonCosTheta = 19;
//	static const double ArrayNBinsMuonCosTheta[NBinsMuonCosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.89,0.95,1.};

	static const int NBinsMuonCosTheta = 18;
	static const double ArrayNBinsMuonCosTheta[NBinsMuonCosTheta+1] = { -1.,-0.85,-0.7,-0.57,-0.45,-0.32,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.72,0.84,0.95,1.};

//	static const int NBinsProtonMomentum = 11; static const double ArrayNBinsProtonMomentum[NBinsProtonMomentum+1] = {0.3,0.38,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.87,1.};
	static const int NBinsProtonMomentum = 10; static const double ArrayNBinsProtonMomentum[NBinsProtonMomentum+1] = {0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.};

//	static const int NBinsProtonPhi = 18; static const double ArrayNBinsProtonPhi[NBinsProtonPhi+1] = {-180.,-160.,-140.,-120.,-100.,-80,-60,-40,-20.,0.,20.,40.,60.,80.,100.,120.,140.,160.,180.}; 
	static const int NBinsProtonPhi = 10; static const double ArrayNBinsProtonPhi[NBinsProtonPhi+1] = { -180.,-144.,-108.,-72.,-36.,0.,36.,72.,108.,144.,180.}; 

//	static const int NBinsProtonCosTheta = 12;
//	static const double ArrayNBinsProtonCosTheta[NBinsProtonCosTheta+1] = { -1.,-0.73,-0.43,-0.18,0.05,0.2,0.34,0.47,0.6,0.7,0.8,0.9,1. };
	static const int NBinsProtonCosTheta = 11;
	static const double ArrayNBinsProtonCosTheta[NBinsProtonCosTheta+1] = { -1.,-0.73,-0.43,-0.18,0.05,0.2,0.37,0.54,0.7,0.8,0.9,1.};
	
	// -----------------------------------------------------------------------------------------------------------------------------------------
	
	// MCC8 CCQE analysis binning 
	
	static const int CCQENBinsECal = 7; 
	static const double CCQEArrayNBinsECal[CCQENBinsECal+1] = {0.3, 0.51428571, 0.72857143, 0.94285714, 1.15714286,1.37142857, 1.58571429, 1.8}; 
	static const int CCQENBinsQ2 = 7; 
	static const double CCQEArrayNBinsQ2[CCQENBinsQ2+1] = {0., 0.17142857, 0.34285714, 0.51428571, 0.68571429, 0.85714286,1.02857143,1.2}; 

	static const int CCQENBinsMuonMomentum = 7; // PRL binning
	static const double CCQEArrayNBinsMuonMomentum[CCQENBinsMuonMomentum+1] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3,1.5}; // PRL binnning

	static const int CCQENBinsMuonPhi = 7; 
	static const double CCQEArrayNBinsMuonPhi[CCQENBinsMuonPhi+1] = {-180.,-128.57142857,-77.14285714,-25.71428571,25.71428571,77.14285714,128.57142857,180.};
	static const int CCQENBinsMuonCosTheta = 7; 
	static const double CCQEArrayNBinsMuonCosTheta[CCQENBinsMuonCosTheta+1] = {-0.65, -0.4083,-0.167,0.075,0.317,0.5583, 0.8, 0.95};

	static const int CCQENBinsProtonMomentum = 7; // PRL binning
	static const double CCQEArrayNBinsProtonMomentum[CCQENBinsProtonMomentum+1] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};  // PRL binning

	static const int CCQENBinsProtonPhi = 7; 
	static const double CCQEArrayNBinsProtonPhi[CCQENBinsProtonPhi+1] ={-180.,-128.57142857,-77.14285714,-25.71428571,25.71428571,77.14285714,128.57142857,180.};
	static const int CCQENBinsProtonCosTheta = 7; 
	static const double CCQEArrayNBinsProtonCosTheta[CCQENBinsProtonCosTheta+1] = {0.15,0.2583,0.366,0.475 ,0.583,0.692,0.8,0.95}; 

	// -----------------------------------------------------------------------------------------------------------------------------------------

	static const int NBinsLLRPID = 40; static const double MinLLRPID = -1., MaxLLRPID = 1.;
	static const int NBinsChi2 = 30; static const double MinChi2 = 0., MaxChi2 = 300.;
	static const int NBinsThreePlaneChi2 = 25; static const double MinThreePlaneChi2 = 0., MaxThreePlaneChi2 = 500.;
	static const int NBinsThreePlaneChi2LogLikelihood = 40;
	static const double MinThreePlaneChi2LogLikelihood = -5., MaxThreePlaneChi2LogLikelihood = 5.;
	static const int NBinsDeltaTheta = 18; static const double MinDeltaTheta = 0., MaxDeltaTheta = 180.;
	static const int NBinsDeltaPhi = 15; static const double MinDeltaPhi = 0., MaxDeltaPhi = 360.;
	static const int NBinsdYZ = 20; static const double MindYZ = 0., MaxdYZ = 500.;
	static const int NBinsdXYZ = 30; static const double MindXYZ = 0., MaxdXYZ = 600.;
	static const int NBinsMuonLength = 50; static const double MinMuonLength = 0., MaxMuonLength = 400.;
	static const int NBinsProtonLength = 50; static const double MinProtonLength = 0., MaxProtonLength = 100.;
	static const int NBinsLengthDifference = 35; static const double MinLengthDifference = -150., MaxLengthDifference = 550.;
	static const int NBinsNPE = 50; static const double MinNPE = 0., MaxNPE = 3000.;
	static const int NBinsDistance = 22; static const double MinDistance = 0., MaxDistance = 11.;
	static const int NBinsEv = 36; static const double MinEv = 0.25, MaxEv = 1.6;
	static const int NBinsNu = 30; static const double MinNu = 0., MaxNu = 0.8;
	static const int NBinsNuScore = 25; static const double MinNuScore = 0., MaxNuScore = 1.;
	static const int NBinsFlashScore = 25; static const double MinFlashScore = 0., MaxFlashScore = 50.; 
	static const int NBinsTrackScore = 25; static const double MinTrackScore = 0.5, MaxTrackScore = 1.;
	static const int NBinsMuonVertexDistance = 25; static const double MinMuonVertexDistance = 0., MaxMuonVertexDistance = 2.;
	static const int NBinsProtonVertexDistance = 25; static const double MinProtonVertexDistance = 0., MaxProtonVertexDistance = 5.;
	static const int NBinsVertexActivity = 25; static const double MinVertexActivity = 0., MaxVertexActivity = 1.;
	static const int NBinsVertexX = 15; static const double MinVertexX = 10., MaxVertexX = 246.;
	static const int NBinsVertexY = 15; static const double MinVertexY = -105., MaxVertexY = 105.;
	static const int NBinsVertexZ = 25; static const double MinVertexZ = 10., MaxVertexZ = 1026.;
	
	// -----------------------------------------------------------------------------------------------------------------------------------------------

	// Labels for 1D plots
	
	static TString LabelXAxisECal = ";E^{Cal} [GeV]"; static TString LabelXAxisTrueECal = ";True E^{Cal} [GeV]";
	static TString LabelXAxisEQE = ";E^{QE} [GeV]"; static TString LabelXAxisTrueEQE = ";True E^{QE} [GeV]";
	static TString LabelXAxisQ2 = ";Q^{2} [GeV^{2}/c^{2}]"; static TString LabelXAxisTrueQ2 = ";True Q^{2} [GeV^{2}/c^{2}]";
		
	static TString LabelXAxisMuonMomentum = ";p_{#mu} [GeV/c]"; static TString LabelXAxisTrueMuonMomentum = ";True p_{#mu} [GeV/c]";
	static TString LabelXAxisProtonMomentum = ";p_{p} [GeV/c]"; static TString LabelXAxisTrueProtonMomentum = ";True p_{p} [GeV/c]";
	static TString LabelXAxisMuonPhi = ";#phi_{#mu} [deg]"; static TString LabelXAxisTrueMuonPhi = ";True #phi_{#mu} [deg]";
	static TString LabelXAxisProtonPhi = ";#phi_{p} [deg]"; static TString LabelXAxisTrueProtonPhi = ";True #phi_{p} [deg]";
	static TString LabelXAxisMuonCosTheta = ";cos#theta_{#mu}"; static TString LabelXAxisTrueMuonCosTheta = ";True cos#theta_{#mu}";
	static TString LabelXAxisProtonCosTheta = ";cos#theta_{p}"; static TString LabelXAxisTrueProtonCosTheta = ";True cos#theta_{p}";

	static TString LabelXAxiskMiss = ";k_{Miss} [GeV/c]"; static TString LabelXAxisTruekMiss = ";True k_{miss} [GeV/c]";
	static TString LabelXAxisPMissMinus = ";P^{-}_{Miss} [GeV/c]"; static TString LabelXAxisTruePMissMinus = ";True P^{-}_{miss} [GeV/c]";
	static TString LabelXAxisPMiss = ";P_{Miss} [GeV/c]"; static TString LabelXAxisTruePMiss = ";True P_{miss} [GeV/c]";

	static TString LabelXAxisDeltaPL = ";#deltap_{L} [GeV/c]"; static TString LabelXAxisTrueDeltaPL = ";True #deltap_{L} [GeV/c]";
	static TString LabelXAxisDeltaPn = ";p_{n} [GeV/c]"; static TString LabelXAxisTrueDeltaPn = ";True p_{n} [GeV/c]";
	static TString LabelXAxisDeltaPnPerp = ";p_{n#perp}  [GeV/c]"; static TString LabelXAxisTrueDeltaPnPerp = ";True p_{n#perp}  [GeV/c]";
	static TString LabelXAxisDeltaPnPerpx = ";p_{n#perp ,x}  [GeV/c]"; static TString LabelXAxisTrueDeltaPnPerpx = ";True p_{n#perp ,x}  [GeV/c]";
	static TString LabelXAxisDeltaPnPerpy = ";p_{n#perp ,y}  [GeV/c]"; static TString LabelXAxisTrueDeltaPnPerpy = ";True p_{n#perp ,y}  [GeV/c]";		
	static TString LabelXAxisDeltaPnPar = ";p_{n#parallel} [GeV/c]"; static TString LabelXAxisTrueDeltaPnPar = ";True p_{n#parallel} [GeV/c]";		
	static TString LabelXAxisDeltaPtx = ";#deltap_{T,x} [GeV/c]"; static TString LabelXAxisTrueDeltaPtx = ";True #deltap_{T,x} [GeV/c]";
	static TString LabelXAxisDeltaPty = ";#deltap_{T,y} [GeV/c]"; static TString LabelXAxisTrueDeltaPty = ";True #deltap_{T,y} [GeV/c]";
	static TString LabelXAxisA = ";#alpha [c]"; static TString LabelXAxisTrueA = ";True #alpha [c]";
	
	static TString LabelXAxisDeltaPT = ";#deltap_{T} [GeV/c]"; static TString LabelXAxisTrueDeltaPT = ";True #deltap_{T} [GeV/c]";
	static TString LabelXAxisDeltaAlphaT = ";#delta#alpha_{T} [deg]"; static TString LabelXAxisTrueDeltaAlphaT = ";True #delta#alpha_{T} [deg]";
	static TString LabelXAxisDeltaAlpha3Dq = ";#alpha_{3D} [deg]"; static TString LabelXAxisTrueDeltaAlpha3Dq = ";True #alpha_{3D} [deg]";
	static TString LabelXAxisDeltaAlpha3DMu = ";#alpha_{3D,#mu} [deg]"; static TString LabelXAxisTrueDeltaAlpha3DMu = ";True #alpha_{3D,#mu} [deg]";	
	static TString LabelXAxisDeltaPhiT = ";#delta#phi_{T} [deg]"; static TString LabelXAxisTrueDeltaPhiT = ";True #delta#phi_{T} [deg]";
	static TString LabelXAxisDeltaPhi3D = ";#phi_{3D} [deg]"; static TString LabelXAxisTrueDeltaPhi3D = ";True #phi_{3D} [deg]";	

	TString RecoLabelXAxisLLRPID = ";Candidate LLR PID Score";		
	TString RecoLabelXAxisMuonLLRPID = ";Muon Candidate LLR PID Score";
	TString RecoLabelXAxisProtonLLRPID = ";Proton Candidate LLR PID Score";
	TString RecoLabelXAxisThreePlaneChi2MuonCandidateLogLikelihood = ";3-Plane LogLikelihood #mu-Candidate";
	TString RecoLabelXAxisThreePlaneChi2ProtonCandidateLogLikelihood = ";3-Plane LogLikelihood p-Candidate";
	TString RecoLabelXAxisChi2 = "Y-Plane ;#chi^{2}_{p}";
	TString RecoLabelXAxisThreePlaneChi2 = ";3-plane #chi^{2}_{p}";
	TString RecoLabelXAxisDeltaTheta = ";#delta#theta_{#mu,p}";
	TString RecoLabelXAxisThreePlaneChi2LogLikelihood = ";3-Plane LogLikelihood";
	TString RecoLabelXAxisDeltaPhi = ";#delta#phi_{#mu,p}";
	TString RecoLabelXAxisdYZ = ";d_{YZ} [cm]";
	TString RecoLabelXAxisdXYZ = ";d_{XYZ} [cm]";
	TString RecoLabelXAxisNPE = ";# PE";
	TString RecoLabelXAxisDistance = ";#mu - p distance [cm]";
	TString RecoLabelXAxisNuScore = ";#nu Score";
	TString RecoLabelXAxisEv = ";E_{#nu} [GeV]";
	TString RecoLabelXAxisNu = ";Energy Transfer [GeV]";
	TString RecoLabelXAxisFlashScore = ";Flash Score";
	TString RecoLabelXAxisMuonTrackScore = ";Muon Track Score";
	TString RecoLabelXAxisMuonVertexDistanceTrackScore = ";Muon Vertex Distance [cm]";
	TString RecoLabelXAxisProtonVertexDistanceTrackScore = ";Proton Vertex Distance [cm]";
	TString RecoLabelXAxisVertexActivity = ";Vertex Activity";
	TString RecoLabelXAxisVertexX = ";Vertex x [cm]";
	TString RecoLabelXAxisVertexY = ";Vertex y [cm]";
	TString RecoLabelXAxisVertexZ = ";Vertex z [cm]";
	TString RecoLabelXAxisProtonTrackScore = ";Proton Track Score";
	TString RecoLabelXAxisLengthDifference = ";l_{#mu} - l_{p} [cm]";
	TString RecoLabelXAxisMuonLength = ";l_{#mu} [cm]";
	TString RecoLabelXAxisProtonLength = ";l_{p} [cm]";
	
	// ------------------------------------------------------------------------------------------------------------------------------------

	// Labels for 2D Plots

	static TString LabelXAxisDeltaPT2D = LabelXAxisTrueDeltaPT+";Reco #deltap_{T} [GeV/c]";
	static TString LabelXAxisDeltaAlphaT2D = LabelXAxisTrueDeltaAlphaT+";Reco #delta#alpha_{T} [deg]";
	static TString LabelXAxisDeltaAlpha3Dq2D = LabelXAxisTrueDeltaAlpha3Dq+";Reco #alpha_{3D} [deg]";
	static TString LabelXAxisDeltaAlpha3DMu2D = LabelXAxisTrueDeltaAlpha3DMu+";Reco #alpha_{3D,#mu} [deg]";		
	static TString LabelXAxisDeltaPhiT2D = LabelXAxisTrueDeltaPhiT+";Reco #delta#phi_{T} [deg]";
	static TString LabelXAxisDeltaPhi3D2D = LabelXAxisTrueDeltaPhi3D+";Reco #phi_{3D} [deg]";	

	static TString LabelXAxisDeltaPL2D = LabelXAxisTrueDeltaPL+";Reco #deltap_{L} [GeV/c]";
	static TString LabelXAxisDeltaPn2D = LabelXAxisTrueDeltaPn+";Reco p_{n} [GeV/c]";
	static TString LabelXAxisDeltaPnPerp2D = LabelXAxisTrueDeltaPnPerp+";Reco p_{n#perp}  [GeV/c]";
	static TString LabelXAxisDeltaPnPerpx2D = LabelXAxisTrueDeltaPnPerpx+";Reco p_{n#perp ,x}  [GeV/c]";
	static TString LabelXAxisDeltaPnPerpy2D = LabelXAxisTrueDeltaPnPerpy+";Reco p_{n#perp ,y}  [GeV/c]";		
	static TString LabelXAxisDeltaPnPar2D = LabelXAxisTrueDeltaPnPar+";Reco p_{n#parallel} [GeV/c]";		
	static TString LabelXAxisDeltaPtx2D = LabelXAxisTrueDeltaPtx+";Reco #deltap_{T,x} [GeV/c]";
	static TString LabelXAxisDeltaPty2D = LabelXAxisTrueDeltaPty+";Reco #deltap_{T,y} [GeV/c]";
	static TString LabelXAxisA2D = LabelXAxisTrueA+";Reco #alpha [c]";

	static TString LabelXAxiskMiss2D = LabelXAxisTruekMiss+";Reco k_{Miss} [GeV/c]";
	static TString LabelXAxisPMiss2D = LabelXAxisTruePMiss+";Reco P_{Miss} [GeV/c]";
	static TString LabelXAxisPMissMinus2D = LabelXAxisTruePMissMinus+";Reco P^{-}_{Miss} [GeV/c]";

	static TString LabelXAxisMuonMomentum2D = LabelXAxisTrueMuonMomentum+";Reco P_{#mu} [GeV/c]";
	static TString LabelXAxisMuonPhi2D = LabelXAxisTrueMuonPhi+";Reco #phi_{#mu} [deg]";
	static TString LabelXAxisMuonCosTheta2D = LabelXAxisTrueMuonCosTheta+";Reco cos#theta_{#mu}";

	static TString LabelXAxisProtonMomentum2D = LabelXAxisTrueProtonMomentum+";Reco P_{p} [GeV/c]";
	static TString LabelXAxisProtonPhi2D = LabelXAxisTrueProtonPhi+";Reco #phi_{p} [deg]";
	static TString LabelXAxisProtonCosTheta2D = LabelXAxisTrueProtonCosTheta+";Reco cos#theta_{p}";

	static TString LabelXAxisECal2D = ";True E^{Cal} [GeV];Reco E^{Cal} [GeV]";
	static TString LabelXAxisEQE2D = ";True E^{QE} [GeV];Reco E^{QE} [GeV]";
	static TString LabelXAxisQ22D = ";True Q^{2} [GeV^{2}/c^{2}];Reco Q^{2} [GeV^{2}/c^{2}]";
	
	static TString RecoLabelXAxisVertexX2D = ";True Vertex x [cm];Reco Vertex x [cm]";
	static TString RecoLabelXAxisVertexY2D = ";True Vertex y [cm];Reco Vertex y [cm]";
	static TString RecoLabelXAxisVertexZ2D = ";True Vertex z [cm];Reco Vertex z [cm]";			

	// --------------------------------------------------------------------------------------------------------------------------------------

	static int OverlayColor = kAzure-4;
	static int BeamOnColor = kBlack;
	static int GenieOverlayColor = kBlack;
	static int GenieColor = 610;
	static int Geniev3OutOfTheBoxColor = kMagenta;
//	static int GenieOverlayColor = kOrange+7;	
	static int GiBUUColor = kGreen+1; // 410
//	static int NuWroColor = kOrange+7;	
	static int NuWroColor = kOrange+10;
	static int NEUTColor = kRed+1;
	static int GENIEv2Color = kBlue;	
	static int GENIEv3_0_4_Color = kGray+1;
//	static int SuSav2Color = 429;	
	static int SuSav2Color = kOrange-3;	
	
	static int NormUncBandColor = kGray+3;	
	
	static int G18LineStyle = kSolid;
	
	static int GiBUULineStyle = 2;
	static int NEUTLineStyle = 7;
	static int NuWroLineStyle = 5;
	
	static int UntunedLineStyle = 2;
	static int G21LineStyle = 7;
	static int Gv2LineStyle = 5;							
	
	// ---------------------------------------------------------------------------------------------------------------------------------------

	// Constants, Cuts & Thresholds

	static const double NPE = 10;

	static const double MuonChi2Cut = 100.;
	static const double ProtonChi2Cut = 80.;

	static const double MuonThreePlaneChi2LogLikelihoodCut = -1.;
	static const double ProtonThreePlaneChi2LogLikelihoodCut = -0.5;
	static const double ProtonLLRPIDScore = 0.05;

	static const double CosmicPID = -999.;
	static const int CosmicPdg = -99;

	static const double BE = 0.04; // GeV

	static const double DeltaRCut = 3; // cm
	static const double ThresholdPurity = 0.1; // 10%
	static const double UBSpaceReso = 0.3; // cm // 3mm spacing between the wires

	static const double DeltaThetaCut = 180.; // deg
	static const double DeltaPhiCentralValue = 180.; // deg
	static const double DeltaPhiOpeningAngle = 35.; // deg
	static const double MaxTransMissMomentum = 0.35; // deg

	static const double MaxMuPDistance = 2.; // cm

	static const double MaxPFParticleTrackDistance = 20.; // cm

	static const double YZBeamFlashVertexMaxDistance = 500.; // cm
	static const double BeamFlashPEThreshold = 100; // PEs

	static const double ChargedPionMomentumThres = 0.07;
	static const double NeutralPionMomentumThres = 0.0;

	static const double PurityThreshold = 0.1;

	static const double MinimumNuScore = 0.;

	static double POTUncertainty = 0.02; // 2% POT Uncertainty		

	static double NTargetUncertainty = 0.01; // 1% NTarget Uncertainty

	constexpr float TRACK_SCORE_CUT = 0.5;	

	// Quality cuts

	const TString MuRangeMCSAgree = "0.25"; const double MuRangeMCSAgreeValue = 0.25;
	const TString MuMinMCSLengthThres = "100";  const double MuMinMCSLengthThresValue = 100;
	const TString MuMinSumHits = "0.7";  const double MuMinSumHitsValue = 0.7;
	const TString PMinSumHits = "1.";  const double PMinSumHitsValue = 1.;	

	// --------------------------------------------------------------------------------------------------------------------------------------

	static const double NuScoreStep = (MaxNuScore - MinNuScore) / double(NBinsNuScore);
	static const double LLPStep = (MaxThreePlaneChi2LogLikelihood - MinThreePlaneChi2LogLikelihood) / double(NBinsThreePlaneChi2LogLikelihood);
	static const double MuonLengthStep = (MaxMuonLength - MinMuonLength) / double(NBinsMuonLength);
	static const double LLRPIDStep = (MaxLLRPID - MinLLRPID) / double(NBinsLLRPID);

	// --------------------------------------------------------------------------------------------------------------------------------

	// Selection  & Quality cuts

	// To recalibrate or not, that is the question

	TString Recalibrate = "";
//	TString Recalibrate = "_Recalibrate"; 

	// Reco level

	const TString CC1p = "CC1p && CandidateMu_MCParticle_Pdg == 13 && CandidateP_MCParticle_Pdg == 2212 && NumberPi0 == 0"; 
	
	const TString CC1p1pi = "CandidateMu_MCParticle_Pdg == 211 && CandidateP_MCParticle_Pdg == 2212 && NumberPi0 == 0";	

	const TString Containment = "CandidateMu_StartContainment == 1 && CandidateP_StartContainment == 1 && CandidateP_EndContainment == 1 && True_CandidateMu_StartContainment == 1";

	const TString ContainmentData = "CandidateMu_StartContainment == 1 && CandidateP_StartContainment == 1 && CandidateP_EndContainment == 1";

	const TString ProtonMom = "CandidateP_P_Range" + Recalibrate + " > " + TString(std::to_string(ArrayNBinsProtonMomentum[0])) + " && CandidateP_P_Range" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsProtonMomentum[NBinsProtonMomentum])) + " && True_CandidateP_P > "+ TString(std::to_string(ArrayNBinsProtonMomentum[0]))\
		+ " && True_CandidateP_P < "+ TString(std::to_string(ArrayNBinsProtonMomentum[NBinsProtonMomentum]));

	const TString ProtonMomData = "CandidateP_P_Range" + Recalibrate + " > " + TString(std::to_string(ArrayNBinsProtonMomentum[0])) + " && CandidateP_P_Range" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsProtonMomentum[NBinsProtonMomentum]));
 
	const TString DeltaPhiT = "Reco_DeltaPhiT" + Recalibrate + " > " + TString(std::to_string(ArrayNBinsDeltaPhiT[0])) + " && Reco_DeltaPhiT" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]))\
		+ " && True_DeltaPhiT > " + TString(std::to_string(ArrayNBinsDeltaPhiT[0])) + " && True_DeltaPhiT < " + TString(std::to_string(ArrayNBinsDeltaPhiT[NBinsDeltaPhiT])); 

	const TString DeltaPhiTData = "Reco_DeltaPhiT" + Recalibrate + " > " + TString(std::to_string(ArrayNBinsDeltaPhiT[0])) + " && Reco_DeltaPhiT" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsDeltaPhiT[NBinsDeltaPhiT])); 

	const TString MuonMom = "True_CandidateMu_P > "+ TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && True_CandidateMu_P < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum]))\
		+ " && ( ( CandidateMu_EndContainment == 1 && CandidateMu_P_Range" + Recalibrate + " > "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && CandidateMu_P_Range" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum])) +" ) || ( CandidateMu_EndContainment == 0 && CandidateMu_P_MCS" + Recalibrate + " > "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && CandidateMu_P_MCS" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum])) +" ) )" ;

	const TString MuonMomData = "( ( CandidateMu_EndContainment == 1 && CandidateMu_P_Range" + Recalibrate + " > "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && CandidateMu_P_Range" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum])) +" ) || ( CandidateMu_EndContainment == 0 && CandidateMu_P_MCS" + Recalibrate + " > "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && CandidateMu_P_MCS" + Recalibrate + " < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum])) +" ) )" ;

	const TString QualityCut = " ( (CandidateMu_EndContainment == 1 && TMath::Abs(CandidateMu_P_Range" + Recalibrate\
		+ "-CandidateMu_P_MCS" + Recalibrate + ")/CandidateMu_P_Range" + Recalibrate + " < "\
		+ MuRangeMCSAgree +" ) || (CandidateMu_EndContainment == 0 && CandidateMu_Length > " + MuMinMCSLengthThres + ") )";

	const TString OnlyContainedMuQualityCut = " ( ( CandidateMu_EndContainment == 1 && TMath::Abs(CandidateMu_P_Range" + Recalibrate + "-CandidateMu_P_MCS" + Recalibrate + ")/CandidateMu_P_Range" + Recalibrate + " < "+ MuRangeMCSAgree +" ) || (CandidateMu_EndContainment == 0) )";

	const TString MinHitsMu = "CandidateMu_P_Range > 0";

	const TString MinHitsP = "CandidateP_P_Range > 0";

	const TString MinMuonVertexDist = "CandidateMuStartVertexDistance < CandidateMuEndVertexDistance";

	const TString MinProtonVertexDist = "CandidatePStartVertexDistance < CandidatePEndVertexDistance";

	const TString AvoidFlippedTracks = "StartToStartDistance < EndToEndDistance";
	
	// --------------------------------------------------------------------------------------------------------------------------------
	
	// Qualifiers using reco info

	const TString qualifier = CC1p +" && "+ Containment + " && " + ProtonMom + " && " + DeltaPhiT + " && " + MuonMom + " && " + QualityCut + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;
		
	const TString qualifierNoPCut = CC1p +" && "+ Containment + " && " + DeltaPhiT + " && " + QualityCut + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;

	const TString qualifierPPNoCut = CC1p +" && "+ Containment + " && " + DeltaPhiT + " && " + MuonMom + " && " + QualityCut + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;				

	const TString qualifierData = ContainmentData +" && " + ProtonMomData + " && " + DeltaPhiTData + " && " + MuonMomData + " && " + QualityCut + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;

	const TString qualifierNoMuRangeCut = CC1p +" && "+ Containment +" && " + ProtonMom + " && " + DeltaPhiT + " && " + MuonMom + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;

	const TString qualifierNoHitSumCut = CC1p +" && "+ Containment +" && " + ProtonMom + " && " + DeltaPhiT + " && " + MuonMom + " && " + QualityCut\
		+ " && " + MinMuonVertexDist + " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;

	const TString qualifierNoExitMuQC = CC1p +" && "+ Containment +" && " + ProtonMom + " && " + DeltaPhiT + " && " + MuonMom + " && " + OnlyContainedMuQualityCut\
		+ " && " + MinMuonVertexDist + " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;
		
	const TString qualifierCC1p1pi = CC1p1pi +" && "+ Containment +" && " + ProtonMom + " && " + DeltaPhiT + " && " + MuonMom + " && " + QualityCut + " && " + MinMuonVertexDist\
		+ " && " + MinProtonVertexDist + " && " + AvoidFlippedTracks;		

	//----------------------------------------//
	
	// Qualifiers using truth level info

	const TString TrueCC1p = "CC1p && NumberPi0 == 0";

	const TString TrueProtonMom = "Proton_MCParticle_Mom > " + TString(std::to_string(ArrayNBinsProtonMomentum[0])) + " && Proton_MCParticle_Mom < "\
		+ TString(std::to_string(ArrayNBinsProtonMomentum[NBinsProtonMomentum]));

	const TString TrueMuonMom = "Muon_MCParticle_Mom > " + TString(std::to_string(ArrayNBinsMuonMomentum[0])) + " && Muon_MCParticle_Mom < "\
		+ TString(std::to_string(ArrayNBinsMuonMomentum[NBinsMuonMomentum])); 

	//const TString TrueDeltaPT = "True_Pt > " + TString(std::to_string(ArrayNBinsDeltaPT[0])) + " && True_Pt < " + TString(std::to_string(ArrayNBinsDeltaPT[NBinsDeltaPT]));

	const TString TrueDeltaAlphaT = "True_DeltaAlphaT > " + TString(std::to_string(ArrayNBinsDeltaAlphaT[0])) + " && True_DeltaAlphaT < "\
		+ TString(std::to_string(ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]));

	const TString TrueDeltaPhiT = "True_DeltaPhiT > " + TString(std::to_string(ArrayNBinsDeltaPhiT[0])) + " && True_DeltaPhiT < " + TString(std::to_string(ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]));

	const TString TrueQualifier = TrueCC1p + " && " + TrueProtonMom + " && " + TrueMuonMom + " && " + TrueDeltaAlphaT + " && " + TrueDeltaPhiT;
	
	const TString TrueQualifierNoPCut = TrueCC1p + " && " + TrueDeltaAlphaT + " && " + TrueDeltaPhiT;	
	
	const TString TrueQualifierNoPPCut = TrueCC1p + " && " + TrueMuonMom + " && " + TrueDeltaAlphaT + " && " + TrueDeltaPhiT;	

	//----------------------------------------//
	
	// And moving on towards a 2D analysis
		
	static const int TwoDNBinsDeltaPT = 3; std::vector<double> TwoDArrayNBinsDeltaPT{0.0,0.2,0.4,1.0};
	static const int TwoDNBinsDeltaPn = 3; std::vector<double> TwoDArrayNBinsDeltaPn{0.0,0.2,0.4,1.0};	
	static const int TwoDNBinsDeltaAlphaT = 4; std::vector<double> TwoDArrayNBinsDeltaAlphaT{0.0,45.0,90.0,135.0,180.0};
	static const int TwoDNBinsDeltaAlpha3Dq = 4; std::vector<double> TwoDArrayNBinsDeltaAlpha3Dq{0.0,45.0,90.0,135.0,180.0};
	static const int TwoDNBinsDeltaAlpha3DMu = 4; std::vector<double> TwoDArrayNBinsDeltaAlpha3DMu{0.0,45.0,90.0,135.0,180.0};		
	static const int TwoDNBinsMuonCosTheta = 4; std::vector<double> TwoDArrayNBinsMuonCosTheta{-1.,0.0,0.5,0.75,1.0};
	static const int TwoDNBinsProtonCosTheta = 4; std::vector<double> TwoDArrayNBinsProtonCosTheta{-1.,0.0,0.5,0.75,1.0};
	static const int TwoDNBinsDeltaPtx = 3; std::vector<double> TwoDArrayNBinsDeltaPtx{-0.55,-0.15,0.15,0.55};	
	static const int TwoDNBinsDeltaPty = 3; std::vector<double> TwoDArrayNBinsDeltaPty{-0.75,-0.15,0.15,0.45};
	static const int TwoDNBinsMuonMomentum = 3; std::vector<double> TwoDArrayNBinsMuonMomentum{0.1,0.4,0.6,1.2};
	static const int TwoDNBinsProtonMomentum = 3; std::vector<double> TwoDArrayNBinsProtonMomentum{0.3,0.5,0.7,1.};		

	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInProtonMomentumSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.3 < Pp < 0.5 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.5 < Pp < 0.7 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}  // 0.7 < Pp < 1.0 GeV/c
												
											};

	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInMuonMomentumSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.1 < Pmu < 0.4 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.4 < Pmu < 0.6 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}  // 0.6 < Pmu < 1.2 GeV/c
												
											};
	
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.}, // DeltaPT < 0.2 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}, // 0.2 < DeltaPT < 0.4 GeV/c
												{0.,22.,44.,66.,88.,110.,145.,180.}  // DeltaPT > 0.4 GeV/c
												
											};
												
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices{ 
//												{0.,22.,44.,66.,88.,110.,145.,180.},      // DeltaPn < 0.2 GeV/c
//												{0.,30.,50.,70.,90.,110.,130.,155.,180.}, // 0.2 < DeltaPn < 0.4 GeV/c
//												{0.,30.,50.,70.,90.,110.,130.,150.,180.}  // DeltaPn > 0.4 GeV/c
	
												{0., 35., 70., 92., 114., 136., 158.,180.},      // DeltaPn < 0.2 GeV/c
												{0., 25., 50., 70., 90., 110., 130., 150.,180.}, // 0.2 < DeltaPn < 0.4 GeV/c
												{0., 35., 70., 90., 110., 130., 145.,180.}       // DeltaPn > 0.4 GeV/c
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.},      // DeltaPn < 0.2 GeV/c
												{0.,30.,50.,70.,90.,110.,130.,155.,180.}, // 0.2 < DeltaPn < 0.4 GeV/c
												{0.,30.,50.,70.,90.,110.,130.,150.,180.}  // DeltaPn > 0.4 GeV/c	
											};											
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsProtonMomentumInDeltaAlphaTSlices{ 
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.}, // DeltaAlphaT < 45 deg
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.}, // 45 < DeltaAlphaT < 90 deg
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.}, // 90 < DeltaAlphaT < 135 deg
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.}  // DeltaAlphaT > 135 deg		
											};

	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices{ 
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.9}, // DeltaAlphaT < 45 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.9}, // 45 < DeltaAlphaT < 90 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, // 90 < DeltaAlphaT < 135 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}  // DeltaAlphaT > 135 deg		
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices{ 
												{0,0.1,0.2,0.3,0.4,0.85},                               // DeltaAlpha3Dq < 45 deg
												{0,0.1,0.2,0.3,0.4,0.55,0.85},                          // 45 < DeltaAlpha3Dq < 90 deg
												{0,0.08,0.15,0.23,0.3,0.38,0.45,0.85},                  // 90 < DeltaAlpha3Dq < 135 deg
												{0,0.08,0.15,0.23,0.3,0.4,0.85}                         // DeltaAlpha3Dq > 135 deg		
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices{ 
												{0,0.1,0.2,0.3,0.4,0.85},                               // DeltaAlpha3DMu < 45 deg
												{0,0.1,0.2,0.3,0.4,0.47,0.55,0.85},                     // 45 < DeltaAlpha3DMu < 90 deg
												{0,0.08,0.15,0.23,0.3,0.38,0.45,0.85},                  // 90 < DeltaAlpha3DMu < 135 deg
												{0,0.08,0.15,0.23,0.3,0.4,0.85}                         // DeltaAlpha3DMu > 135 deg		
											};																						
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPnInDeltaAlphaTSlices{ 
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.85}, // DeltaAlphaT < 45 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.85}, // 45 < DeltaAlphaT < 90 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.85}, // 90 < DeltaAlphaT < 135 deg
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.85}  // DeltaAlphaT > 135 deg		
											};																						

	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPhiTInDeltaPTSlices{ 
												{0.,12.5,25.,37.5,180.},
												{0.,20.,40.,60.,80.,100.,120.,140.,160.,180.},
												{0.,25.,50.,75.,105.,145.,180.}	
											};	
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPnInDeltaPTSlices{ 
												{0,0.07,0.14,0.2,0.3,0.4,0.85},
												{0,0.3,0.4,0.47,0.85},
												{0,0.47,0.55,0.65,0.75,0.85}	
											};	
											
	//----------------------------------------//																															
									
	std::vector< std::vector<double> > TwoDArrayNBinsECalInDeltaPTSlices{ 
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6}	
											};
											
	//----------------------------------------//																															
									
	std::vector< std::vector<double> > TwoDArrayNBinsECalInDeltaPtxSlices{ 
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6}	
											};
											
	//----------------------------------------//																															
									
	std::vector< std::vector<double> > TwoDArrayNBinsECalInDeltaPtySlices{ 
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.32,1.6}	
											};																																	
														
	//----------------------------------------//																															
									
	std::vector< std::vector<double> > TwoDArrayNBinsECalInDeltaAlphaTSlices{ 
												{0.2,0.38,0.55,0.73,0.9,1.08,1.23,1.4,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.39,1.6},
												{0.2,0.34,0.48,0.62,0.76,0.9,1.04,1.18,1.39,1.6},						
												{0.2,0.38,0.55,0.73,0.9,1.08,1.25,1.43,1.6}
											};																															
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.},
												{0.,22.,44.,66.,88.,110.,145.,180.},
												{0.,22.,44.,66.,88.,110.,145.,180.},				
												{0.,22.,44.,66.,88.,110.,145.,180.}
											};	
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaAlphaTInProtonCosThetaSlices{ 
												{0.,22.,44.,66.,88.,110.,145.,180.},
												{0.,22.,44.,66.,88.,110.,145.,180.},
												{0.,22.,44.,66.,88.,110.,145.,180.},				
												{0.,22.,44.,66.,88.,110.,145.,180.}
											};											
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPTInMuonCosThetaSlices{ 
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, // -1    < cosθμ < 0
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, //  0    < cosθμ < 0.5
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, //  0.5  < cosθμ < 0.75
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}  //  0.75 < cosθμ < 1
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPTInProtonCosThetaSlices{ 
												{0,0.1,0.2,0.3,0.4,0.55,0.65,0.75,0.9}, // -1    < cosθp < 0
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, //  0    < cosθp < 0.5
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.75,0.9}, //  0.5  < cosθp < 0.75
												{0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.47,0.55,0.65,0.9}  //  0.75 < cosθp < 1
											};																						
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsMuonMomentumInMuonCosThetaSlices{ 
												{ 0.1,0.2,0.3,0.4,1.2},                           // -1    < cosθμ < 0
												{ 0.1,0.2,0.3,0.4,0.5,1.2},                       //  0    < cosθμ < 0.5
												{ 0.1,0.2,0.3,0.4,0.5,0.64,0.77,1.2},             //  0.5  < cosθμ < 0.75
												{ 0.1,0.2,0.3,0.4,0.5,0.64,0.77,0.9,1.,1.1,1.2}   //  0.75 < cosθμ < 1
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsProtonCosThetaInMuonCosThetaSlices{ 
												{ -1.,0.05,0.2,0.37,0.54,0.7,0.8,0.9,1.},                   // -1    < cosθμ < 0
												{ -1.,-0.73,-0.43,-0.18,0.05,0.2,0.37,0.54,0.7,0.8,1.}, //  0    < cosθμ < 0.5
												{ -1.,-0.73,-0.43,-0.18,0.05,0.2,0.37,0.54,0.7,0.8,0.9,1.}, //  0.5  < cosθμ < 0.75
												{ -1.,-0.73,-0.43,-0.18,0.05,0.2,0.37,0.54,0.7,0.8,0.9,1.}  //  0.75 < cosθμ < 1
											};											

	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsProtonMomentumInProtonCosThetaSlices{ 
												{0.3,0.38,0.45,0.5,0.55,0.625,1.},
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,1.},
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.},
												{0.3,0.38,0.45,0.5,0.55,0.625,0.7,0.75,0.8,0.87,1.}
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPtyInDeltaPtxSlices{ 
										{-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45},
										{-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45},
										{-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45}
											};
											
	//----------------------------------------//
	
	std::vector< std::vector<double> > TwoDArrayNBinsDeltaPtxInDeltaPtySlices{ 
										{-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55},
										{-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55},
										{-0.55,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.55}
											};																																	
																																										
											
	//----------------------------------------//	
	
	std::vector< std::vector< std::vector<double> > > TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices{	{ // DeltaPT < 0.2 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // DeltaAphaT < 45 deg
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // 45 < DeltaAphaT < 90 deg 
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // 90 < DeltaAlphaT < 135 deg
														{0.2,0.4,0.6,0.8,1.,1.2,1.6} // 135 < DeltaAlphaT < 180 deg
														},	
														
														{ // 0.2 < DeltaPT < 0.4 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // DeltaAphaT < 45 deg
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // 45 < DeltaAphaT < 90 deg
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // 90 < DeltaAlphaT < 135 deg		
														{0.2,0.4,0.6,0.8,1.,1.2,1.6} // 135 < DeltaAlphaT < 180 deg
														},
														
														{ // DeltaPT > 0.4 GeV/c
														{0.2,0.5,0.8,1.2,1.6}, // DeltaAphaT < 45 deg
														{0.2,0.5,0.8,1.2,1.6}, // 45 < DeltaAphaT < 90 deg
														{0.2,0.5,0.8,1.2,1.6}, // 90 < DeltaAlphaT < 135 deg		
														{0.2,0.5,0.8,1.2,1.6} // 135 < DeltaAlphaT < 180 deg
														}		
																
													};
													
	std::vector< std::vector< std::vector<double> > > TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices{	{ // DeltaPtx < -0.15 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // DeltaPty < -0.15 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // -0.15 < DeltaPty < 0.15 GeV/c 
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}  // DeltaPty > 0.15 GeV/c
														},	
														
														{ // -0.15 < DeltaPtx < 0.15 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // DeltaPty < -0.15 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}, // -0.15 < DeltaPty < 0.15 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6}  // DeltaPty > 0.15 GeV/c
														},
														
														{ // DeltaPtx > 0.15 GeV/c
														{0.2,0.5,0.8,1.2,1.6}, // DeltaPty < -0.15 GeV/c
														{0.2,0.5,0.8,1.2,1.6}, // -0.15 < DeltaPty < 0.15 GeV/c
														{0.2,0.6,1.1,1.6}  // DeltaPty > 0.15 GeV/c
														}		
																
													};													
													
	//----------------------------------------//	
	
	std::vector< std::vector< std::vector<double> > > TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices{
														{ // -1 < cosθμ < 0
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,1.6},      // 0.40 < Pμ < 0.60 GeV/c
														{0.2,1.6}           // 0.60 < Pμ < 1.20 GeV/c
														},	
														
														{ // 0 < cosθμ < 0.5
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,0.8,1.6},  // 0.40 < Pμ < 0.6 GeV/c
														{0.2,1.6}   // 0.6 < Pμ < 1.20 GeV/c
														},
														
														{ // 0.5 < cosθμ < 0.75
														{0.2,0.4,0.6,1.6},  // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.6,0.8,1.6},  // 0.40 < Pμ < 0.6 GeV/c
														{0.2,0.8,1.,1.6}    // 0.6 < Pμ < 1.20 GeV/c
														},														
														
														{ // 0.75 < cosθμ < 1
														{0.2,0.5,1.6},      // 0.10 < Pμ < 0.40 GeV/c
														{0.2,0.8,1.6},      // 0.40 < Pμ < 0.6 GeV/c
														{0.2,0.8,1.2,1.6}   // 0.6 < Pμ < 1.20 GeV/c
														}		
																
													};
													
	//----------------------------------------//	
	
	std::vector< std::vector< std::vector<double> > > TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices{
														{ // -1 < cosθp < 0
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.30 < Pp < 0.50 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.50 < Pp < 0.70 GeV/c
														{0.2,1.6}                         // 0.70 < Pp < 1.00 GeV/c
														},	
														
														{ // 0 < cosθp < 0.5
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.30 < Pp < 0.50 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.50 < Pp < 0.70 GeV/c
														{0.2,0.6,0.9,1.2,1.6}             // 0.70 < Pp < 1.00 GeV/c
														},
														
														{ // 0.5 < cosθp < 0.75
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.30 < Pp < 0.50 GeV/c
														{0.2,0.4,0.6,0.8,1.,1.2,1.6},     // 0.50 < Pp < 0.70 GeV/c
														{0.2,0.6,0.8,1.,1.2,1.6}  // 0.70 < Pp < 1.00 GeV/c
														},														
														
														{ // 0.75 < cosθμ < 1
														{0.2,0.8,1.6},                    // 0.30 < Pp < 0.50 GeV/c
														{0.2,0.5,0.8,1.2,1.6},            // 0.50 < Pp < 0.70 GeV/c
														{0.2,0.5,0.8,1.2,1.6}             // 0.70 < Pp < 1.00 GeV/c
														}		
																
													};																										
			
	//----------------------------------------//	

	// Scaling factor for multi dimensional analysis

	static std::map<TString,double> MultiDimScaleFactor =
	{
		{ "DeltaPTPlot", 1. },
		{ "DeltaAlphaTPlot", 1. },
		{ "DeltaAlpha3DqPlot", 1. },
		{ "DeltaAlpha3DMuPlot", 1. },				
		{ "DeltaPhiTPlot", 1. },
		{ "DeltaPhi3DPlot", 1. },		
		{ "MuonMomentumPlot", 1. },
		{ "MuonCosThetaPlot", 1. },
		{ "MuonCosThetaSingleBinPlot", 1. },
		{ "MuonPhiPlot", 1. },
		{ "ProtonMomentumPlot", 1. },
		{ "ProtonCosThetaPlot", 1. },
		{ "ProtonPhiPlot", 1. },
		{ "ECalPlot", 1. },
		{ "ECal_DeltaPT_0_00To0_20Plot", TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) },	
		{ "ECal_DeltaPT_0_20To0_40Plot", TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) },	
		{ "ECal_DeltaPT_0_40To1_00Plot", TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) },
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15Plot", TwoDArrayNBinsDeltaPtx.at(1) - TwoDArrayNBinsDeltaPtx.at(0) },
		{ "ECal_DeltaPtx_Minus0_15To0_15Plot", TwoDArrayNBinsDeltaPtx.at(2) - TwoDArrayNBinsDeltaPtx.at(1) },
		{ "ECal_DeltaPtx_0_15To0_55Plot", TwoDArrayNBinsDeltaPtx.at(3) - TwoDArrayNBinsDeltaPtx.at(2) },
		{ "ECal_DeltaPty_Minus0_75ToMinus0_15Plot", TwoDArrayNBinsDeltaPty.at(1) - TwoDArrayNBinsDeltaPty.at(0) },
		{ "ECal_DeltaPty_Minus0_15To0_15Plot", TwoDArrayNBinsDeltaPty.at(2) - TwoDArrayNBinsDeltaPty.at(1) },
		{ "ECal_DeltaPty_0_15To0_45Plot", TwoDArrayNBinsDeltaPty.at(3) - TwoDArrayNBinsDeltaPty.at(2) },				
		{ "EQEPlot", 1. },
		{ "Q2Plot", 1. },	
		{ "CCQEMuonMomentumPlot", 1. },
		{ "CCQEMuonCosThetaPlot", 1. },
		{ "CCQEProtonMomentumPlot", 1. },
		{ "CCQEProtonCosThetaPlot", 1. },
		{ "CCQEECalPlot", 1. },
		{ "CCQEQ2Plot", 1. },	
		{ "DeltaPLPlot", 1. },
		{ "DeltaPnPlot", 1. },
		{ "DeltaPnPerpPlot", 1. },
		{ "DeltaPnPerpxPlot", 1. },	
		{ "DeltaPnPerpyPlot", 1. },			
		{ "DeltaPnParPlot", 1. },				
		{ "DeltaPtxPlot", 1. },
		{ "DeltaPtyPlot", 1. },	
		{ "APlot", 1. },	
		{ "kMissPlot", 1. },
		{ "PMissPlot", 1. },	
		{ "PMissMinusPlot", 1. },
		{ "DeltaAlphaT_DeltaPT_0_00To0_20Plot", TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) },
		{ "DeltaAlphaT_DeltaPT_0_20To0_40Plot", TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) },
		{ "DeltaAlphaT_DeltaPT_0_40To1_00Plot", TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) },
		{ "DeltaAlpha3Dq_DeltaPn_0_00To0_20Plot", TwoDArrayNBinsDeltaPn.at(1) - TwoDArrayNBinsDeltaPn.at(0) },
		{ "DeltaAlpha3Dq_DeltaPn_0_20To0_40Plot", TwoDArrayNBinsDeltaPn.at(2) - TwoDArrayNBinsDeltaPn.at(1) },
		{ "DeltaAlpha3Dq_DeltaPn_0_40To1_00Plot", TwoDArrayNBinsDeltaPn.at(3) - TwoDArrayNBinsDeltaPn.at(2) },
		{ "DeltaAlpha3DMu_DeltaPn_0_00To0_20Plot", TwoDArrayNBinsDeltaPn.at(1) - TwoDArrayNBinsDeltaPn.at(0) },
		{ "DeltaAlpha3DMu_DeltaPn_0_20To0_40Plot", TwoDArrayNBinsDeltaPn.at(2) - TwoDArrayNBinsDeltaPn.at(1) },
		{ "DeltaAlpha3DMu_DeltaPn_0_40To1_00Plot", TwoDArrayNBinsDeltaPn.at(3) - TwoDArrayNBinsDeltaPn.at(2) },				
		{ "DeltaAlphaT_MuonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) },
		{ "DeltaAlphaT_MuonCosTheta_0_00To0_50Plot", TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) },
		{ "DeltaAlphaT_MuonCosTheta_0_50To0_75Plot", TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) },
		{ "DeltaAlphaT_MuonCosTheta_0_75To1_00Plot", TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) },
		{ "DeltaAlphaT_MuonMomentum_0_10To0_40Plot", TwoDArrayNBinsMuonMomentum.at(1) - TwoDArrayNBinsMuonMomentum.at(0) },
		{ "DeltaAlphaT_MuonMomentum_0_40To0_60Plot", TwoDArrayNBinsMuonMomentum.at(2) - TwoDArrayNBinsMuonMomentum.at(1) },
		{ "DeltaAlphaT_MuonMomentum_0_60To1_20Plot", TwoDArrayNBinsMuonMomentum.at(3) - TwoDArrayNBinsMuonMomentum.at(2) },		
		{ "DeltaAlphaT_ProtonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) },
		{ "DeltaAlphaT_ProtonCosTheta_0_00To0_50Plot", TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) },
		{ "DeltaAlphaT_ProtonCosTheta_0_50To0_75Plot", TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) },
		{ "DeltaAlphaT_ProtonCosTheta_0_75To1_00Plot", TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) },
		{ "DeltaAlphaT_ProtonMomentum_0_30To0_50Plot", TwoDArrayNBinsProtonMomentum.at(1) - TwoDArrayNBinsProtonMomentum.at(0) },
		{ "DeltaAlphaT_ProtonMomentum_0_50To0_70Plot", TwoDArrayNBinsProtonMomentum.at(2) - TwoDArrayNBinsProtonMomentum.at(1) },
		{ "DeltaAlphaT_ProtonMomentum_0_70To1_00Plot", TwoDArrayNBinsProtonMomentum.at(3) - TwoDArrayNBinsProtonMomentum.at(2) },		
		{ "DeltaPT_MuonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) },
		{ "DeltaPT_MuonCosTheta_0_00To0_50Plot", TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) },
		{ "DeltaPT_MuonCosTheta_0_50To0_75Plot", TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) },
		{ "DeltaPT_MuonCosTheta_0_75To1_00Plot", TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) },
		{ "DeltaPT_ProtonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) },
		{ "DeltaPT_ProtonCosTheta_0_00To0_50Plot", TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) },
		{ "DeltaPT_ProtonCosTheta_0_50To0_75Plot", TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) },
		{ "DeltaPT_ProtonCosTheta_0_75To1_00Plot", TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) },
		{ "ProtonMomentum_DeltaAlphaT_0_00To45_00Plot", TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) },
		{ "ProtonMomentum_DeltaAlphaT_45_00To90_00Plot", TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) },
		{ "ProtonMomentum_DeltaAlphaT_90_00To135_00Plot", TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) },
		{ "ProtonMomentum_DeltaAlphaT_135_00To180_00Plot", TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) },
		{ "DeltaPT_DeltaAlphaT_0_00To45_00Plot", TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) },
		{ "DeltaPT_DeltaAlphaT_45_00To90_00Plot", TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) },
		{ "DeltaPT_DeltaAlphaT_90_00To135_00Plot", TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) },
		{ "DeltaPT_DeltaAlphaT_135_00To180_00Plot", TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) },
		{ "DeltaPn_DeltaAlphaT_0_00To45_00Plot", TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) },
		{ "DeltaPn_DeltaAlphaT_45_00To90_00Plot", TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) },
		{ "DeltaPn_DeltaAlphaT_90_00To135_00Plot", TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) },
		{ "DeltaPn_DeltaAlphaT_135_00To180_00Plot", TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) },	
		{ "DeltaPn_DeltaAlpha3Dq_0_00To45_00Plot", TwoDArrayNBinsDeltaAlpha3Dq.at(1) - TwoDArrayNBinsDeltaAlpha3Dq.at(0) },
		{ "DeltaPn_DeltaAlpha3Dq_45_00To90_00Plot", TwoDArrayNBinsDeltaAlpha3Dq.at(2) - TwoDArrayNBinsDeltaAlpha3Dq.at(1) },
		{ "DeltaPn_DeltaAlpha3Dq_90_00To135_00Plot", TwoDArrayNBinsDeltaAlpha3Dq.at(3) - TwoDArrayNBinsDeltaAlpha3Dq.at(2) },
		{ "DeltaPn_DeltaAlpha3Dq_135_00To180_00Plot", TwoDArrayNBinsDeltaAlpha3Dq.at(4) - TwoDArrayNBinsDeltaAlpha3Dq.at(3) },
		{ "DeltaPn_DeltaAlpha3DMu_0_00To45_00Plot", TwoDArrayNBinsDeltaAlpha3DMu.at(1) - TwoDArrayNBinsDeltaAlpha3DMu.at(0) },
		{ "DeltaPn_DeltaAlpha3DMu_45_00To90_00Plot", TwoDArrayNBinsDeltaAlpha3DMu.at(2) - TwoDArrayNBinsDeltaAlpha3DMu.at(1) },
		{ "DeltaPn_DeltaAlpha3DMu_90_00To135_00Plot", TwoDArrayNBinsDeltaAlpha3DMu.at(3) - TwoDArrayNBinsDeltaAlpha3DMu.at(2) },
		{ "DeltaPn_DeltaAlpha3DMu_135_00To180_00Plot", TwoDArrayNBinsDeltaAlpha3DMu.at(4) - TwoDArrayNBinsDeltaAlpha3DMu.at(3) },					
		{ "ECal_DeltaAlphaT_0_00To45_00Plot", TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) },
		{ "ECal_DeltaAlphaT_45_00To90_00Plot", TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) },
		{ "ECal_DeltaAlphaT_90_00To135_00Plot", TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) },
		{ "ECal_DeltaAlphaT_135_00To180_00Plot", TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) },						
		{ "MuonMomentum_MuonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) },
		{ "MuonMomentum_MuonCosTheta_0_00To0_50Plot", TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) },
		{ "MuonMomentum_MuonCosTheta_0_50To0_75Plot", TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) },
		{ "MuonMomentum_MuonCosTheta_0_75To1_00Plot", TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) },
		{ "ProtonCosTheta_MuonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) },
		{ "ProtonCosTheta_MuonCosTheta_0_00To0_50Plot", TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) },
		{ "ProtonCosTheta_MuonCosTheta_0_50To0_75Plot", TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) },
		{ "ProtonCosTheta_MuonCosTheta_0_75To1_00Plot", TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) },		
		{ "ProtonMomentum_ProtonCosTheta_Minus1_00To0_00Plot", TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) },
		{ "ProtonMomentum_ProtonCosTheta_0_00To0_50Plot", TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) },
		{ "ProtonMomentum_ProtonCosTheta_0_50To0_75Plot", TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) },
		{ "ProtonMomentum_ProtonCosTheta_0_75To1_00Plot", TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) },	
		{ "DeltaPhiT_DeltaPT_0_00To0_20Plot", TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) },
		{ "DeltaPhiT_DeltaPT_0_20To0_40Plot", TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) },				
		{ "DeltaPhiT_DeltaPT_0_40To1_00Plot", TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) },	
//		{ "DeltaPn_DeltaPT_0_00To0_20Plot", TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) },
//		{ "DeltaPn_DeltaPT_0_20To0_40Plot", TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) },				
//		{ "DeltaPn_DeltaPT_0_40To1_00Plot", TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) },
		{ "DeltaPty_DeltaPtx_Minus0_55ToMinus0_15Plot", TwoDArrayNBinsDeltaPtx.at(1) - TwoDArrayNBinsDeltaPtx.at(0) },
		{ "DeltaPty_DeltaPtx_Minus0_15To0_15Plot", TwoDArrayNBinsDeltaPtx.at(2) - TwoDArrayNBinsDeltaPtx.at(1) },
		{ "DeltaPty_DeltaPtx_0_15To0_55Plot", TwoDArrayNBinsDeltaPtx.at(3) - TwoDArrayNBinsDeltaPtx.at(2) },
		{ "DeltaPtx_DeltaPty_Minus0_75ToMinus0_15Plot", TwoDArrayNBinsDeltaPty.at(1) - TwoDArrayNBinsDeltaPty.at(0) },
		{ "DeltaPtx_DeltaPty_Minus0_15To0_15Plot", TwoDArrayNBinsDeltaPty.at(2) - TwoDArrayNBinsDeltaPty.at(1) },
		{ "DeltaPtx_DeltaPty_0_15To0_45Plot",  TwoDArrayNBinsDeltaPty.at(3) - TwoDArrayNBinsDeltaPty.at(2) },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_0_00To45_00Plot", ( TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) ) * ( TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) ) },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_45_00To90_00Plot", ( TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) ) * ( TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) ) },
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_90_00To135_00Plot", ( TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) ) * ( TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) ) },	
		{ "ECal_DeltaPT_0_00To0_20_DeltaAlphaT_135_00To180_00Plot", ( TwoDArrayNBinsDeltaPT.at(1) - TwoDArrayNBinsDeltaPT.at(0) ) * ( TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) ) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_0_00To45_00Plot", ( TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) )  * ( TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) ) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_45_00To90_00Plot", ( TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) ) * ( TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) ) },
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_90_00To135_00Plot", ( TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) ) * ( TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) ) },	
		{ "ECal_DeltaPT_0_20To0_40_DeltaAlphaT_135_00To180_00Plot", ( TwoDArrayNBinsDeltaPT.at(2) - TwoDArrayNBinsDeltaPT.at(1) ) * ( TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) ) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_0_00To45_00Plot", ( TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) )  * ( TwoDArrayNBinsDeltaAlphaT.at(1) - TwoDArrayNBinsDeltaAlphaT.at(0) ) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_45_00To90_00Plot", ( TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) ) * ( TwoDArrayNBinsDeltaAlphaT.at(2) - TwoDArrayNBinsDeltaAlphaT.at(1) ) },
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_90_00To135_00Plot", ( TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) ) * ( TwoDArrayNBinsDeltaAlphaT.at(3) - TwoDArrayNBinsDeltaAlphaT.at(2) ) },	
		{ "ECal_DeltaPT_0_40To1_00_DeltaAlphaT_135_00To180_00Plot",  ( TwoDArrayNBinsDeltaPT.at(3) - TwoDArrayNBinsDeltaPT.at(2) ) * ( TwoDArrayNBinsDeltaAlphaT.at(4) - TwoDArrayNBinsDeltaAlphaT.at(3) ) },	
		
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_75ToMinus0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(1) - TwoDArrayNBinsDeltaPtx.at(0) ) * ( TwoDArrayNBinsDeltaPty.at(1) - TwoDArrayNBinsDeltaPty.at(0) ) },	
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_Minus0_15To0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(1) - TwoDArrayNBinsDeltaPtx.at(0) ) * ( TwoDArrayNBinsDeltaPty.at(2) - TwoDArrayNBinsDeltaPty.at(1) ) },	
		{ "ECal_DeltaPtx_Minus0_55ToMinus0_15_DeltaPty_0_15To0_45Plot", ( TwoDArrayNBinsDeltaPtx.at(1) - TwoDArrayNBinsDeltaPtx.at(0) ) * ( TwoDArrayNBinsDeltaPty.at(3) - TwoDArrayNBinsDeltaPty.at(2) ) },
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_75ToMinus0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(2) - TwoDArrayNBinsDeltaPtx.at(1) ) * ( TwoDArrayNBinsDeltaPty.at(1) - TwoDArrayNBinsDeltaPty.at(0) ) },	
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_Minus0_15To0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(2) - TwoDArrayNBinsDeltaPtx.at(1) ) * ( TwoDArrayNBinsDeltaPty.at(2) - TwoDArrayNBinsDeltaPty.at(1) ) },	
		{ "ECal_DeltaPtx_Minus0_15To0_15_DeltaPty_0_15To0_45Plot", ( TwoDArrayNBinsDeltaPtx.at(2) - TwoDArrayNBinsDeltaPtx.at(1) ) * ( TwoDArrayNBinsDeltaPty.at(3) - TwoDArrayNBinsDeltaPty.at(2) ) },
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_75ToMinus0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(3) - TwoDArrayNBinsDeltaPtx.at(2) ) * ( TwoDArrayNBinsDeltaPty.at(1) - TwoDArrayNBinsDeltaPty.at(0) ) },	
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_Minus0_15To0_15Plot", ( TwoDArrayNBinsDeltaPtx.at(3) - TwoDArrayNBinsDeltaPtx.at(2) ) * ( TwoDArrayNBinsDeltaPty.at(2) - TwoDArrayNBinsDeltaPty.at(1) ) },	
		{ "ECal_DeltaPtx_0_15To0_55_DeltaPty_0_15To0_45Plot", ( TwoDArrayNBinsDeltaPtx.at(3) - TwoDArrayNBinsDeltaPtx.at(2) ) * ( TwoDArrayNBinsDeltaPty.at(3) - TwoDArrayNBinsDeltaPty.at(2) ) },											
			
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_10To0_40Plot",  ( TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) ) * ( TwoDArrayNBinsMuonMomentum.at(1) - TwoDArrayNBinsMuonMomentum.at(0) ) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_10To0_40Plot",  ( TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) ) * ( TwoDArrayNBinsMuonMomentum.at(1) - TwoDArrayNBinsMuonMomentum.at(0) ) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_10To0_40Plot", ( TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) ) * ( TwoDArrayNBinsMuonMomentum.at(1) - TwoDArrayNBinsMuonMomentum.at(0) ) },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_10To0_40Plot", ( TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) ) * ( TwoDArrayNBinsMuonMomentum.at(1) - TwoDArrayNBinsMuonMomentum.at(0) ) },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_40To0_60Plot", ( TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) ) * ( TwoDArrayNBinsMuonMomentum.at(2) - TwoDArrayNBinsMuonMomentum.at(1) ) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_40To0_60Plot", ( TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) ) * ( TwoDArrayNBinsMuonMomentum.at(2) - TwoDArrayNBinsMuonMomentum.at(1) ) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_40To0_60Plot", ( TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) ) * ( TwoDArrayNBinsMuonMomentum.at(2) - TwoDArrayNBinsMuonMomentum.at(1) ) },
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_40To0_60Plot",  ( TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) ) * ( TwoDArrayNBinsMuonMomentum.at(2) - TwoDArrayNBinsMuonMomentum.at(1) ) },
		{ "ECal_MuonCosTheta_Minus1_00To0_00_MuonMomentum_0_60To1_20Plot", ( TwoDArrayNBinsMuonCosTheta.at(1) - TwoDArrayNBinsMuonCosTheta.at(0) ) * ( TwoDArrayNBinsMuonMomentum.at(3) - TwoDArrayNBinsMuonMomentum.at(2) ) },
		{ "ECal_MuonCosTheta_0_00To0_50_MuonMomentum_0_60To1_20Plot", ( TwoDArrayNBinsMuonCosTheta.at(2) - TwoDArrayNBinsMuonCosTheta.at(1) ) * ( TwoDArrayNBinsMuonMomentum.at(3) - TwoDArrayNBinsMuonMomentum.at(2) ) },
		{ "ECal_MuonCosTheta_0_50To0_75_MuonMomentum_0_60To1_20Plot",  ( TwoDArrayNBinsMuonCosTheta.at(3) - TwoDArrayNBinsMuonCosTheta.at(2) ) * ( TwoDArrayNBinsMuonMomentum.at(3) - TwoDArrayNBinsMuonMomentum.at(2) )},
		{ "ECal_MuonCosTheta_0_75To1_00_MuonMomentum_0_60To1_20Plot",  ( TwoDArrayNBinsMuonCosTheta.at(4) - TwoDArrayNBinsMuonCosTheta.at(3) ) * ( TwoDArrayNBinsMuonMomentum.at(3) - TwoDArrayNBinsMuonMomentum.at(2) ) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_30To0_50Plot", ( TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) ) * ( TwoDArrayNBinsProtonMomentum.at(1) - TwoDArrayNBinsProtonMomentum.at(0) ) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_30To0_50Plot", ( TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) ) * ( TwoDArrayNBinsProtonMomentum.at(1) - TwoDArrayNBinsProtonMomentum.at(0) ) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_30To0_50Plot", ( TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) ) * ( TwoDArrayNBinsProtonMomentum.at(1) - TwoDArrayNBinsProtonMomentum.at(0) ) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_30To0_50Plot",  ( TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) ) * ( TwoDArrayNBinsProtonMomentum.at(1) - TwoDArrayNBinsProtonMomentum.at(0) ) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_50To0_70Plot", ( TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) ) * ( TwoDArrayNBinsProtonMomentum.at(2) - TwoDArrayNBinsProtonMomentum.at(1) ) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_50To0_70Plot", ( TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) ) * ( TwoDArrayNBinsProtonMomentum.at(2) - TwoDArrayNBinsProtonMomentum.at(1) ) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_50To0_70Plot", ( TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) ) * ( TwoDArrayNBinsProtonMomentum.at(2) - TwoDArrayNBinsProtonMomentum.at(1) ) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_50To0_70Plot", ( TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) ) * ( TwoDArrayNBinsProtonMomentum.at(2) - TwoDArrayNBinsProtonMomentum.at(1) ) },
		{ "ECal_ProtonCosTheta_Minus1_00To0_00_ProtonMomentum_0_70To1_00Plot", ( TwoDArrayNBinsProtonCosTheta.at(1) - TwoDArrayNBinsProtonCosTheta.at(0) ) * ( TwoDArrayNBinsProtonMomentum.at(3) - TwoDArrayNBinsProtonMomentum.at(2) ) },
		{ "ECal_ProtonCosTheta_0_00To0_50_ProtonMomentum_0_70To1_00Plot", ( TwoDArrayNBinsProtonCosTheta.at(2) - TwoDArrayNBinsProtonCosTheta.at(1) ) * ( TwoDArrayNBinsProtonMomentum.at(3) - TwoDArrayNBinsProtonMomentum.at(2) ) },
		{ "ECal_ProtonCosTheta_0_50To0_75_ProtonMomentum_0_70To1_00Plot", ( TwoDArrayNBinsProtonCosTheta.at(3) - TwoDArrayNBinsProtonCosTheta.at(2) ) * ( TwoDArrayNBinsProtonMomentum.at(3) - TwoDArrayNBinsProtonMomentum.at(2) ) },
		{ "ECal_ProtonCosTheta_0_75To1_00_ProtonMomentum_0_70To1_00Plot", ( TwoDArrayNBinsProtonCosTheta.at(4) - TwoDArrayNBinsProtonCosTheta.at(3) ) * ( TwoDArrayNBinsProtonMomentum.at(3) - TwoDArrayNBinsProtonMomentum.at(2) ) },
		{ "SerialDeltaPT_MuonCosThetaPlot", 1. },
		{ "SerialDeltaPT_ProtonCosThetaPlot", 1. },
		{ "SerialDeltaPT_DeltaAlphaTPlot", 1. },
		{ "SerialDeltaPn_DeltaAlpha3DqPlot", 1. },
		{ "SerialDeltaPn_DeltaAlpha3DMuPlot", 1. },				
		{ "SerialProtonMomentum_DeltaAlphaTPlot", 1. },		
		{ "SerialMuonMomentum_MuonCosThetaPlot", 1. },
		{ "SerialProtonMomentum_ProtonCosThetaPlot", 1. },
		{ "SerialDeltaAlphaT_MuonCosThetaPlot", 1. },
		{ "SerialDeltaAlphaT_MuonMomentumPlot", 1. },									
		{ "SerialDeltaAlphaT_ProtonCosThetaPlot", 1. },
		{ "SerialDeltaAlphaT_ProtonMomentumPlot", 1. },		
		{ "SerialDeltaAlphaT_DeltaPTPlot", 1. },
		{ "SerialDeltaAlpha3Dq_DeltaPnPlot", 1. },
		{ "SerialDeltaAlpha3DMu_DeltaPnPlot", 1. },				
		{ "SerialDeltaPhiT_DeltaPTPlot", 1. },
//		{ "SerialDeltaPn_DeltaPTPlot", 1. },
		{ "SerialDeltaPn_DeltaAlphaTPlot", 1. },		
		{ "SerialECal_DeltaPTPlot", 1. },
		{ "SerialECal_DeltaPtxPlot", 1. },		
		{ "SerialECal_DeltaPtyPlot", 1. },				
		{ "SerialECal_DeltaAlphaTPlot", 1. },				
		{ "SerialProtonCosTheta_MuonCosThetaPlot", 1. },
		{ "SerialDeltaPty_DeltaPtxPlot", 1. },
		{ "SerialDeltaPtx_DeltaPtyPlot", 1. },
		{ "SerialECal_DeltaPTDeltaAlphaTPlot", 1. },
		{ "SerialECal_DeltaPtxDeltaPtyPlot", 1. },				
		{ "SerialECal_MuonCosThetaMuonMomentumPlot", 1. },
		{ "SerialECal_ProtonCosThetaProtonMomentumPlot", 1. },
		{ "VertexXPlot", 1. },
		{ "VertexYPlot", 1. },
		{ "VertexZPlot", 1. },																
	};	
	
	//----------------------------------------//
	
	// Interaction labels
	
	//const std::vector<int> InteBreakColors{kBlack,kBlue-5,kYellow+1,kOrange+7,kRed+1,kBlue};
	const std::vector<int> InteBreakColors{kBlack,kAzure-4,kOrange-3,kGreen+1,kRed+1,kBlue};		
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};
	const int NInte = InteractionLabels.size();
	
	const std::vector<std::string> knobList = {"AGKYpT1pi_UBGenie","AGKYxF1pi_UBGenie","AhtBY_UBGenie","AxFFCCQEshape_UBGenie","BhtBY_UBGenie",
                 "CV1uBY_UBGenie","CV2uBY_UBGenie","DecayAngMEC_UBGenie","EtaNCEL_UBGenie","FrAbs_N_UBGenie",
                 "FrAbs_pi_UBGenie","FrCEx_N_UBGenie","FrCEx_pi_UBGenie","FrInel_N_UBGenie","FrInel_pi_UBGenie",
                 "FrPiProd_N_UBGenie","FrPiProd_pi_UBGenie","FracDelta_CCMEC_UBGenie","FracPN_CCMEC_UBGenie","MFP_N_UBGenie",
                 "MFP_pi_UBGenie","MaCCQE_UBGenie","MaCCRES_UBGenie","MaNCEL_UBGenie","MaNCRES_UBGenie",
                 "MvCCRES_UBGenie","MvNCRES_UBGenie","NonRESBGvbarnCC1pi_UBGenie","NonRESBGvbarnCC2pi_UBGenie","NonRESBGvbarnNC1pi_UBGenie",
                 "NonRESBGvbarnNC2pi_UBGenie","NonRESBGvbarpCC1pi_UBGenie","NonRESBGvbarpCC2pi_UBGenie","NonRESBGvbarpNC1pi_UBGenie","NonRESBGvbarpNC2pi_UBGenie",
                 "NonRESBGvnCC1pi_UBGenie","NonRESBGvnCC2pi_UBGenie","NonRESBGvnNC1pi_UBGenie","NonRESBGvnNC2pi_UBGenie","NonRESBGvpCC1pi_UBGenie",
                 "NonRESBGvpCC2pi_UBGenie","NonRESBGvpNC1pi_UBGenie","NonRESBGvpNC2pi_UBGenie","NormCCMEC_UBGenie","NormNCMEC_UBGenie",
                 "RDecBR1eta_UBGenie","RDecBR1gamma_UBGenie","RPA_CCQE_UBGenie","Theta_Delta2Npi_UBGenie","TunedCentralValue_UBGenie",
                 "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","splines_general_Spline"};	
	
	//----------------------------------------//		

}
#endif
