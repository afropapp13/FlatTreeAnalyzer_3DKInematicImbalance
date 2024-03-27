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

#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Constants.h"
#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.h"
#include "/exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.h"

using namespace std;
using namespace Constants;

//Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h, double SF = 1.);

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

        // Output root & txt files

	// Root

	TString FileNameAndPath = "OutputFiles/" + fweights + "FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	// Txt

	ofstream OutputTxt;
	OutputTxt.open("OutputFiles/EventDump_" + fOutputFile + ".txt");

	//----------------------------------------//

	// Plot declaration

	// Post FSI

	// 1D Fine Binning

	TH1D* TrueFineBinMuonCosThetaPlot[NInte];
	TH1D* TrueFineBinDeltaPTPlot[NInte];
	TH1D* TrueFineBinDeltaPnPlot[NInte];
	TH1D* TrueFineBinDeltaAlpha3DMuPlot[NInte];
	TH1D* TrueFineBinDeltaAlpha3DqPlot[NInte];
	TH1D* TrueFineBinDeltaAlphaTPlot[NInte];
	TH1D* TrueFineBinDeltaPhiTPlot[NInte];
	TH1D* TrueFineBinDeltaPhi3DPlot[NInte];
	TH1D* TrueFineBinDeltaPnPerpPlot[NInte];
	TH1D* TrueFineBinDeltaPnPerpxPlot[NInte];
	TH1D* TrueFineBinDeltaPnPerpyPlot[NInte];
	TH1D* TrueFineBinDeltaPtxPlot[NInte];
	TH1D* TrueFineBinDeltaPtyPlot[NInte];
	TH1D* TrueFineBinDeltaPnParPlot[NInte];

	TH2D* TrueFineBinECal2DPlot[NInte];
	TH2D* TrueFineBinq2DPlot[NInte];

	// 2D Fine Binning Uncorrelated

	TH1D* TrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
	TH1D* TrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaAlpha3Dq];
	TH1D* TrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[NInte][TwoDNBinsDeltaAlpha3DMu];
	TH1D* TrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
	TH1D* TrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];

	// 2D Fine Binning Correlated

	TH1D* SerialTrueFineBinDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueFineBinDeltaPn_InDeltaAlpha3DqPlot[NInte];
	TH1D* SerialTrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[NInte];
	TH1D* SerialTrueFineBinDeltaAlphaT_InDeltaPTPlot[NInte];
	TH1D* SerialTrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[NInte];
	TH1D* SerialTrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[NInte];

	// 1D Nominal Binning

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueDeltaPTPlot[NInte];
	TH1D* TrueDeltaPnPlot[NInte];
	TH1D* TrueDeltaAlpha3DMuPlot[NInte];
	TH1D* TrueDeltaAlpha3DqPlot[NInte];
	TH1D* TrueDeltaAlphaTPlot[NInte];
	TH1D* TrueDeltaPhiTPlot[NInte];
	TH1D* TrueDeltaPhi3DPlot[NInte];
	TH1D* TrueDeltaPnPerpPlot[NInte];
	TH1D* TrueDeltaPnPerpxPlot[NInte];
	TH1D* TrueDeltaPnPerpyPlot[NInte];
	TH1D* TrueDeltaPtxPlot[NInte];
	TH1D* TrueDeltaPtyPlot[NInte];
	TH1D* TrueDeltaPnParPlot[NInte];

	// 2D Nominal Binning Uncorrelated

	TH1D* TrueDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
	TH1D* TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaAlpha3Dq];
	TH1D* TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[NInte][TwoDNBinsDeltaAlpha3DMu];
	TH1D* TrueDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
	TH1D* TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
	TH1D* TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];

	// 2D Nominal Binning Correlated

	TH1D* SerialTrueDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlpha3DqPlot[NInte];
	TH1D* SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[NInte];
	TH1D* SerialTrueDeltaAlphaT_InDeltaPTPlot[NInte];
	TH1D* SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[NInte];
	TH1D* SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[NInte];

	// 3D Nominal Binning Uncorrelated

	TH1D* TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InDeltaPtxDeltaPtyTwoDPlot[NInte][TwoDNBinsDeltaPT][TwoDNBinsDeltaAlphaT];
	TH1D* TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[NInte][TwoDNBinsMuonCosTheta][TwoDNBinsMuonMomentum];
	TH1D* TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[NInte][TwoDNBinsProtonCosTheta][TwoDNBinsProtonMomentum];

	// 3D Nominal Binning Correlated
	
	TH1D* SerialTrueECal_InDeltaPTDeltaAlphaTPlot[NInte];
	TH1D* SerialTrueECal_InDeltaPtxDeltaPtyPlot[NInte];
	TH1D* SerialTrueECal_InMuonCosThetaMuonMomentumPlot[NInte];
	TH1D* SerialTrueECal_InProtonCosThetaProtonMomentumPlot[NInte];

	//----------------------------------------//

	// Pre FSI

	// 1D Fine Binning

	TH1D* NoFSITrueFineBinMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPTPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPnPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaAlpha3DMuPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaAlpha3DqPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaAlphaTPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPhiTPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPhi3DPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPnPerpPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPnPerpxPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPnPerpyPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPtxPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPtyPlot[NInte];
	TH1D* NoFSITrueFineBinDeltaPnParPlot[NInte];

	TH2D* NoFSITrueFineBinECal2DPlot[NInte];
	TH2D* NoFSITrueFineBinq2DPlot[NInte];

	// 2D Fine Binning Uncorrelated

        TH1D* NoFSITrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
        TH1D* NoFSITrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaAlpha3Dq];
        TH1D* NoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[NInte][TwoDNBinsDeltaAlpha3DMu];
        TH1D* NoFSITrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
        TH1D* NoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
        TH1D* NoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];

	// 2D Fine Binning Correlated

	TH1D* SerialNoFSITrueFineBinDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DqPlot[NInte];
	TH1D* SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[NInte];
	TH1D* SerialNoFSITrueFineBinDeltaAlphaT_InDeltaPTPlot[NInte];
	TH1D* SerialNoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[NInte];
	TH1D* SerialNoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[NInte];

	// 1D Nominal Binning

	TH1D* NoFSITrueMuonCosThetaPlot[NInte];
	TH1D* NoFSITrueDeltaPTPlot[NInte];
	TH1D* NoFSITrueDeltaPnPlot[NInte];
	TH1D* NoFSITrueDeltaAlpha3DMuPlot[NInte];
	TH1D* NoFSITrueDeltaAlpha3DqPlot[NInte];
	TH1D* NoFSITrueDeltaAlphaTPlot[NInte];
	TH1D* NoFSITrueDeltaPhiTPlot[NInte];
	TH1D* NoFSITrueDeltaPhi3DPlot[NInte];
	TH1D* NoFSITrueDeltaPnPerpPlot[NInte];
	TH1D* NoFSITrueDeltaPnPerpxPlot[NInte];
	TH1D* NoFSITrueDeltaPnPerpyPlot[NInte];
	TH1D* NoFSITrueDeltaPtxPlot[NInte];
	TH1D* NoFSITrueDeltaPtyPlot[NInte];
	TH1D* NoFSITrueDeltaPnParPlot[NInte];

	// 2D Nominal Binning Uncorrelated

        TH1D* NoFSITrueDeltaPT_InDeltaAlphaTTwoDPlot[NInte][TwoDNBinsDeltaAlphaT];
        TH1D* NoFSITrueDeltaPn_InDeltaAlpha3DqTwoDPlot[NInte][TwoDNBinsDeltaAlpha3Dq];
        TH1D* NoFSITrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[NInte][TwoDNBinsDeltaAlpha3DMu];
        TH1D* NoFSITrueDeltaAlphaT_InDeltaPTTwoDPlot[NInte][TwoDNBinsDeltaPT];
        TH1D* NoFSITrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];
        TH1D* NoFSITrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[NInte][TwoDNBinsDeltaPn];

	// 2D Nominal Binning Correlated

	TH1D* SerialNoFSITrueDeltaPT_InDeltaAlphaTPlot[NInte];
	TH1D* SerialNoFSITrueDeltaPn_InDeltaAlpha3DqPlot[NInte];
	TH1D* SerialNoFSITrueDeltaPn_InDeltaAlpha3DMuPlot[NInte];
	TH1D* SerialNoFSITrueDeltaAlphaT_InDeltaPTPlot[NInte];
	TH1D* SerialNoFSITrueDeltaAlpha3Dq_InDeltaPnPlot[NInte];
	TH1D* SerialNoFSITrueDeltaAlpha3DMu_InDeltaPnPlot[NInte];

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  // Post FSI

	  // 1D Fine Binning

	  TrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinMuonCosThetaPlot",";cos(#theta_{#mu})",20,-1.,1.);
	  TrueFineBinDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPTPlot",LabelXAxisDeltaPT,25,ArrayNBinsDeltaPT[0],ArrayNBinsDeltaPT[NBinsDeltaPT]);
	  TrueFineBinDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPnPlot",LabelXAxisDeltaPn,25,ArrayNBinsDeltaPn[0],ArrayNBinsDeltaPn[NBinsDeltaPn]);
	  TrueFineBinDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaAlpha3DMuPlot",LabelXAxisDeltaAlpha3DMu,25,ArrayNBinsDeltaAlpha3DMu[0],ArrayNBinsDeltaAlpha3DMu[NBinsDeltaAlpha3DMu]);
	  TrueFineBinDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,25,ArrayNBinsDeltaAlpha3Dq[0],ArrayNBinsDeltaAlpha3Dq[NBinsDeltaAlpha3Dq]);
	  TrueFineBinDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,25,ArrayNBinsDeltaAlphaT[0],ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]);
	  TrueFineBinDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPhiTPlot",LabelXAxisDeltaPhiT,40,ArrayNBinsDeltaPhiT[0],ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]);
	  TrueFineBinDeltaPhi3DPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,40,ArrayNBinsDeltaPhi3D[0],ArrayNBinsDeltaPhi3D[NBinsDeltaPhi3D]);
	  TrueFineBinDeltaPnPerpPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPnPerpPlot",LabelXAxisDeltaPnPerp,30,ArrayNBinsDeltaPnPerp[0],ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]);
	  TrueFineBinDeltaPnPerpxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPnPerpxPlot",LabelXAxisDeltaPnPerpx,30,ArrayNBinsDeltaPnPerpx[0],ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]);
	  TrueFineBinDeltaPnPerpyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPnPerpyPlot",LabelXAxisDeltaPnPerpy,30,ArrayNBinsDeltaPnPerpy[0],ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]);
	  TrueFineBinDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPtxPlot",LabelXAxisDeltaPtx,30,ArrayNBinsDeltaPtx[0],ArrayNBinsDeltaPtx[NBinsDeltaPtx]);
	  TrueFineBinDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPtyPlot",LabelXAxisDeltaPty,30,ArrayNBinsDeltaPty[0],ArrayNBinsDeltaPty[NBinsDeltaPty]);
	  TrueFineBinDeltaPnParPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueFineBinDeltaPnParPlot",LabelXAxisDeltaPnPar,30,-0.7,0.5);

	  TrueFineBinECal2DPlot[inte] = new TH2D(InteractionLabels[inte]+"TrueFineBinECal2DPlot",";E_{#nu}^{true} [GeV];E_{Cal} [GeV]",25,0.,1.5,20,0.,1.5);
	  TrueFineBinq2DPlot[inte] = new TH2D(InteractionLabels[inte]+"TrueFineBinq2DPlot",";q^{true} [GeV];q^{reco} [GeV/c]",25,0.,1.5,20,0.,1.5);

	  // 1D Nominal Binning

	  TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	  TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	  TrueDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn);
	  TrueDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlpha3DMuPlot",LabelXAxisDeltaAlpha3DMu,NBinsDeltaAlpha3DMu,ArrayNBinsDeltaAlpha3DMu);
	  TrueDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq);
	  TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	  TrueDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT);
	  TrueDeltaPhi3DPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D);
	  TrueDeltaPnPerpPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPerpPlot",LabelXAxisDeltaPnPerp,NBinsDeltaPnPerp,ArrayNBinsDeltaPnPerp);
	  TrueDeltaPnPerpxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPerpxPlot",LabelXAxisDeltaPnPerpx,NBinsDeltaPnPerpx,ArrayNBinsDeltaPnPerpx);
	  TrueDeltaPnPerpyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnPerpyPlot",LabelXAxisDeltaPnPerpy,NBinsDeltaPnPerpy,ArrayNBinsDeltaPnPerpy);
	  TrueDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx);
	  TrueDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty);
	  TrueDeltaPnParPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPnParPlot",LabelXAxisDeltaPnPar,NBinsDeltaPnPar,ArrayNBinsDeltaPnPar);

	  // 2D & 3D Fine & Nominal Binning Uncorrelated

	  for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

	    TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";

	    TrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,20,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][ TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1 ]);

	    TrueDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

	  }

	  for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {

	    TString DeltaPnTwoDInDeltaAlpha3DqLabel = "DeltaPn_DeltaAlpha3Dq_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq+1])+"Plot";

	    TrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaPnTwoDInDeltaAlpha3DqLabel,LabelXAxisDeltaPn,20,TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][0],TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][ TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq].size()-1 ]);

	    TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlpha3DqLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][0]);

	  }

	  for (int WhichDeltaAlpha3DMu = 0; WhichDeltaAlpha3DMu < TwoDNBinsDeltaAlpha3DMu; WhichDeltaAlpha3DMu++) {

	    TString DeltaPnTwoDInDeltaAlpha3DMuLabel = "DeltaPn_DeltaAlpha3DMu_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu+1])+"Plot";

	    TrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaPnTwoDInDeltaAlpha3DMuLabel,LabelXAxisDeltaPn,20,TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][0],TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][ TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu].size()-1 ]);

	    TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu] = new TH1D(InteractionLabels[inte]+"True"+DeltaPnTwoDInDeltaAlpha3DMuLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][0]);

	  }

	  for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

	    TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";

	    TrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,20,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][ TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1 ]);

	    TrueDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

	    for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

	      TString ECalTwoDInDeltaPTDeltaAlphaTLabel = "ECal_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";
	      TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[inte][WhichDeltaPT][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPTDeltaAlphaTLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices[WhichDeltaPT][WhichDeltaAlphaT][0]);

	    }

	  }

	  for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

	    TString DeltaAlpha3DqTwoDInDeltaPnLabel = "DeltaAlpha3Dq_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";

	    TrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaAlpha3DqTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3Dq,20,TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][0],TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][ TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn].size()-1 ]);

	    TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlpha3DqTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3Dq,TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][0]);

	    TString DeltaAlpha3DMuTwoDInDeltaPnLabel = "DeltaAlpha3DMu_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";

	    TrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"TrueFineBin"+DeltaAlpha3DMuTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3DMu,20,TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][0],TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][ TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn].size()-1  ]);

	    TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"True"+DeltaAlpha3DMuTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3DMu,TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][0]);

	  }

	  for (int WhichMuonCosTheta = 0; WhichMuonCosTheta < TwoDNBinsMuonCosTheta; WhichMuonCosTheta++) {

	    for (int WhichMuonMomentum = 0; WhichMuonMomentum < TwoDNBinsMuonMomentum; WhichMuonMomentum++) {

	      TString ECalTwoDInMuonCosThetaMuonMomentumLabel = "ECal_MuonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonCosTheta[WhichMuonCosTheta+1])+"_MuonMomentum_"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsMuonMomentum[WhichMuonMomentum+1])+"Plot";
	      TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[inte][WhichMuonCosTheta][WhichMuonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInMuonCosThetaMuonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum].size()-1,&TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices[WhichMuonCosTheta][WhichMuonMomentum][0]);

	    }

	  }

	  for (int WhichProtonCosTheta = 0; WhichProtonCosTheta < TwoDNBinsProtonCosTheta; WhichProtonCosTheta++) {

	    for (int WhichProtonMomentum = 0; WhichProtonMomentum < TwoDNBinsProtonMomentum; WhichProtonMomentum++) {

	      TString ECalTwoDInProtonCosThetaProtonMomentumLabel = "ECal_ProtonCosTheta_"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonCosTheta[WhichProtonCosTheta+1])+"_ProtonMomentum_"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum])+"To"+tools.ConvertToString(TwoDArrayNBinsProtonMomentum[WhichProtonMomentum+1])+"Plot";
	      TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[inte][WhichProtonCosTheta][WhichProtonMomentum] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInProtonCosThetaProtonMomentumLabel,LabelXAxisECal,TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum].size()-1,&TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices[WhichProtonCosTheta][WhichProtonMomentum][0]);

	    }

	  }

	  for (int WhichDeltaPtx = 0; WhichDeltaPtx < TwoDNBinsDeltaPtx; WhichDeltaPtx++) {

	    for (int WhichDeltaPty = 0; WhichDeltaPty < TwoDNBinsDeltaPty; WhichDeltaPty++) {

	      TString ECalTwoDInDeltaPtxDeltaPtyLabel = "ECal_DeltaPtx_"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPtx[WhichDeltaPtx+1])+"_DeltaPty_"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPty[WhichDeltaPty+1])+"Plot";
	      TrueECal_InDeltaPtxDeltaPtyTwoDPlot[inte][WhichDeltaPtx][WhichDeltaPty] = new TH1D(InteractionLabels[inte]+"True"+ECalTwoDInDeltaPtxDeltaPtyLabel,LabelXAxisECal,TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[WhichDeltaPtx][WhichDeltaPty].size()-1,&TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices[WhichDeltaPtx][WhichDeltaPty][0]);

	    }

	  }

	  // 2D Fine Binning Correlated

	  SerialTrueFineBinDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
	  SerialTrueFineBinDeltaPn_InDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaPn_DeltaAlpha3DqPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices)[0]);
	  SerialTrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaPn_DeltaAlpha3DMuPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices)[0]);
	  SerialTrueFineBinDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
	  SerialTrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaAlpha3Dq_DeltaPnPlot",LabelXAxisDeltaAlpha3Dq,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices)[0]);
	  SerialTrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueFineBinDeltaAlpha3DMu_DeltaPnPlot",LabelXAxisDeltaAlpha3DMu,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices)[0]);

	  // 2D Nominal Binning Correlated

	  SerialTrueDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
	  SerialTrueDeltaPn_InDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaPn_DeltaAlpha3DqPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices)[0]);
	  SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaPn_DeltaAlpha3DMuPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices)[0]);
	  SerialTrueDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
	  SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaAlpha3Dq_DeltaPnPlot",LabelXAxisDeltaAlpha3Dq,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices)[0]);
	  SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueDeltaAlpha3DMu_DeltaPnPlot",LabelXAxisDeltaAlpha3DMu,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices)[0]);

	  // 3D Nominal Binning Correlated
	  
	  SerialTrueECal_InDeltaPTDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueECal_DeltaPTDeltaAlphaTPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices)[0]);
	  SerialTrueECal_InDeltaPtxDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueECal_DeltaPtxDeltaPtyPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices)[0]);
	  SerialTrueECal_InMuonCosThetaMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueECal_MuonCosThetaMuonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices)[0]);
	  SerialTrueECal_InProtonCosThetaProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialTrueECal_ProtonCosThetaProtonMomentumPlot",LabelXAxisECal,tools.Return3DNBins(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices),&tools.Return3DBinIndices(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices)[0]);


	  //--------------------------------------------------//

	  // Pre FSI

	  // 1D Fine Binning

	  NoFSITrueFineBinMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinMuonCosThetaPlot",";cos(#theta_{#mu})",20,-1.,1.);
	  NoFSITrueFineBinDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPTPlot",LabelXAxisDeltaPT,25,ArrayNBinsDeltaPT[0],ArrayNBinsDeltaPT[NBinsDeltaPT]);
	  NoFSITrueFineBinDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPnPlot",LabelXAxisDeltaPn,25,ArrayNBinsDeltaPn[0],ArrayNBinsDeltaPn[NBinsDeltaPn]);
	  NoFSITrueFineBinDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaAlpha3DMuPlot",LabelXAxisDeltaAlpha3DMu,25,ArrayNBinsDeltaAlpha3DMu[0],ArrayNBinsDeltaAlpha3DMu[NBinsDeltaAlpha3DMu]);
	  NoFSITrueFineBinDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,25,ArrayNBinsDeltaAlpha3Dq[0],ArrayNBinsDeltaAlpha3Dq[NBinsDeltaAlpha3Dq]);
	  NoFSITrueFineBinDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,25,ArrayNBinsDeltaAlphaT[0],ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]);
	  NoFSITrueFineBinDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPhiTPlot",LabelXAxisDeltaPhiT,40,ArrayNBinsDeltaPhiT[0],ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]);
	  NoFSITrueFineBinDeltaPhi3DPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,40,ArrayNBinsDeltaPhi3D[0],ArrayNBinsDeltaPhi3D[NBinsDeltaPhi3D]);
	  NoFSITrueFineBinDeltaPnPerpPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPnPerpPlot",LabelXAxisDeltaPnPerp,30,ArrayNBinsDeltaPnPerp[0],ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]);
	  NoFSITrueFineBinDeltaPnPerpxPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPnPerpxPlot",LabelXAxisDeltaPnPerpx,30,ArrayNBinsDeltaPnPerpx[0],ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]);
	  NoFSITrueFineBinDeltaPnPerpyPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPnPerpyPlot",LabelXAxisDeltaPnPerpy,30,ArrayNBinsDeltaPnPerpy[0],ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]);
	  NoFSITrueFineBinDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPtxPlot",LabelXAxisDeltaPtx,30,ArrayNBinsDeltaPtx[0],ArrayNBinsDeltaPtx[NBinsDeltaPtx]);
	  NoFSITrueFineBinDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPtyPlot",LabelXAxisDeltaPty,30,ArrayNBinsDeltaPty[0],ArrayNBinsDeltaPty[NBinsDeltaPty]);
	  NoFSITrueFineBinDeltaPnParPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBinDeltaPnParPlot",LabelXAxisDeltaPnPar,30,-0.7,0.5);

	  NoFSITrueFineBinECal2DPlot[inte] = new TH2D(InteractionLabels[inte]+"NoFSITrueFineBinECal2DPlot",";E_{#nu}^{true} [GeV];E_{Cal} [GeV]",25,0.,1.5,20,0.,1.5);
	  NoFSITrueFineBinq2DPlot[inte] = new TH2D(InteractionLabels[inte]+"NoFSITrueFineBinq2DPlot",";q^{true} [GeV];q^{reco} [GeV/c]",25,0.,1.5,20,0.,1.5);

	  // 1D Nominal Binning

	  NoFSITrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueMuonCosThetaPlot",";cos(#theta_{#mu})",NBinsMuonCosTheta,ArrayNBinsMuonCosTheta[0],ArrayNBinsMuonCosTheta[NBinsMuonCosTheta]);
	  NoFSITrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPTPlot",LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT[0],ArrayNBinsDeltaPT[NBinsDeltaPT]);
	  NoFSITrueDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPnPlot",LabelXAxisDeltaPn,NBinsDeltaPn,ArrayNBinsDeltaPn[0],ArrayNBinsDeltaPn[NBinsDeltaPn]);
	  NoFSITrueDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaAlpha3DMuPlot",LabelXAxisDeltaAlpha3DMu,NBinsDeltaAlpha3DMu,ArrayNBinsDeltaAlpha3DMu[0],ArrayNBinsDeltaAlpha3DMu[NBinsDeltaAlpha3DMu]);
	  NoFSITrueDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaAlpha3DqPlot",LabelXAxisDeltaAlpha3Dq,NBinsDeltaAlpha3Dq,ArrayNBinsDeltaAlpha3Dq[0],ArrayNBinsDeltaAlpha3Dq[NBinsDeltaAlpha3Dq]);
	  NoFSITrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaAlphaTPlot",LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT[0],ArrayNBinsDeltaAlphaT[NBinsDeltaAlphaT]);
	  NoFSITrueDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPhiTPlot",LabelXAxisDeltaPhiT,NBinsDeltaPhiT,ArrayNBinsDeltaPhiT[0],ArrayNBinsDeltaPhiT[NBinsDeltaPhiT]);
	  NoFSITrueDeltaPhi3DPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPhi3DPlot",LabelXAxisDeltaPhi3D,NBinsDeltaPhi3D,ArrayNBinsDeltaPhi3D[0],ArrayNBinsDeltaPhi3D[NBinsDeltaPhi3D]);
	  NoFSITrueDeltaPnPerpPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPnPerpPlot",LabelXAxisDeltaPnPerp,NBinsDeltaPnPerp,ArrayNBinsDeltaPnPerp[0],ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]);
	  NoFSITrueDeltaPnPerpxPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPnPerpxPlot",LabelXAxisDeltaPnPerpx,NBinsDeltaPnPerpx,ArrayNBinsDeltaPnPerpx[0],ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]);
	  NoFSITrueDeltaPnPerpyPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPnPerpyPlot",LabelXAxisDeltaPnPerpy,NBinsDeltaPnPerpy,ArrayNBinsDeltaPnPerpy[0],ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]);
	  NoFSITrueDeltaPtxPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPtxPlot",LabelXAxisDeltaPtx,NBinsDeltaPtx,ArrayNBinsDeltaPtx[0],ArrayNBinsDeltaPtx[NBinsDeltaPtx]);
	  NoFSITrueDeltaPtyPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPtyPlot",LabelXAxisDeltaPty,NBinsDeltaPty,ArrayNBinsDeltaPty[0],ArrayNBinsDeltaPty[NBinsDeltaPty]);
	  NoFSITrueDeltaPnParPlot[inte] = new TH1D(InteractionLabels[inte]+"NoFSITrueDeltaPnParPlot",LabelXAxisDeltaPnPar,NBinsDeltaPnPar,ArrayNBinsDeltaPnPar[0],ArrayNBinsDeltaPnPar[NBinsDeltaPnPar]);

	  // 2D Fine & Nominal Binning Uncorrelated

	  for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

	    TString DeltaPTTwoDInDeltaAlphaTLabel = "DeltaPT_DeltaAlphaT_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlphaT[WhichDeltaAlphaT+1])+"Plot";

	    NoFSITrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,20,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0],TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][ TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1 ]);

	    NoFSITrueDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaPTTwoDInDeltaAlphaTLabel,LabelXAxisDeltaPT,TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT].size()-1,&TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices[WhichDeltaAlphaT][0]);

	  }

	  for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {

	    TString DeltaPnTwoDInDeltaAlpha3DqLabel = "DeltaPn_DeltaAlpha3Dq_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3Dq[WhichDeltaAlpha3Dq+1])+"Plot";

	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaPnTwoDInDeltaAlpha3DqLabel,LabelXAxisDeltaPn,20,TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][0],TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][ TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq].size()-1 ]);

	    NoFSITrueDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaPnTwoDInDeltaAlpha3DqLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices[WhichDeltaAlpha3Dq][0]);

	  }

	  for (int WhichDeltaAlpha3DMu = 0; WhichDeltaAlpha3DMu < TwoDNBinsDeltaAlpha3DMu; WhichDeltaAlpha3DMu++) {

	    TString DeltaPnTwoDInDeltaAlpha3DMuLabel = "DeltaPn_DeltaAlpha3DMu_"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaAlpha3DMu[WhichDeltaAlpha3DMu+1])+"Plot";

	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaPnTwoDInDeltaAlpha3DMuLabel,LabelXAxisDeltaPn,20,TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][0],TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][ TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu].size()-1 ]);

	    NoFSITrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaPnTwoDInDeltaAlpha3DMuLabel,LabelXAxisDeltaPn,TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu].size()-1,&TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices[WhichDeltaAlpha3DMu][0]);

	  }

	  for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

	    TString DeltaAlphaTTwoDInDeltaPTLabel = "DeltaAlphaT_DeltaPT_"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPT[WhichDeltaPT+1])+"Plot";

	    NoFSITrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,20,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0],TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][ TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1 ]);

	    NoFSITrueDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaAlphaTTwoDInDeltaPTLabel,LabelXAxisDeltaAlphaT,TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT].size()-1,&TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices[WhichDeltaPT][0]);

	  }

	  for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

	    TString DeltaAlpha3DqTwoDInDeltaPnLabel = "DeltaAlpha3Dq_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";

	    NoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaAlpha3DqTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3Dq,20,TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][0],TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][ TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn].size()-1 ]);

	    NoFSITrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaAlpha3DqTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3Dq,TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices[WhichDeltaPn][0]);

	    TString DeltaAlpha3DMuTwoDInDeltaPnLabel = "DeltaAlpha3DMu_DeltaPn_"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn])+"To"+tools.ConvertToString(TwoDArrayNBinsDeltaPn[WhichDeltaPn+1])+"Plot";

	    NoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"NoFSITrueFineBin"+DeltaAlpha3DMuTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3DMu,20,TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][0],TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][ TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn].size()-1 ]);

	    NoFSITrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn] = new TH1D(InteractionLabels[inte]+"NoFSITrue"+DeltaAlpha3DMuTwoDInDeltaPnLabel,LabelXAxisDeltaAlpha3DMu,TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn].size()-1,&TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices[WhichDeltaPn][0]);

	  }

	  // 2D Fine Binning Correlated

	  SerialNoFSITrueFineBinDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
	  SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaPn_DeltaAlpha3DqPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices)[0]);
	  SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaPn_DeltaAlpha3DMuPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices)[0]);
	  SerialNoFSITrueFineBinDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
	  SerialNoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaAlpha3Dq_DeltaPnPlot",LabelXAxisDeltaAlpha3Dq,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices)[0]);
	  SerialNoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueFineBinDeltaAlpha3DMu_DeltaPnPlot",LabelXAxisDeltaAlpha3DMu,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices)[0]);

	  // 2D Nominal Binning Correlated

	  SerialNoFSITrueDeltaPT_InDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaPT_DeltaAlphaTPlot",LabelXAxisDeltaPT,tools.Return2DNBins(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices)[0]);
	  SerialNoFSITrueDeltaPn_InDeltaAlpha3DqPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaPn_DeltaAlpha3DqPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices)[0]);
	  SerialNoFSITrueDeltaPn_InDeltaAlpha3DMuPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaPn_DeltaAlpha3DMuPlot",LabelXAxisDeltaPn,tools.Return2DNBins(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices)[0]);
	  SerialNoFSITrueDeltaAlphaT_InDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaAlphaT_DeltaPTPlot",LabelXAxisDeltaAlphaT,tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices)[0]);
	  SerialNoFSITrueDeltaAlpha3Dq_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaAlpha3Dq_DeltaPnPlot",LabelXAxisDeltaAlpha3Dq,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices)[0]);
	  SerialNoFSITrueDeltaAlpha3DMu_InDeltaPnPlot[inte] = new TH1D(InteractionLabels[inte]+"SerialNoFSITrueDeltaAlpha3DMu_DeltaPnPlot",LabelXAxisDeltaAlpha3DMu,tools.Return2DNBins(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices),&tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices)[0]);
	  
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
	  if (fOutputFile == "ACHILLES") { weight = fScaleFactor*Units*Weight; }
	  //if (jentry%1000 == 0) { std::cout << "fScaleFactor = " << fScaleFactor << ", weight = " << weight << std::endl; }

	  //	  cout << "fScaleFactor = " << fScaleFactor << ", Weight = " << Weight << endl;

	  //double weight = 1.;	
	  if (fOutputFile == "GiBUU_2021") { weight = weight/100.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2021_Inclusive") { weight = weight/10.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2023_ME_DOW") { weight = weight/5.; } // To increase the stats, the GiBUU sample has been produced in 100 sample
	  if (fOutputFile == "GiBUU_2023") { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 50 samples	
	  if (fOutputFile == "GiBUU_2023_Patch1") { weight = weight/500.; } // To increase the stats, the GiBUU sample has been produced in 500 samples	
	  if (fOutputFile == "GiBUU_2021_NoFSI") { weight = weight/100.; } // To increase the stats, Ben Bogart produced the no fsi GiBUU sample in 100 samples
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
	    double DeltaPL = reco_stv_tool.ReturnPL();
	    double DeltaPtx = reco_stv_tool.ReturnPtx();
	    double DeltaPty = reco_stv_tool.ReturnPty();
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double DeltaAlpha3Dq = reco_stv_tool.ReturnDeltaAlpha3Dq();
	    double DeltaAlpha3DMu = reco_stv_tool.ReturnDeltaAlpha3DMu();
	    double DeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
	    double DeltaPhi3D = reco_stv_tool.ReturnDeltaPhi3D();
	    double DeltaPn = reco_stv_tool.ReturnPn();
	    double DeltaPnPerp = reco_stv_tool.ReturnPnPerp();
	    double DeltaPnPerpx = reco_stv_tool.ReturnPnPerpx();
	    double DeltaPnPerpy = reco_stv_tool.ReturnPnPerpy();
	    double DeltaPnPar = reco_stv_tool.ReturnPnPar();
	    double ECal = reco_stv_tool.ReturnECalMB();
	    double qtrue = Enu_true - Muon4Vector.E();
	    double qreco = ECal - Muon4Vector.E();


	    double DefaultMuonCosTheta = Muon4Vector.CosTheta();
	    double DefaultDeltaPT = reco_stv_tool.ReturnPt();
	    double DefaultDeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double DefaultDeltaAlpha3Dq = reco_stv_tool.ReturnDeltaAlpha3Dq();
	    double DefaultDeltaAlpha3DMu = reco_stv_tool.ReturnDeltaAlpha3DMu();
	    double DefaultDeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
	    double DefaultDeltaPhi3D = reco_stv_tool.ReturnDeltaPhi3D();
	    double DefaultDeltaPn = reco_stv_tool.ReturnPn();
	    double DefaultDeltaPnPerp = reco_stv_tool.ReturnPnPerp();
	    double DefaultDeltaPnPerpx = reco_stv_tool.ReturnPnPerpx();
	    double DefaultDeltaPnPerpy = reco_stv_tool.ReturnPnPerpy();
	    double DefaultDeltaPtx = reco_stv_tool.ReturnPtx();
	    double DefaultDeltaPty = reco_stv_tool.ReturnPty();
	    double DefaultDeltaPnPar = reco_stv_tool.ReturnPnPar();

	    //----------------------------------------//

	    // Dumping values in txt file 

	    OutputTxt << "Event = " << ientry << endl;
	    OutputTxt << "DeltaAlpha3Dq = " << DeltaAlpha3Dq << ", DeltaPn = " << DeltaPn << ", DeltaPL = " << DeltaPL << endl;
	    OutputTxt << "DeltaAlphaT = " << DeltaAlphaT << ", DeltaPT = " << DeltaPT << ", ECal = " << ECal << endl;
	    OutputTxt << "Muon P = " << MuonMomentum << ", Px = " << Muon4Vector.X() << ", Py = " << Muon4Vector.Y() << ", Pz = " << Muon4Vector.Z() << endl;
	    OutputTxt << "Proton P = " << ProtonMomentum << ", Px = " << Proton4Vector.X() << ", Py = " << Proton4Vector.Y() << ", Pz = " << Proton4Vector.Z() << endl << endl;

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    TrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueFineBinDeltaPTPlot[0]->Fill(DeltaPT,weight);
	    TrueFineBinDeltaPnPlot[0]->Fill(DeltaPn,weight);
	    TrueFineBinDeltaAlpha3DMuPlot[0]->Fill(DeltaAlpha3DMu,weight);
	    TrueFineBinDeltaAlpha3DqPlot[0]->Fill(DeltaAlpha3Dq,weight);
	    TrueFineBinDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
	    TrueFineBinDeltaPhiTPlot[0]->Fill(DeltaPhiT,weight);
	    TrueFineBinDeltaPhi3DPlot[0]->Fill(DeltaPhi3D,weight);
	    TrueFineBinDeltaPnPerpPlot[0]->Fill(DeltaPnPerp,weight);
	    TrueFineBinDeltaPnPerpxPlot[0]->Fill(DeltaPnPerpx,weight);
	    TrueFineBinDeltaPnPerpyPlot[0]->Fill(DeltaPnPerpy,weight);
	    TrueFineBinDeltaPtxPlot[0]->Fill(DeltaPtx,weight);
	    TrueFineBinDeltaPtyPlot[0]->Fill(DeltaPty,weight);
	    TrueFineBinDeltaPnParPlot[0]->Fill(DeltaPnPar,weight);

	    TrueFineBinECal2DPlot[0]->Fill(Enu_true,ECal,weight);
	    TrueFineBinq2DPlot[0]->Fill(qtrue,qreco,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    TrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueFineBinDeltaPTPlot[genie_mode]->Fill(DeltaPT,weight);
	    TrueFineBinDeltaPnPlot[genie_mode]->Fill(DeltaPn,weight);
	    TrueFineBinDeltaAlpha3DMuPlot[genie_mode]->Fill(DeltaAlpha3DMu,weight);
	    TrueFineBinDeltaAlpha3DqPlot[genie_mode]->Fill(DeltaAlpha3Dq,weight);
	    TrueFineBinDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
	    TrueFineBinDeltaPhiTPlot[genie_mode]->Fill(DeltaPhiT,weight);
	    TrueFineBinDeltaPhi3DPlot[genie_mode]->Fill(DeltaPhi3D,weight);
	    TrueFineBinDeltaPnPerpPlot[genie_mode]->Fill(DeltaPnPerp,weight);
	    TrueFineBinDeltaPnPerpxPlot[genie_mode]->Fill(DeltaPnPerpx,weight);
	    TrueFineBinDeltaPnPerpyPlot[genie_mode]->Fill(DeltaPnPerpy,weight);
	    TrueFineBinDeltaPtxPlot[genie_mode]->Fill(DeltaPtx,weight);
	    TrueFineBinDeltaPtyPlot[genie_mode]->Fill(DeltaPty,weight);
	    TrueFineBinDeltaPnParPlot[genie_mode]->Fill(DeltaPnPar,weight);

	    TrueFineBinECal2DPlot[genie_mode]->Fill(Enu_true,ECal,weight);
            TrueFineBinq2DPlot[genie_mode]->Fill(qtrue,qreco,weight);

	    //----------------------------------------//

	    // Overflow & underflow bins

	    // Overflow
	    // Affects

	    // DeltaPT                                                                                                                                                                        
	    // DeltaPn                                                                                                                                                             
	    // DeltaPnPerp                                                                                                                                                                    
	    // DeltaPnPerpx                                                                                                                                                                    
	    // DeltaPnPerpy                                                                                                                                                                    
	    // DeltaPtx                                                                                                                                                                    
	    // DeltaPty                                                                                                                                                                    
	    // DeltaPnPar                                                                                                                                                                     
	    if (DeltaPT > ArrayNBinsDeltaPT[NBinsDeltaPT]) { DeltaPT = 0.9999 * ArrayNBinsDeltaPT[NBinsDeltaPT]; }
	    if (DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { DeltaPn = 0.9999 * ArrayNBinsDeltaPn[NBinsDeltaPn]; }
	    if (DeltaPnPerp > ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]) { DeltaPnPerp = 0.9999 * ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]; }
	    if (DeltaPnPerpx > ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]) { DeltaPnPerpx = 0.9999 * ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]; }
	    if (DeltaPnPerpy > ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]) { DeltaPnPerpy = 0.9999 * ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]; }
	    if (DeltaPtx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { DeltaPtx = 0.9999 * ArrayNBinsDeltaPtx[NBinsDeltaPtx]; }
	    if (DeltaPty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { DeltaPty = 0.9999 * ArrayNBinsDeltaPty[NBinsDeltaPty]; }
	    if (DeltaPnPar > ArrayNBinsDeltaPnPar[NBinsDeltaPnPar]) { DeltaPnPar = 0.9999 * ArrayNBinsDeltaPnPar[NBinsDeltaPnPar]; } \

	    //----------------------------------------//

	    // Underflow
	    // Affects

	    // DeltaPnPerp                                                                                                                                                                    
	    // DeltaPnPerpx                                                                                                                                                                    
	    // DeltaPnPerpy                                                                                                                                                                    
	    // DeltaPtx                                                                                                                                                                    
	    // DeltaPty                                                                                                                                                                    
	    // DeltaPnPar                                                                                                                                                                     

	    if (DeltaPnPerp < ArrayNBinsDeltaPnPerp[0]) { DeltaPnPerp = 1.0001 * ArrayNBinsDeltaPnPerp[0]; }
	    if (DeltaPnPerpx < ArrayNBinsDeltaPnPerpx[0]) { DeltaPnPerpx = 1.0001 * ArrayNBinsDeltaPnPerpx[0]; }
	    if (DeltaPnPerpy < ArrayNBinsDeltaPnPerpy[0]) { DeltaPnPerpy = 1.0001 * ArrayNBinsDeltaPnPerpy[0]; }
	    if (DeltaPtx < ArrayNBinsDeltaPtx[0]) { DeltaPtx = 1.0001 * ArrayNBinsDeltaPtx[0]; }
	    if (DeltaPty < ArrayNBinsDeltaPty[0]) { DeltaPty = 1.0001 * ArrayNBinsDeltaPty[0]; }
	    if (DeltaPnPar < ArrayNBinsDeltaPnPar[0]) { DeltaPnPar = 1.0001 * ArrayNBinsDeltaPnPar[0]; }	\

	    //----------------------------------------//

	    // Indices for 2D analysis

	    int DeltaPTTwoDIndex = tools.ReturnIndex(DeltaPT, TwoDArrayNBinsDeltaPT); 
	    if (DeltaPTTwoDIndex < 0) { continue; }
	    int DeltaPnTwoDIndex = tools.ReturnIndex(DeltaPn, TwoDArrayNBinsDeltaPn);
	    if (DeltaPnTwoDIndex < 0) { continue; }
	    int DeltaAlphaTTwoDIndex = tools.ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
	    if (DeltaAlphaTTwoDIndex < 0) { continue; }
	    int DeltaAlpha3DqTwoDIndex = tools.ReturnIndex(DeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
	    if (DeltaAlpha3DqTwoDIndex < 0) { continue; }
	    int DeltaAlpha3DMuTwoDIndex = tools.ReturnIndex(DeltaAlpha3DMu, TwoDArrayNBinsDeltaAlpha3DMu);
	    if (DeltaAlpha3DMuTwoDIndex < 0) { continue; }

	    int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,DeltaAlphaT);
	    int SerialDeltaAlpha3DqInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices,DeltaPnTwoDIndex,DeltaAlpha3Dq);
	    int SerialDeltaAlpha3DMuInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices,DeltaPnTwoDIndex,DeltaAlpha3DMu);
	    int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,DeltaPT);
	    int SerialDeltaPnInDeltaAlpha3DqIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices,DeltaAlpha3DqTwoDIndex,DeltaPn);
	    int SerialDeltaPnInDeltaAlpha3DMuIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices,DeltaAlpha3DMuTwoDIndex,DeltaPn);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 2D Fine Binning Uncorrelated

	    TrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(DefaultDeltaAlphaT,weight);
	    TrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3Dq,weight);
	    TrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3DMu,weight);
	    TrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(DefaultDeltaPT,weight);
	    TrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[0][DeltaAlpha3DqTwoDIndex]->Fill(DefaultDeltaPn,weight);
	    TrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[0][DeltaAlpha3DMuTwoDIndex]->Fill(DefaultDeltaPn,weight);

	    // 2D Fine Binning Correlated

	    SerialTrueFineBinDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialTrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialTrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialTrueFineBinDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialTrueFineBinDeltaPn_InDeltaAlpha3DqPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialTrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    // 1D Nominal Binning

	    TrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    TrueDeltaPTPlot[0]->Fill(DeltaPT,weight);
	    TrueDeltaPnPlot[0]->Fill(DeltaPn,weight);
	    TrueDeltaAlpha3DMuPlot[0]->Fill(DeltaAlpha3DMu,weight);
	    TrueDeltaAlpha3DqPlot[0]->Fill(DeltaAlpha3Dq,weight);
	    TrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
	    TrueDeltaPhiTPlot[0]->Fill(DeltaPhiT,weight);
	    TrueDeltaPhi3DPlot[0]->Fill(DeltaPhi3D,weight);
	    TrueDeltaPnPerpPlot[0]->Fill(DeltaPnPerp,weight);
	    TrueDeltaPnPerpxPlot[0]->Fill(DeltaPnPerpx,weight);
	    TrueDeltaPnPerpyPlot[0]->Fill(DeltaPnPerpy,weight);
	    TrueDeltaPtxPlot[0]->Fill(DeltaPtx,weight);
	    TrueDeltaPtyPlot[0]->Fill(DeltaPty,weight);
	    TrueDeltaPnParPlot[0]->Fill(DeltaPnPar,weight);

	    // 2D Nominal Binning Uncorrelated

	    TrueDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);
	    TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DeltaAlpha3Dq,weight);
	    TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DeltaAlpha3DMu,weight);
	    TrueDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(DeltaPT,weight);
	    TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[0][DeltaAlpha3DqTwoDIndex]->Fill(DeltaPn,weight);
	    TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[0][DeltaAlpha3DMuTwoDIndex]->Fill(DeltaPn,weight);

	    // 2D Nominal Binning Correlated

	    SerialTrueDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialTrueDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialTrueDeltaPn_InDeltaAlpha3DqPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    //----------------------------------------//

	    // filling in the histo based on the interaction mode

	    // 2D Fine Binning Uncorrelated

	    TrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(DefaultDeltaAlphaT,weight);
	    TrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3Dq,weight);
	    TrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3DMu,weight);
	    TrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(DefaultDeltaPT,weight);
	    TrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[genie_mode][DeltaAlpha3DqTwoDIndex]->Fill(DefaultDeltaPn,weight);
	    TrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[genie_mode][DeltaAlpha3DMuTwoDIndex]->Fill(DefaultDeltaPn,weight);

	    // 2D Fine Binning Correlated

	    SerialTrueFineBinDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialTrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialTrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialTrueFineBinDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialTrueFineBinDeltaPn_InDeltaAlpha3DqPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialTrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    // 1D Nominal Binning

	    TrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    TrueDeltaPTPlot[genie_mode]->Fill(DeltaPT,weight);
	    TrueDeltaPnPlot[genie_mode]->Fill(DeltaPn,weight);
	    TrueDeltaAlpha3DMuPlot[genie_mode]->Fill(DeltaAlpha3DMu,weight);
	    TrueDeltaAlpha3DqPlot[genie_mode]->Fill(DeltaAlpha3Dq,weight);
	    TrueDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
	    TrueDeltaPhiTPlot[genie_mode]->Fill(DeltaPhiT,weight);
	    TrueDeltaPhi3DPlot[genie_mode]->Fill(DeltaPhi3D,weight);
	    TrueDeltaPnPerpPlot[genie_mode]->Fill(DeltaPnPerp,weight);
	    TrueDeltaPnPerpxPlot[genie_mode]->Fill(DeltaPnPerpx,weight);
	    TrueDeltaPnPerpyPlot[genie_mode]->Fill(DeltaPnPerpy,weight);
	    TrueDeltaPtxPlot[genie_mode]->Fill(DeltaPtx,weight);
	    TrueDeltaPtyPlot[genie_mode]->Fill(DeltaPty,weight);
	    TrueDeltaPnParPlot[genie_mode]->Fill(DeltaPnPar,weight);

	    // 2D Nominal Binning Uncorrelated

	    TrueDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);
	    TrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DeltaAlpha3Dq,weight);
	    TrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DeltaAlpha3DMu,weight);
	    TrueDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(DeltaPT,weight);
	    TrueDeltaPn_InDeltaAlpha3DqTwoDPlot[genie_mode][DeltaAlpha3DqTwoDIndex]->Fill(DeltaPn,weight);
	    TrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[genie_mode][DeltaAlpha3DMuTwoDIndex]->Fill(DeltaPn,weight);

	    // 2D Nominal Binning Correlated

	    SerialTrueDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialTrueDeltaAlpha3Dq_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialTrueDeltaAlpha3DMu_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialTrueDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialTrueDeltaPn_InDeltaAlpha3DqPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialTrueDeltaPn_InDeltaAlpha3DMuPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    //----------------------------------------//

	    // ECal overflow / underflow bins
	    if (DeltaPtx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { DeltaPtx = 0.5 * (ArrayNBinsDeltaPtx[NBinsDeltaPtx] + ArrayNBinsDeltaPtx[NBinsDeltaPtx-1]); }
	    if (DeltaPtx < ArrayNBinsDeltaPtx[0]) { DeltaPtx = 0.5 * (ArrayNBinsDeltaPtx[0] + ArrayNBinsDeltaPtx[1]); }
	    if (DeltaPty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { DeltaPty = 0.5 * (ArrayNBinsDeltaPty[NBinsDeltaPty] + ArrayNBinsDeltaPty[NBinsDeltaPty-1]); }
	    if (DeltaPty < ArrayNBinsDeltaPty[0]) { DeltaPty = 0.5 * (ArrayNBinsDeltaPty[0] + ArrayNBinsDeltaPty[1]); }
	    if (ECal > ArrayNBinsECal[NBinsECal]) { ECal = 0.5 * (ArrayNBinsECal[NBinsECal] + ArrayNBinsECal[NBinsECal-1]); }
	    if (ECal < ArrayNBinsECal[0]) { ECal = 0.5 * (ArrayNBinsECal[0] + ArrayNBinsECal[1]); }

	    // 3D indices

	    int MuonMomentumTwoDIndex = tools.ReturnIndex(MuonMomentum, TwoDArrayNBinsMuonMomentum);
	    if (MuonMomentumTwoDIndex < 0) { cout << "MuonMomentumTwoDIndex = " << MuonMomentumTwoDIndex << endl; }
	    int ProtonMomentumTwoDIndex = tools.ReturnIndex(ProtonMomentum, TwoDArrayNBinsProtonMomentum);
	    if (ProtonMomentumTwoDIndex < 0) { cout << "ProtonMomentumTwoDIndex = " << ProtonMomentumTwoDIndex << endl; }
	    int MuonCosThetaTwoDIndex = tools.ReturnIndex(MuonCosTheta, TwoDArrayNBinsMuonCosTheta);
	    if (MuonMomentumTwoDIndex < 0) { cout << "MuonMomentumTwoDIndex = " << MuonMomentumTwoDIndex << endl; }
	    int ProtonCosThetaTwoDIndex = tools.ReturnIndex(ProtonCosTheta, TwoDArrayNBinsProtonCosTheta);
	    if (ProtonMomentumTwoDIndex < 0) { cout << "ProtonMomentumTwoDIndex = " << ProtonMomentumTwoDIndex << endl; }
	    int DeltaPtxTwoDIndex = tools.ReturnIndex(DeltaPtx, TwoDArrayNBinsDeltaPtx);
	    if (DeltaPtxTwoDIndex < 0) { cout << "DeltaPtxTwoDIndex = " << DeltaPtxTwoDIndex << endl; }
	    int DeltaPtyTwoDIndex = tools.ReturnIndex(DeltaPty, TwoDArrayNBinsDeltaPty);
	    if (DeltaPtyTwoDIndex < 0) { cout << "DeltaPtyTwoDIndex = " << DeltaPtyTwoDIndex << endl; }

	    int SerialECalInDeltaPTDeltaAlphaTIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPTDeltaAlphaTSlices,DeltaPTTwoDIndex,DeltaAlphaTTwoDIndex,ECal);
	    int SerialECalInDeltaPtxDeltaPtyIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInDeltaPtxDeltaPtySlices,DeltaPtxTwoDIndex,DeltaPtyTwoDIndex,ECal);
	    int SerialECalInMuonCosThetaMuonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInMuonCosThetaMuonMomentumSlices,MuonCosThetaTwoDIndex,MuonMomentumTwoDIndex,ECal);
	    int SerialECalInProtonCosThetaProtonMomentumIndex = tools.ReturnIndexIn3DList(TwoDArrayNBinsECalInProtonCosThetaProtonMomentumSlices,ProtonCosThetaTwoDIndex,ProtonMomentumTwoDIndex,ECal);

	    // 3D filling for all events

	    TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[0][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InDeltaPtxDeltaPtyTwoDPlot[0][DeltaPtxTwoDIndex][DeltaPtyTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[0][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[0][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);

	    SerialTrueECal_InDeltaPTDeltaAlphaTPlot[0]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
	    SerialTrueECal_InDeltaPtxDeltaPtyPlot[0]->Fill(SerialECalInDeltaPtxDeltaPtyIndex,weight);
	    SerialTrueECal_InMuonCosThetaMuonMomentumPlot[0]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
	    SerialTrueECal_InProtonCosThetaProtonMomentumPlot[0]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

	    // 3D filling based on interaction mode

	    TrueECal_InDeltaPTDeltaAlphaTTwoDPlot[genie_mode][DeltaPTTwoDIndex][DeltaAlphaTTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InDeltaPtxDeltaPtyTwoDPlot[genie_mode][DeltaPtxTwoDIndex][DeltaPtyTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InMuonCosThetaMuonMomentumTwoDPlot[genie_mode][MuonCosThetaTwoDIndex][MuonMomentumTwoDIndex]->Fill(ECal,weight);
	    TrueECal_InProtonCosThetaProtonMomentumTwoDPlot[genie_mode][ProtonCosThetaTwoDIndex][ProtonMomentumTwoDIndex]->Fill(ECal,weight);

	    SerialTrueECal_InDeltaPTDeltaAlphaTPlot[genie_mode]->Fill(SerialECalInDeltaPTDeltaAlphaTIndex,weight);
	    SerialTrueECal_InDeltaPtxDeltaPtyPlot[genie_mode]->Fill(SerialECalInDeltaPtxDeltaPtyIndex,weight);
	    SerialTrueECal_InMuonCosThetaMuonMomentumPlot[genie_mode]->Fill(SerialECalInMuonCosThetaMuonMomentumIndex,weight);
	    SerialTrueECal_InProtonCosThetaProtonMomentumPlot[genie_mode]->Fill(SerialECalInProtonCosThetaProtonMomentumIndex,weight);

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
	    double DeltaPT = reco_stv_tool.ReturnPt();
	    double DeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double DeltaAlpha3Dq = reco_stv_tool.ReturnDeltaAlpha3Dq();
	    double DeltaAlpha3DMu = reco_stv_tool.ReturnDeltaAlpha3DMu();
	    double DeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
	    double DeltaPhi3D = reco_stv_tool.ReturnDeltaPhi3D();
	    double DeltaPn = reco_stv_tool.ReturnPn();
	    double DeltaPnPerp = reco_stv_tool.ReturnPnPerp();
	    double DeltaPnPerpx = reco_stv_tool.ReturnPnPerpx();
	    double DeltaPnPerpy = reco_stv_tool.ReturnPnPerpy();
	    double DeltaPtx = reco_stv_tool.ReturnPtx();
	    double DeltaPty = reco_stv_tool.ReturnPty();
	    double DeltaPnPar = reco_stv_tool.ReturnPnPar();
            double ECal = reco_stv_tool.ReturnECalMB();
            double qtrue = Enu_true - Muon4Vector.E();
            double qreco = ECal - Muon4Vector.E();

	    double DefaultMuonCosTheta = Muon4Vector.CosTheta();
	    double DefaultDeltaPT = reco_stv_tool.ReturnPt();
	    double DefaultDeltaAlphaT = reco_stv_tool.ReturnDeltaAlphaT();
	    double DefaultDeltaAlpha3Dq = reco_stv_tool.ReturnDeltaAlpha3Dq();
	    double DefaultDeltaAlpha3DMu = reco_stv_tool.ReturnDeltaAlpha3DMu();
	    double DefaultDeltaPhiT = reco_stv_tool.ReturnDeltaPhiT();
	    double DefaultDeltaPhi3D = reco_stv_tool.ReturnDeltaPhi3D();
	    double DefaultDeltaPn = reco_stv_tool.ReturnPn();
	    double DefaultDeltaPnPerp = reco_stv_tool.ReturnPnPerp();
	    double DefaultDeltaPnPerpx = reco_stv_tool.ReturnPnPerpx();
	    double DefaultDeltaPnPerpy = reco_stv_tool.ReturnPnPerpy();
	    double DefaultDeltaPtx = reco_stv_tool.ReturnPtx();
	    double DefaultDeltaPty = reco_stv_tool.ReturnPty();
	    double DefaultDeltaPnPar = reco_stv_tool.ReturnPnPar();

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinDeltaPTPlot[0]->Fill(DeltaPT,weight);
	    NoFSITrueFineBinDeltaPnPlot[0]->Fill(DeltaPn,weight);
	    NoFSITrueFineBinDeltaAlpha3DMuPlot[0]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueFineBinDeltaAlpha3DqPlot[0]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueFineBinDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
	    NoFSITrueFineBinDeltaPhiTPlot[0]->Fill(DeltaPhiT,weight);
	    NoFSITrueFineBinDeltaPhi3DPlot[0]->Fill(DeltaPhi3D,weight);
	    NoFSITrueFineBinDeltaPnPerpPlot[0]->Fill(DeltaPnPerp,weight);
	    NoFSITrueFineBinDeltaPnPerpxPlot[0]->Fill(DeltaPnPerpx,weight);
	    NoFSITrueFineBinDeltaPnPerpyPlot[0]->Fill(DeltaPnPerpy,weight);
	    NoFSITrueFineBinDeltaPtxPlot[0]->Fill(DeltaPtx,weight);
	    NoFSITrueFineBinDeltaPtyPlot[0]->Fill(DeltaPty,weight);
	    NoFSITrueFineBinDeltaPnParPlot[0]->Fill(DeltaPnPar,weight);

	    NoFSITrueFineBinECal2DPlot[0]->Fill(Enu_true,ECal,weight);
	    NoFSITrueFineBinq2DPlot[0]->Fill(qtrue,qreco,weight);

	    // filling in the histo based on the interaction mode

	    // 1D Fine Binning

	    NoFSITrueFineBinMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    NoFSITrueFineBinDeltaPTPlot[genie_mode]->Fill(DeltaPT,weight);
	    NoFSITrueFineBinDeltaPnPlot[genie_mode]->Fill(DeltaPn,weight);
	    NoFSITrueFineBinDeltaAlpha3DMuPlot[genie_mode]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueFineBinDeltaAlpha3DqPlot[genie_mode]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueFineBinDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
	    NoFSITrueFineBinDeltaPhiTPlot[genie_mode]->Fill(DeltaPhiT,weight);
	    NoFSITrueFineBinDeltaPhi3DPlot[genie_mode]->Fill(DeltaPhi3D,weight);
	    NoFSITrueFineBinDeltaPnPerpPlot[genie_mode]->Fill(DeltaPnPerp,weight);
	    NoFSITrueFineBinDeltaPnPerpxPlot[genie_mode]->Fill(DeltaPnPerpx,weight);
	    NoFSITrueFineBinDeltaPnPerpyPlot[genie_mode]->Fill(DeltaPnPerpy,weight);
	    NoFSITrueFineBinDeltaPtxPlot[genie_mode]->Fill(DeltaPtx,weight);
	    NoFSITrueFineBinDeltaPtyPlot[genie_mode]->Fill(DeltaPty,weight);
	    NoFSITrueFineBinDeltaPnParPlot[genie_mode]->Fill(DeltaPnPar,weight);

	    NoFSITrueFineBinECal2DPlot[genie_mode]->Fill(Enu_true,ECal,weight);
	    NoFSITrueFineBinq2DPlot[genie_mode]->Fill(qtrue,qreco,weight);

	    //----------------------------------------//

	    // Overflow & underflow bins

	    // Overflow
	    // Affects

	    // DeltaPT                                                                                                                                                                        
	    // DeltaPn                                                                                                                                                                        
	    // DeltaPnPerp                                                                                                                                                                    
	    // DeltaPnPerpx                                                                                                                                                                    
	    // DeltaPnPerpy                                                                                                                                                                    
	    // DeltaPtx                                                                                                                                                                    
	    // DeltaPty                                                                                                                                                                    
	    // DeltaPnPar                                                                                                                                                                     

	    if (DeltaPT > ArrayNBinsDeltaPT[NBinsDeltaPT]) { DeltaPT = 0.9999 * ArrayNBinsDeltaPT[NBinsDeltaPT]; }
	    if (DeltaPn > ArrayNBinsDeltaPn[NBinsDeltaPn]) { DeltaPn = 0.9999 * ArrayNBinsDeltaPn[NBinsDeltaPn]; }
	    if (DeltaPnPerp > ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]) { DeltaPnPerp = 0.9999 * ArrayNBinsDeltaPnPerp[NBinsDeltaPnPerp]; }
	    if (DeltaPnPerpx > ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]) { DeltaPnPerpx = 0.9999 * ArrayNBinsDeltaPnPerpx[NBinsDeltaPnPerpx]; }
	    if (DeltaPnPerpy > ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]) { DeltaPnPerpy = 0.9999 * ArrayNBinsDeltaPnPerpy[NBinsDeltaPnPerpy]; }
	    if (DeltaPtx > ArrayNBinsDeltaPtx[NBinsDeltaPtx]) { DeltaPtx = 0.9999 * ArrayNBinsDeltaPtx[NBinsDeltaPtx]; }
	    if (DeltaPty > ArrayNBinsDeltaPty[NBinsDeltaPty]) { DeltaPty = 0.9999 * ArrayNBinsDeltaPty[NBinsDeltaPty]; }
	    if (DeltaPnPar > ArrayNBinsDeltaPnPar[NBinsDeltaPnPar]) { DeltaPnPar = 0.9999 * ArrayNBinsDeltaPnPar[NBinsDeltaPnPar]; } \

	    //----------------------------------------//

	    // Underflow
	    // Affects

	    // DeltaPnPerp                                                                                                                                                                    
	    // DeltaPnPerpx                                                                                                                                                                    
	    // DeltaPnPerpy                                                                                                                                                                    
	    // DeltaPtx                                                                                                                                                                    
	    // DeltaPty                                                                                                                                                                    
	    // DeltaPnPar                                                                                                                                                                     

	    if (DeltaPnPerp < ArrayNBinsDeltaPnPerp[0]) { DeltaPnPerp = 1.0001 * ArrayNBinsDeltaPnPerp[0]; }
	    if (DeltaPnPerpx < ArrayNBinsDeltaPnPerpx[0]) { DeltaPnPerpx = 1.0001 * ArrayNBinsDeltaPnPerpx[0]; }
	    if (DeltaPnPerpy < ArrayNBinsDeltaPnPerpy[0]) { DeltaPnPerpx = 1.0001 * ArrayNBinsDeltaPnPerpy[0]; }
	    if (DeltaPtx < ArrayNBinsDeltaPtx[0]) { DeltaPtx = 1.0001 * ArrayNBinsDeltaPtx[0]; }
	    if (DeltaPty < ArrayNBinsDeltaPty[0]) { DeltaPtx = 1.0001 * ArrayNBinsDeltaPty[0]; }
	    if (DeltaPnPar < ArrayNBinsDeltaPnPar[0]) { DeltaPnPar = 1.0001 * ArrayNBinsDeltaPnPar[0]; }	

	    //----------------------------------------//

	    // Indices for 2D analysis

	    int DeltaPTTwoDIndex = tools.ReturnIndex(DeltaPT, TwoDArrayNBinsDeltaPT); 
	    if (DeltaPTTwoDIndex < 0) { continue; }
	    int DeltaPnTwoDIndex = tools.ReturnIndex(DeltaPn, TwoDArrayNBinsDeltaPn);
	    if (DeltaPnTwoDIndex < 0) { continue; }
	    int DeltaAlphaTTwoDIndex = tools.ReturnIndex(DeltaAlphaT, TwoDArrayNBinsDeltaAlphaT);
	    if (DeltaAlphaTTwoDIndex < 0) { continue; }
	    int DeltaAlpha3DqTwoDIndex = tools.ReturnIndex(DeltaAlpha3Dq, TwoDArrayNBinsDeltaAlpha3Dq);
	    if (DeltaAlpha3DqTwoDIndex < 0) { continue; }
	    int DeltaAlpha3DMuTwoDIndex = tools.ReturnIndex(DeltaAlpha3DMu, TwoDArrayNBinsDeltaAlpha3DMu);
	    if (DeltaAlpha3DMuTwoDIndex < 0) { continue; }

	    int SerialDeltaAlphaTInDeltaPTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlphaTInDeltaPTSlices,DeltaPTTwoDIndex,DeltaAlphaT);
	    int SerialDeltaAlpha3DqInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DqInDeltaPnSlices,DeltaPnTwoDIndex,DeltaAlpha3Dq);
	    int SerialDeltaAlpha3DMuInDeltaPnIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaAlpha3DMuInDeltaPnSlices,DeltaPnTwoDIndex,DeltaAlpha3DMu);
	    int SerialDeltaPTInDeltaAlphaTIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPTInDeltaAlphaTSlices,DeltaAlphaTTwoDIndex,DeltaPT);
	    int SerialDeltaPnInDeltaAlpha3DqIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DqSlices,DeltaAlpha3DqTwoDIndex,DeltaPn);
	    int SerialDeltaPnInDeltaAlpha3DMuIndex = tools.ReturnIndexIn2DList(TwoDArrayNBinsDeltaPnInDeltaAlpha3DMuSlices,DeltaAlpha3DMuTwoDIndex,DeltaPn);

	    //----------------------------------------//

	    // filling in the histo regardless of interaction mode

	    // 2D Fine Binning Uncorrelated

	    NoFSITrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(DefaultDeltaAlphaT,weight);
	    NoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3Dq,weight);
	    NoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3DMu,weight);
	    NoFSITrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(DefaultDeltaPT,weight);
	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[0][DeltaAlpha3DqTwoDIndex]->Fill(DefaultDeltaPn,weight);
	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[0][DeltaAlpha3DMuTwoDIndex]->Fill(DefaultDeltaPn,weight);

	    // 2D Fine Binning Correlated

	    SerialNoFSITrueFineBinDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialNoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialNoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialNoFSITrueFineBinDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DqPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    // 1D Nominal Binning

	    NoFSITrueMuonCosThetaPlot[0]->Fill(MuonCosTheta,weight);
	    NoFSITrueDeltaPTPlot[0]->Fill(DeltaPT,weight);
	    NoFSITrueDeltaPnPlot[0]->Fill(DeltaPn,weight);
	    NoFSITrueDeltaAlpha3DMuPlot[0]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueDeltaAlpha3DqPlot[0]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
	    NoFSITrueDeltaPhiTPlot[0]->Fill(DeltaPhiT,weight);
	    NoFSITrueDeltaPhi3DPlot[0]->Fill(DeltaPhi3D,weight);
	    NoFSITrueDeltaPnPerpPlot[0]->Fill(DeltaPnPerp,weight);
	    NoFSITrueDeltaPnPerpxPlot[0]->Fill(DeltaPnPerpx,weight);
	    NoFSITrueDeltaPnPerpyPlot[0]->Fill(DeltaPnPerpy,weight);
	    NoFSITrueDeltaPtxPlot[0]->Fill(DeltaPtx,weight);
	    NoFSITrueDeltaPtyPlot[0]->Fill(DeltaPty,weight);
	    NoFSITrueDeltaPnParPlot[0]->Fill(DeltaPnPar,weight);

	    // 2D Nominal Binning Uncorrelated

	    NoFSITrueDeltaAlphaT_InDeltaPTTwoDPlot[0][DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);
	    NoFSITrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[0][DeltaPnTwoDIndex]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueDeltaPT_InDeltaAlphaTTwoDPlot[0][DeltaAlphaTTwoDIndex]->Fill(DeltaPT,weight);
	    NoFSITrueDeltaPn_InDeltaAlpha3DqTwoDPlot[0][DeltaAlpha3DqTwoDIndex]->Fill(DeltaPn,weight);
	    NoFSITrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[0][DeltaAlpha3DMuTwoDIndex]->Fill(DeltaPn,weight);

	    // 2D Nominal Binning Correlated

	    SerialNoFSITrueDeltaAlphaT_InDeltaPTPlot[0]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialNoFSITrueDeltaAlpha3Dq_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialNoFSITrueDeltaAlpha3DMu_InDeltaPnPlot[0]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialNoFSITrueDeltaPT_InDeltaAlphaTPlot[0]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialNoFSITrueDeltaPn_InDeltaAlpha3DqPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialNoFSITrueDeltaPn_InDeltaAlpha3DMuPlot[0]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    //----------------------------------------//

	    // filling in the histo based on the interaction mode

	    // 2D Fine Binning Uncorrelated

	    NoFSITrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(DefaultDeltaAlphaT,weight);
	    NoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3Dq,weight);
	    NoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DefaultDeltaAlpha3DMu,weight);
	    NoFSITrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(DefaultDeltaPT,weight);
	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[genie_mode][DeltaAlpha3DqTwoDIndex]->Fill(DefaultDeltaPn,weight);
	    NoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[genie_mode][DeltaAlpha3DMuTwoDIndex]->Fill(DefaultDeltaPn,weight);

	    // 2D Fine Binning Correlated

	    SerialNoFSITrueFineBinDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialNoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialNoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialNoFSITrueFineBinDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DqPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialNoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);

	    // 1D Nominal Binning

	    NoFSITrueMuonCosThetaPlot[genie_mode]->Fill(MuonCosTheta,weight);
	    NoFSITrueDeltaPTPlot[genie_mode]->Fill(DeltaPT,weight);
	    NoFSITrueDeltaPnPlot[genie_mode]->Fill(DeltaPn,weight);
	    NoFSITrueDeltaAlpha3DMuPlot[genie_mode]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueDeltaAlpha3DqPlot[genie_mode]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
	    NoFSITrueDeltaPhiTPlot[genie_mode]->Fill(DeltaPhiT,weight);
	    NoFSITrueDeltaPhi3DPlot[genie_mode]->Fill(DeltaPhi3D,weight);
	    NoFSITrueDeltaPnPerpPlot[genie_mode]->Fill(DeltaPnPerp,weight);
	    NoFSITrueDeltaPnPerpxPlot[genie_mode]->Fill(DeltaPnPerpx,weight);
	    NoFSITrueDeltaPnPerpyPlot[genie_mode]->Fill(DeltaPnPerpy,weight);
	    NoFSITrueDeltaPtxPlot[genie_mode]->Fill(DeltaPtx,weight);
	    NoFSITrueDeltaPtyPlot[genie_mode]->Fill(DeltaPty,weight);
	    NoFSITrueDeltaPnParPlot[genie_mode]->Fill(DeltaPnPar,weight);

	    // 2D Nominal Binning Uncorrelated

	    NoFSITrueDeltaAlphaT_InDeltaPTTwoDPlot[genie_mode][DeltaPTTwoDIndex]->Fill(DeltaAlphaT,weight);
	    NoFSITrueDeltaAlpha3Dq_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DeltaAlpha3Dq,weight);
	    NoFSITrueDeltaAlpha3DMu_InDeltaPnTwoDPlot[genie_mode][DeltaPnTwoDIndex]->Fill(DeltaAlpha3DMu,weight);
	    NoFSITrueDeltaPT_InDeltaAlphaTTwoDPlot[genie_mode][DeltaAlphaTTwoDIndex]->Fill(DeltaPT,weight);
	    NoFSITrueDeltaPn_InDeltaAlpha3DqTwoDPlot[genie_mode][DeltaAlpha3DqTwoDIndex]->Fill(DeltaPn,weight);
	    NoFSITrueDeltaPn_InDeltaAlpha3DMuTwoDPlot[genie_mode][DeltaAlpha3DMuTwoDIndex]->Fill(DeltaPn,weight);

	    // 2D Nominal Binning Correlated

	    SerialNoFSITrueDeltaAlphaT_InDeltaPTPlot[genie_mode]->Fill(SerialDeltaAlphaTInDeltaPTIndex,weight);
	    SerialNoFSITrueDeltaAlpha3Dq_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DqInDeltaPnIndex,weight);
	    SerialNoFSITrueDeltaAlpha3DMu_InDeltaPnPlot[genie_mode]->Fill(SerialDeltaAlpha3DMuInDeltaPnIndex,weight);
	    SerialNoFSITrueDeltaPT_InDeltaAlphaTPlot[genie_mode]->Fill(SerialDeltaPTInDeltaAlphaTIndex,weight);
	    SerialNoFSITrueDeltaPn_InDeltaAlpha3DqPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DqIndex,weight);
	    SerialNoFSITrueDeltaPn_InDeltaAlpha3DMuPlot[genie_mode]->Fill(SerialDeltaPnInDeltaAlpha3DMuIndex,weight);


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
		Reweight(TrueFineBinDeltaPTPlot[inte]);
		Reweight(TrueFineBinDeltaPnPlot[inte]);
		Reweight(TrueFineBinDeltaPnPerpPlot[inte]);
		Reweight(TrueFineBinDeltaPnPerpxPlot[inte]);
		Reweight(TrueFineBinDeltaPnPerpyPlot[inte]);
		Reweight(TrueFineBinDeltaPtxPlot[inte]);
		Reweight(TrueFineBinDeltaPtyPlot[inte]);
		Reweight(TrueFineBinDeltaPnParPlot[inte]);
		Reweight(TrueFineBinDeltaAlphaTPlot[inte]);
		Reweight(TrueFineBinDeltaAlpha3DqPlot[inte]);
		Reweight(TrueFineBinDeltaAlpha3DMuPlot[inte]);
		Reweight(TrueFineBinDeltaPhiTPlot[inte]);
		Reweight(TrueFineBinDeltaPhi3DPlot[inte]);

	        // Pre FSI

	        // 1D Fine Binning

		Reweight(NoFSITrueFineBinMuonCosThetaPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPTPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPnPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPnPerpPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPnPerpxPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPnPerpyPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPtxPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPtyPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPnParPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaAlphaTPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaAlpha3DqPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaAlpha3DMuPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPhiTPlot[inte]);
		Reweight(NoFSITrueFineBinDeltaPhi3DPlot[inte]);

		//----------------------------------------//

		// 2D Fine Binning Uncorrelated

		for (int WhichDeltaPn = 0; WhichDeltaPn < TwoDNBinsDeltaPn; WhichDeltaPn++) {

		  Reweight(TrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPn+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPn) ) );
		  Reweight(TrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPn+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPn) ) );

		  Reweight(NoFSITrueFineBinDeltaAlpha3Dq_InDeltaPnTwoDPlot[inte][WhichDeltaPn], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPn+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPn) ) );
		  Reweight(NoFSITrueFineBinDeltaAlpha3DMu_InDeltaPnTwoDPlot[inte][WhichDeltaPn], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPn+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPn) ) );

		}

		for (int WhichDeltaPT = 0; WhichDeltaPT < TwoDNBinsDeltaPT; WhichDeltaPT++) {

		  Reweight(TrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPT+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPT) ) );

		  Reweight(NoFSITrueFineBinDeltaAlphaT_InDeltaPTTwoDPlot[inte][WhichDeltaPT], 1. / ( TwoDArrayNBinsDeltaPn.at(WhichDeltaPT+1) - TwoDArrayNBinsDeltaPn.at(WhichDeltaPT) ) );

		}

		for (int WhichDeltaAlpha3Dq = 0; WhichDeltaAlpha3Dq < TwoDNBinsDeltaAlpha3Dq; WhichDeltaAlpha3Dq++) {

		  Reweight(TrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq], 1. / ( TwoDArrayNBinsDeltaAlpha3Dq.at(WhichDeltaAlpha3Dq+1) - TwoDArrayNBinsDeltaAlpha3Dq.at(WhichDeltaAlpha3Dq) ) );

		  Reweight(NoFSITrueFineBinDeltaPn_InDeltaAlpha3DqTwoDPlot[inte][WhichDeltaAlpha3Dq], 1. / ( TwoDArrayNBinsDeltaAlpha3Dq.at(WhichDeltaAlpha3Dq+1) - TwoDArrayNBinsDeltaAlpha3Dq.at(WhichDeltaAlpha3Dq) ) );

		}

		for (int WhichDeltaAlpha3DMu = 0; WhichDeltaAlpha3DMu < TwoDNBinsDeltaAlpha3DMu; WhichDeltaAlpha3DMu++) {

		  Reweight(TrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu], 1. / ( TwoDArrayNBinsDeltaAlpha3DMu.at(WhichDeltaAlpha3DMu+1) - TwoDArrayNBinsDeltaAlpha3DMu.at(WhichDeltaAlpha3DMu) ) );

		  Reweight(NoFSITrueFineBinDeltaPn_InDeltaAlpha3DMuTwoDPlot[inte][WhichDeltaAlpha3DMu], 1. / ( TwoDArrayNBinsDeltaAlpha3DMu.at(WhichDeltaAlpha3DMu+1) - TwoDArrayNBinsDeltaAlpha3DMu.at(WhichDeltaAlpha3DMu) ) );

		}

		for (int WhichDeltaAlphaT = 0; WhichDeltaAlphaT < TwoDNBinsDeltaAlphaT; WhichDeltaAlphaT++) {

		  Reweight(TrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT], 1. / ( TwoDArrayNBinsDeltaAlphaT.at(WhichDeltaAlphaT+1) - TwoDArrayNBinsDeltaAlphaT.at(WhichDeltaAlphaT) ) );

		Reweight(NoFSITrueFineBinDeltaPT_InDeltaAlphaTTwoDPlot[inte][WhichDeltaAlphaT], 1. / ( TwoDArrayNBinsDeltaAlphaT.at(WhichDeltaAlphaT+1) - TwoDArrayNBinsDeltaAlphaT.at(WhichDeltaAlphaT) ) );

		}
	       

		//----------------------------------------//

	} // End of the loop over the interaction processes		

	//----------------------------------------//		

	// index 0 = all events
	int NBins = TrueFineBinDeltaPTPlot[0]->GetXaxis()->GetNbins();

	// Loop over all the non-QE processes and merge the plots
	// Clone plot with index 2 = MEC
	TH1D* DeltaPTBkg = (TH1D*)(TrueFineBinDeltaPTPlot[2]->Clone());
	TH1D* DeltaPnBkg = (TH1D*)(TrueFineBinDeltaPnPlot[2]->Clone());

	// Add RES, DIS, COH
	for (int inte = 3; inte < NInte; inte++) {
	  
	  DeltaPTBkg->Add(TrueFineBinDeltaPTPlot[inte]);
	  DeltaPnBkg->Add(TrueFineBinDeltaPnPlot[inte]);

	}

	double cutvalue[NBins + 1];

	double signalEffPT[NBins + 1];
	double bkgEffPT[NBins + 1];
	double productPT[NBins + 1];
	double purityPT[NBins + 1];

	double signalEffPn[NBins + 1];
	double bkgEffPn[NBins + 1];
	double productPn[NBins + 1];
	double purityPn[NBins + 1];

	for (int ibin = 0; ibin < NBins; ibin++) {

	  cutvalue[ibin] = TrueFineBinDeltaPTPlot[0]->GetBinLowEdge(ibin+1);

	  signalEffPT[ibin] = TrueFineBinDeltaPTPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPTPlot[1]->Integral() ; // QE 
	  bkgEffPT[ibin] = 1. - DeltaPTBkg->Integral(0,ibin) / DeltaPTBkg->Integral() ; // bkg: MEC, RES, DIS, COH
	  productPT[ibin] = signalEffPT[ibin] * bkgEffPT[ibin] ;
	  purityPT[ibin] = TrueFineBinDeltaPTPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPTPlot[0]->Integral(0,ibin) * 100.; // QE purity (%)

	  signalEffPn[ibin] = TrueFineBinDeltaPnPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPnPlot[1]->Integral() ; // QE 
	  bkgEffPn[ibin] = 1. - DeltaPnBkg->Integral(0,ibin) / DeltaPnBkg->Integral() ; // bkg: MEC, RES, DIS, COH
	  productPn[ibin] = signalEffPn[ibin] * bkgEffPn[ibin] ;
	  purityPn[ibin] = TrueFineBinDeltaPnPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPnPlot[0]->Integral(0,ibin) * 100.; // QE purity (%)

	  // DeltaPn
	  if ( TMath::Abs(cutvalue[ibin] - 0.3) < 0.019 ) { 

	    cout << cutvalue[ibin] << " / 0.3 GeV/c DeltaPn cut event fractional loss = " << 1 - (TrueFineBinDeltaPnPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPnPlot[1]->Integral() )  << endl; 

	    cout << cutvalue[ibin] << " / 0.3 GeV/c DeltaPn cut purity = " << purityPn[ibin] << " %" << endl; 

	  }

	  
	  // DeltaPT
	  if ( TMath::Abs(cutvalue[ibin] - 0.2) < 0.019 ) { 

	    cout << cutvalue[ibin] << " / 0.2 GeV/c DeltaPT cut event fractional loss = " << 1 - (TrueFineBinDeltaPTPlot[1]->Integral(0,ibin) / TrueFineBinDeltaPTPlot[1]->Integral() )  << endl; 

	    cout << cutvalue[ibin] << " / 0.3 GeV/c DeltaPT cut purity = " << purityPT[ibin] << " %" << endl; 

	  }

	}

	signalEffPT[NBins] = 1.;
	bkgEffPT[NBins] = 0.;
	productPT[NBins] = 0.;

	signalEffPn[NBins] = 1.;
	bkgEffPn[NBins] = 0.;
	productPn[NBins] = 0.;

	TGraph* gDeltaPT = new TGraph(NBins, signalEffPT, bkgEffPT);
	TGraph* gDeltaPn = new TGraph(NBins, signalEffPn, bkgEffPn);

	TGraph* gCutDeltaPT = new TGraph(NBins, cutvalue, productPT);
	TGraph* gCutDeltaPn = new TGraph(NBins, cutvalue, productPn);

	file->cd();
	gDeltaPT->Write("gDeltaPT");
	gDeltaPn->Write("gDeltaPn");

	gCutDeltaPT->Write("gCutDeltaPT");
	gCutDeltaPn->Write("gCutDeltaPn");

	//----------------------------------------//		
		
	file->cd();
	file->Write();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created" << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	OutputTxt.close();

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
