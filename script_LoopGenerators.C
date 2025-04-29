{
	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
	 	
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_6_0_AR23_20i_00_000/dune_atmospherics/14_1000180400_CC_v3_6_0_AR23_20i_00_000.flat.root"); WhichName.push_back("atmo_AR23_20i_00_000");	
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_2_AR23_20i_00_000_CC/14_1000180400_CC_v3_4_2_AR23_20i_00_000.flat.root"); WhichName.push_back("bnb_AR23_20i_00_000");			 

	//----------------------------------------//

    gROOT->ProcessLine(".L Tools.cxx+");
    gROOT->ProcessLine(".L STV_Tools.cxx+");
	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}

	//gROOT->ProcessLine(".q");

};
