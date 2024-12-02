{

  // The T2K Tune weights will be applied

  vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
	
        WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6"); 
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_2_G18_10a_02_11a/14_1000180400_CC_v3_4_2_G18_10a_02_11a.flat.root"); WhichName.push_back("GENIE_v3_4_2_G18_10a_02_11a");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_2_G18_10b_02_11a/14_1000180400_CC_v3_4_2_G18_10b_02_11a.flat.root"); WhichName.push_back("GENIE_v3_4_2_G18_10b_02_11a");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_2_G18_10d_02_11a/14_1000180400_CC_v3_4_2_G18_10d_02_11a.flat.root"); WhichName.push_back("GENIE_v3_4_2_G18_10d_02_11a");

	//----------------------------------------//

        gROOT->ProcessLine(".L ../myClasses/Tools.cxx+");
        gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\",\"Weights\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
