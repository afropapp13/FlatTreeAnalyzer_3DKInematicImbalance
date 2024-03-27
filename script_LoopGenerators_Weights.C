{

  // The T2K Tune weights will be applied

  vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
	
        WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6"); 

        WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_10a_02_11a.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_02_11a"); 

        WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_10a_02_11b.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_02_11b"); 

        WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_10a_03_320.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_03_320"); 

        WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_10a_03_330.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_03_330"); 

	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_10a_02_11a_RS.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_02_11a_RS");
	
	//WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.ghep_v3_4_0_G18_02a_02_11a.prep.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_02a_02_11a");
	
	//----------------------------------------//

        gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx+");
        gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\",\"Weights\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
