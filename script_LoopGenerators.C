{
	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
		
	//WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6");			 
	WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.nuwro_19_02_1.flat.root"); WhichName.push_back("NuWro_19_02_1");
	//WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.neut_5_4_0_1.flat.root"); WhichName.push_back("NEUT_5_4_0_1");		
        //WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GiBUU_Samples/GiBUU2023_Patch1/GiBUU2023.flat.root"); WhichName.push_back("GiBUU_2023 Patch 1");	

	//----------------------------------------//

        gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx+");
        gROOT->ProcessLine(".L /uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
