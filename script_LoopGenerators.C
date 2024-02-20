{
	vector<TString> WhichSample; vector<TString> WhichName;

	//----------------------------------------//
			
	WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.genie_v3_00_06.flat.root"); WhichName.push_back("GENIE_v3_0_6");			 
	WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.genie_v2_12_10.mec.flat.root"); WhichName.push_back("GENIE_v2_12_10_MEC");
	WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.nuwro_19_02_1.flat.root"); WhichName.push_back("NuWro_19_02_1");
	WhichSample.push_back("/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.neut_5_4_0_1.flat.root"); WhichName.push_back("NEUT_5_4_0_1");		
        WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GiBUU_Samples/GiBUU2023_Patch1/GiBUU2023.flat.root"); WhichName.push_back("GiBUU_2023");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_EffSF.flat.root"); WhichName.push_back("GENIE_v3_0_6_EffSF");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_RFG.flat.root"); WhichName.push_back("GENIE_v3_0_6_RFG");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_NoCoulomb.flat.root"); WhichName.push_back("GENIE_v3_0_6_NoCoulomb");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_hN2018.flat.root"); WhichName.push_back("GENIE_v3_0_6_hN2018");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_2_0_G18_10d_02_11a_CC/14_1000180400_CC_v3_2_0_G18_10d_02_11a.flat.root"); WhichName.push_back("GENIE_v3_2_0_G4");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_NoRPA.flat.root"); WhichName.push_back("GENIE_v3_0_6_NoRPA");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/gntp.0.gprep_G21_11b_00_000.flat.root"); WhichName.push_back("GENIE_v3_0_6_G21_11b_00_000");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_0_AR23_20i_00_000_CC/14_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root"); WhichName.push_back("GENIE_v3_4_0_AR23_20i_00_000");
	WhichSample.push_back("/pnfs/uboone/persistent/users/apapadop/GENIETweakedSamples/v3_4_0_G18_10a_02_11a_CC_HondaAtmo/14_1000180400_CC_v3_4_0_G18_10a_02_11a.flat.root"); WhichName.push_back("GENIE_v3_4_0_G18_10a_02_11a_Honda");


		
	//----------------------------------------//

        gROOT->ProcessLine(".L ../myClasses/Tools.cxx+");
        gROOT->ProcessLine(".L ../myClasses/STV_Tools.cxx+");
	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx+");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

	}
	//gROOT->ProcessLine(".q");
};
