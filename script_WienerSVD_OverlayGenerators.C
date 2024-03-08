{

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Util.C");
	gROOT->ProcessLine(".L WienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("WienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,true)");
	//ACHILLES
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,true)");
	//ANL_SF
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,true)");
}
