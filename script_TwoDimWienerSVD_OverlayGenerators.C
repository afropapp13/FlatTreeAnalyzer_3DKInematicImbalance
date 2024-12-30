{

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Util.C");

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx");

	gROOT->ProcessLine(".L TwoDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,true)");
	//ANL SF
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,true)");
	//Closure test
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,true)");
	//GiBUU
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,true)");

}
