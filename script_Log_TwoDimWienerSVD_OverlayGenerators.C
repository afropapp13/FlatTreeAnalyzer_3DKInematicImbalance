{

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Util.C");

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx");

	gROOT->ProcessLine(".L LogTwoDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("LogTwoDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("LogTwoDimWienerSVD_OverlayGenerators(false,true)");

}
