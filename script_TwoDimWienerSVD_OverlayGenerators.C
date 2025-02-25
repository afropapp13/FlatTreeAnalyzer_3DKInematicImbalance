{

	gROOT->ProcessLine(".L Util.C");

	gROOT->ProcessLine(".L Tools.cxx");

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
