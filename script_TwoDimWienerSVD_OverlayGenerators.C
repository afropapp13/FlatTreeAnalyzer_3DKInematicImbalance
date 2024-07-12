{

	gROOT->ProcessLine(".L ../myClasses/Util.C");

	gROOT->ProcessLine(".L ../myClasses/Tools.cxx");

	gROOT->ProcessLine(".L TwoDimWienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,true)");
	//ANL SF
	//gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,true)");
	//Closure
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,true)");
	//nuclear
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,true)");
	//mec
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//gibuu
	gROOT->ProcessLine("TwoDimWienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");

}
