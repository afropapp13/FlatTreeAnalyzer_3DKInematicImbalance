{

	gROOT->ProcessLine(".L ../myClasses/Util.C");
	gROOT->ProcessLine(".L WienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("WienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,true)");
	//Closure test
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,true)");
	//NuWro closure test
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,true)");
	//sf noah
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//nuclear
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false, false,true)");
	//mec
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false, false,true)");
	//gibuu
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,false, false,true)");
	//tune fsi
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,false, false,false,true)");
	//g18t vs ar23
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,false, false,false,false,true)");



}
