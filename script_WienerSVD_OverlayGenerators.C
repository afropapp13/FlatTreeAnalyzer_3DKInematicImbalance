{

	gROOT->ProcessLine(".L Util.C");
	gROOT->ProcessLine(".L WienerSVD_OverlayGenerators.cpp");
	//GENIE versions
	gROOT->ProcessLine("WienerSVD_OverlayGenerators()");
	//AltGen
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,true)");
	//ACHILLES
	//gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,true)");
	//ANL_SF
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,true)");
	//Closure
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,true)");
	//GiBUU
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,true)");
	//G24
	gROOT->ProcessLine("WienerSVD_OverlayGenerators(false,false,false,false,false,false,true)");
}
