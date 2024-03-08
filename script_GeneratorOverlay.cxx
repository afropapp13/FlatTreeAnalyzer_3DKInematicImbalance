{

	gROOT->ProcessLine(".L GeneratorOverlay.cxx+");
	gROOT->ProcessLine("GeneratorOverlay()");
	gROOT->ProcessLine("GeneratorOverlay(\"GENIE\")");
	gROOT->ProcessLine("GeneratorOverlay(\"GiBUU\")");	
	gROOT->ProcessLine("GeneratorOverlay(\"Closure\")");
	gROOT->ProcessLine("GeneratorOverlay(\"Weights\")");
	gROOT->ProcessLine("GeneratorOverlay(\"DUNE\")");
	gROOT->ProcessLine("GeneratorOverlay(\"EnergyIndependence\")");
	gROOT->ProcessLine("GeneratorOverlay(\"ANL_SF\")");
	gROOT->ProcessLine("GeneratorOverlay(\"Gv2\")");



};
