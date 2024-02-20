{

	gROOT->ProcessLine(".L GeneratorOverlay.cxx+");
	gROOT->ProcessLine("GeneratorOverlay()");
	gROOT->ProcessLine("GeneratorOverlay(\"NuclearModel\")");
	gROOT->ProcessLine("GeneratorOverlay(\"FSI\")");
	gROOT->ProcessLine("GeneratorOverlay(\"NuclearEffects\")");
	gROOT->ProcessLine("GeneratorOverlay(\"Honda\")");

};
