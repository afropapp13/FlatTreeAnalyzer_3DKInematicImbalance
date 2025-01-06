{

	gROOT->ProcessLine(".L GeneratorOverlay.cxx+");
	gROOT->ProcessLine("GeneratorOverlay()");
	gROOT->ProcessLine("GeneratorOverlay(\"NuclearModel\")");
	gROOT->ProcessLine("GeneratorOverlay(\"FSI\")");
	gROOT->ProcessLine("GeneratorOverlay(\"FSIAR23\")");
	gROOT->ProcessLine("GeneratorOverlay(\"Honda\")");
	gROOT->ProcessLine("GeneratorOverlay(\"HondaRw\")");
	gROOT->ProcessLine("GeneratorOverlay(\"Tune\")");
	gROOT->ProcessLine("GeneratorOverlay(\"flux\")");
	gROOT->ProcessLine("GeneratorOverlay(\"neut\")");

};
