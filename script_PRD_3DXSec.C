#include <iostream>
#include <vector>

void script_PRD_3DXSec() {

	gROOT->ProcessLine(".L /exp/uboone/app/users/apapadop/uboonecode_v08_00_00_52/srcs/ubana/ubana/myClasses/Tools.cxx++");
	gROOT->ProcessLine(".L PRD_3DXSec_Ratio.cxx++");
	gROOT->ProcessLine("PRD_3DXSec_Ratio()");

//	gROOT->ProcessLine(".q");		

}
