#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TString.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>

#include "constants.h"

using namespace std;
using namespace constants;

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 1) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());

}

//----------------------------------------//		

void divide_bin_width(TH1D* h, double SF=1.) {

	int NBins = h->GetXaxis()->GetNbins();
  
	for (int i = 0; i < NBins; i++) {
  
	  double CurrentEntry = h->GetBinContent(i+1);
	  double NewEntry = SF * CurrentEntry / h->GetBinWidth(i+1);
  
	  double CurrentError = h->GetBinError(i+1);
	  double NewError = SF * CurrentError / h->GetBinWidth(i+1);
  
	  h->SetBinContent(i+1,NewEntry); 
	  h->SetBinError(i+1,NewError); 
	  //h->SetBinError(i+1,0.000001); 
  
	}
  
  }

  //----------------------------------------//	
