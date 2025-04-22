#include "TFile.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <math.h>


using namespace std;

//---------------------------//

vector<string> split (string s, string delimiter) {

    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;
  
    while ( (pos_end = s.find (delimiter, pos_start)) != string::npos ) {
  
      token = s.substr (pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      res.push_back (token);
  
    }
  
    res.push_back (s.substr (pos_start));
    return res;
  
}

//---------------------------//

void convert_text_to_histo() {

    //---------------------------//

    // I/O files

    TFile *f = new TFile("honda_flux_solmax.root" ,"recreate");
    ifstream input("hms-nu-01-01-n3650.d");

    //TFile *f = new TFile("honda_flux_solmin.root" ,"recreate");
    //ifstream input("hms-nu-01-01-000.d");

    string line;
    int line_counter = 0;

    double min_e = 0.;
    double max_e = 3.;    

    //---------------------------//

    vector<double> energy;
    vector<double> flux;
    vector<double> bin_edges;    

    //---------------------------//  

    // loop over the text file

    while ( std::getline(input, line) ) { 

        //---------------------------//

        vector<string> words = split (line," ");

        double energy_line = std::stod(words[1].c_str());
        // BNB barely has any entries above 3 GeV
        if (energy_line > max_e) { continue; }
        ++line_counter;
        energy.push_back( energy_line );
        //cout << "energy_line = " << energy_line << endl;

        double flux_line = std::stod(words[2].c_str() );
        flux.push_back( flux_line );    
        //cout << "flux_line = " << flux_line << endl;    

        // lower edge
        if (line_counter == 2) { bin_edges.push_back( energy[0] - (energy[1] - energy[0])/2. ); }
        //upper edge
        if (line_counter >=2) { bin_edges.push_back( energy[line_counter-2] + (energy[line_counter-1] - energy[line_counter-2])/2. ); }

        //if (line_counter >=2) { cout << "energy.back() = " << energy.back() << "  bin_edges.back() = " << bin_edges.back() << endl; }

    } // end of the loop over the text lines

    bin_edges.push_back( energy.back() + (energy.back() - energy.at( energy.size()-2 ))/2. );

    //---------------------------//  
    
    TGraph* g = new TGraph(line_counter,&energy[0],&flux[0]);
    //g->Draw("AC*");

    //cout << "bin_edges.size() = " << bin_edges.size() << endl;
    //cout << "energy.size() = " << energy.size() << endl;     

    TH1D* h = new TH1D("h",";E_{#nu} [GeV]",bin_edges.size()-1,bin_edges.at(0),bin_edges.back());

    for (int i = 1; i <= line_counter; i++) {

        //cout << "energy.at(i-1) = " << energy.at(i-1) << endl;
        //int bin = h->FindBin(energy.at(i-1));
        //cout << "bin = " << bin << endl;

        double x,y;
        g->GetPoint(i,x,y);
        h->Fill(x,y);

    }

    //h->Draw("hist");

    f->cd();
    h->Write();

} // end of the program