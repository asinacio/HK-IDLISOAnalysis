// Code for a simple chi2 fit to calibration data to extract extinction lengths

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "toml/toml_helper.h"

#include "DataPMT.hh"
#include "DataRunInfo.hh"

const std::string TAG = "[opticalAnalysis]: ";
const std::string ERR = "[ERROR]: ";
const std::string WAR = "[WARNING]: ";

using namespace std;

double CalcChi2( std::vector<double> wnpe, std::vector<double> d, double intensity, double proposedExL ){

  double chi2 = 0;

  for( int i = 0; i < (int)wnpe.size(); i++ ){
    double model = intensity * exp(-d[i] / proposedExL); // not considering the PMT angular responses yet

    double value = (wnpe[i]-model)*(wnpe[i]-model)/model;

    chi2 = chi2 + value;
  }

  return chi2;

}

void SendHelp();

int main(int argc, char* argv[]){

  cout << "\n";
  cout << "###########################" << endl;
  cout << "###### HK Optics Fit ######" << endl;
  cout << "###########################" << endl;
  cout << "\n";

  
  char *toml = NULL;
  int opt;
  while((opt = getopt(argc, argv, ":t:h")) != -1){
    switch(opt){
    case 't':
      toml = optarg;
      break;
    case '?':
      printf("%sUnknown option: -%c\n", WAR.c_str(), optopt);
      printf("%sPass -h for help.\n", WAR.c_str());
      return 0;
    case ':':
      printf("%sOption -%c requires argument.\n", WAR.c_str(), optopt);
      printf("%sPass -h for help.\n", WAR.c_str());
      return 0;
    case 'h':
      SendHelp();
    default:
      return 0;
    }
  }
  if(toml == NULL){
    std::cout << ERR << "Config file required. Exiting.." << std::endl;
    return -1;
  }
  //Search for file in the /inputs directory
  std::string tomlFile = toml;
  ifstream tomlStream(tomlFile.c_str());
  if(!tomlStream.good()){
    std::cout << ERR << "Cannot find config file: " << tomlFile << std::endl;
    std::cout << ERR << "Check the file path. Exiting.." << std::endl;
    return -1;
  }
    
  const toml::value config = toml::parse(tomlFile);

  const auto inFile = toml::find<std::string>(config, "input");

  std::string inFileName = std::getenv("ANAUP") + std::string("/inputs/") + std::string(inFile);
  std::cout << TAG << "Reading input file: " << inFileName << std::endl;

  DataRunInfo *data = new DataRunInfo();

  data->FillRunInfo( inFileName );
  data->FillPMTInfo( inFileName );
  data->FillmPMTInfo( inFileName );

  // Prints for debugging
  cout << data->GetLambda() << "  " << data->GetNPhotons() << endl;

  map< Int_t, DataPMT >::iterator iLP;
  for ( iLP = data->GetPMTIterBegin(); iLP != data->GetPMTIterEnd(); iLP++ ){
    cout << iLP->first << "  " << ( iLP->second ).GetPMTID() << "  " << ( iLP->second ).GetPMTSourceDist() << "  " << endl;
    // "first" accesses the first element of the map, which is the int of the PMT ID. "second" accesses the PMT object for that ID
  }

  /*
  // Open and read data files
  TFile *inFile = new TFile( inFileName.c_str(), "read" );

  // If files does not exist, print a message and exit
  if( !inFile || inFile->IsZombie() ){
    std::cerr << "Error opening file!" << endl;
    exit(-1);
  }

  cout << "Reading file " << inFileName << "..." << endl;
  inFile->ls(); // print file structure, for debugging*/

  /*
  // Loading the ttrees
  TTree *tsource = (TTree*)inFile->Get("source");
  TTree *tpmt_geom = (TTree*)inFile->Get("pmt_geom");
  TTree *tmpmt_geom = (TTree*)inFile->Get("mpmt_geom");
  TTree *tpmt_hits = (TTree*)inFile->Get("pmt_hits");
  TTree *tmpmt_hits = (TTree*)inFile->Get("mpmt_hits");

  // tpmt_hits->Print(); // print ttree structure, for debugging

  // Loading information about the laser and fibre
  cout << "Loading calibration source information..." << endl;
  double sourcePosx, sourcePosy, sourcePosz, sourceWl, sourceI, sourceAng, groupVel;
  tsource->SetBranchAddress("vtx_x", &sourcePosx);  // Source x coordinate
  tsource->SetBranchAddress("vtx_y", &sourcePosy);  // Source y coordinate
  tsource->SetBranchAddress("vtx_z", &sourcePosz);  // Source z coordinate
  tsource->SetBranchAddress("lambda", &sourceWl);   // Source wavelength [nm]
  tsource->SetBranchAddress("nPhotons", &sourceI);  // Source intensity, number of photons simulated
  tsource->SetBranchAddress("opAng", &sourceAng);   // Source opening angle
  tsource->SetBranchAddress("vg", &groupVel);       // Photon group velocity [nm s-1]

  //Access the single event
  tsource->GetEntry(0);

  // Loading PMT information
  cout << "Loading PMT data..." << endl;
  double dist, xpos, ypos, zpos, costh, cosths, phis, omega, rad, eff;
  int pmtID, mpmtID, isPMTOn;
  tpmt_geom->SetBranchAddress("R", &dist);          // Distance to PMT from source
  tpmt_geom->SetBranchAddress("x", &xpos);          // PMT x coordinate
  tpmt_geom->SetBranchAddress("y", &ypos);          // PMT y coordinate
  tpmt_geom->SetBranchAddress("z", &zpos);          // PMT z coordinate
  tpmt_geom->SetBranchAddress("costh", &costh);     // Photon angle relative to PMT
  tpmt_geom->SetBranchAddress("cosths", &cosths);   // PMT theta relative to source
  tpmt_geom->SetBranchAddress("phis", &phis);       // PMT phi relative to source
  tpmt_geom->SetBranchAddress("omega", &omega);     // Solid angle subtended
  tpmt_geom->SetBranchAddress("radius", &rad);      // PMT radius
  tpmt_geom->SetBranchAddress("eff", &eff);         // PMT efficiency
  tpmt_geom->SetBranchAddress("pmtID", &pmtID);     // PMT ID (enum. from 0)
  tpmt_geom->SetBranchAddress("isPMTOn", &isPMTOn); // PMT status (0 = off, 1 = on) - can add extra status tags if needed

  double nHits, nPE, time_tofcorr, time;
  tpmt_hits->SetBranchAddress("nHits", &nHits);
  tpmt_hits->SetBranchAddress("nPE", &nPE);
  tpmt_hits->SetBranchAddress("time", &time);
  tpmt_hits->SetBranchAddress("time_tofcorr", &time_tofcorr);
  tpmt_hits->SetBranchAddress("pmtID", &pmtID);

  // For testing
  /*cout << "Size tpmt_geom = " << tpmt_geom->GetEntries() << ", Size tpmt_hits = " << tpmt_geom->GetEntries() << endl;

  TH1* h1 = new TH1I("h1", "h1 title", 100, 0.0, 4000.0);
  for(int i = 0; i<tpmt_hits->GetEntries(); i++){
    tpmt_hits->GetEntry(i);
    h1->Fill(time);
  }
  h1->Draw();*/

  // After loading the data, probably want to apply PMT selection cuts
  // For example to exclude PMTs offline, with bad timing/electronics calibration, etc
  /*  cout << "Applying PMT selection cuts..." << endl;

  // Fit Model: NPE = I0 * Omega * A(theta) * exp( -d/L_alpha )
  // With I0 being the diffuser intensity entering the detector
  // Omega is the solid angle subtended by the PMT -- in principle this can be calculated apriori per PMT and used to weigh the measured NPE
  // A(theta) is the angular response of the PMT
  // d is the distance from the diffuser to the PMT
  // L_alpha is the extinction/attenuation length of the detector medium
 
  // Measurements: 
  // -- Total extinction: could be determined by specifically looking at the variation of intensity in the beam spot (NOTE: could also be determined from a 2D minimisation of the NPE and time, keeping in mind that only some of the PMTs see direct light, and other just see exclusively scattered light )
  // -- Scattering: determined by looking at regions away from the beamspot, where only scattered light could have arrived

  // Simple fit to start with: only fit attenuations
  std::vector<double> wNPE; // NPE corrected by solid angle
  std::vector<double> d; // Distance source-PMT
  for(int i = 0; i<tpmt_hits->GetEntries(); i++){
    tpmt_hits->GetEntry(i);
    double npe = nPE;
    tpmt_geom->GetEntry(i);
    npe = npe / omega;
    
    wNPE.push_back(npe);
    d.push_back(dist);
  }

  // Basic Chi2 fit
  // Scanning extinction lengths that minimize Chi2
  double minExL = 0.0;
  double maxExL = 10000.0; // in cm (?)
  double bestChi2 = 10e30; // start with absurdly large value for the best Chi2
  double bestExL = 0.0;
  for( int i = minExL; i < maxExL; i++ ){
    double chi2 = CalcChi2(wNPE,d,sourceI,i);
      if(chi2 < bestChi2){
        bestChi2 = chi2;
        bestExL = i;
      }
  }

  cout << "The best fit extinction length is " << bestExL << endl;*/
  
  // inFile->Close();

  return 0;

}

void SendHelp(){

  std::cout << TAG << "Main analysis program.\n"
	    << "USAGE: opticalAnalysis -t [config toml]\n"
	    << "OPTIONS:\n"
	    << "-t : Config file\n"
	    << "-h : Display this message\n";

}
