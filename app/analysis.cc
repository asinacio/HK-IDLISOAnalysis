#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph2D.h"

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

  // Create data quality plots, save to output file
  TH1 *hHits = new TH1F("hHits","Integrated Hits",100, 0.0, 1000.0);
  TH1 *hPE = new TH1F("hPE","Integrated Photoelectrons",100, 0.0, 100.0);
  TH1 *htime = new TH1F("htime","Time",100, -100.0, 2000.0);
  TH1 *htimeToFCorrected = new TH1F("htimeToFCorrected","Time ToF corrected",100, -100.0, 2000.0);

  TH2 *nHitsvsSourceDist = new TH2F("nHitsvsSourceDist","Hits vs distance to source",100,0.0,1000.0,100,0.0,8000.0);
  TH2 *nHitsvsSourceCos = new TH2F("nHitsvsSourceCos","Hits vs cosTheta relative to source",100,0.0,1000.0,100,-1.0,1.0);
  TH2 *nHitsvsSourcePhi = new TH2F("nHitsvsSourcePhi","Hits vs phi relative to source",100,0.0,1000.0,100,-180.0,180.0);
  
  TGraph2D *goccupancyMap = new TGraph2D;

  int counter = 0;
  map<Int_t, DataPMT>::iterator iPMT;
  for( iPMT = data->GetPMTIterBegin(); iPMT != data->GetPMTIterEnd(); iPMT++ ){

    // Only fill plots for normal PMTs (type 1)
    if( data->GetPMT(iPMT->first).GetPMTType() != 1 ) continue;
    
    hHits->Fill( data->GetPMT(iPMT->first).GetPMTnHits() );
    hPE->Fill( data->GetPMT(iPMT->first).GetPMTnPE() );
    htime->Fill( data->GetPMT(iPMT->first).GetPMTTime() );
    htimeToFCorrected->Fill( data->GetPMT(iPMT->first).GetPMTTime_ToFCorr() );

    nHitsvsSourceDist->Fill( data->GetPMT(iPMT->first).GetPMTnHits(), data->GetPMT(iPMT->first).GetPMTSourceDist() );
    nHitsvsSourceCos->Fill( data->GetPMT(iPMT->first).GetPMTnHits(), data->GetPMT(iPMT->first).GetPMTCosTh() );
    nHitsvsSourcePhi->Fill( data->GetPMT(iPMT->first).GetPMTnHits(), data->GetPMT(iPMT->first).GetPMTPhi()*180.0/TMath::Pi() );

    // Project cylinder into 2D
    TVector3 pos = data->GetPMT(iPMT->first).GetPMTPos();
    double height = pos.Z();
    double azimuth = pos.Phi()*180.0/TMath::Pi();
    if( data->GetPMT(iPMT->first).GetPMTnHits() > 0  ){
      goccupancyMap->SetPoint(counter,azimuth,height,data->GetPMT(iPMT->first).GetPMTnHits());
      counter++;
    }
    
  }
  goccupancyMap->SetMarkerStyle(21);
  goccupancyMap->SetMarkerSize(0.5);

  std::string outrootFileName = "outputs/test.root";
  TFile *outrootFile = new TFile();
  outrootFile = TFile::Open(outrootFileName.c_str(), "RECREATE");
  outrootFile->cd();
  hHits->Write();
  hPE->Write();
  htime->Write();
  htimeToFCorrected->Write();
  nHitsvsSourceDist->Write();
  nHitsvsSourceCos->Write();
  nHitsvsSourcePhi->Write();
  goccupancyMap->Write();
  outrootFile->Close();
  
  delete outrootFile;
  delete hHits;
  delete hPE;
  delete htime;
  delete htimeToFCorrected;

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
