//********************************************
// Script to extract information from WCSim
// root files and put them into library-
// independent trees.
//
// Partially based on code from Ka Ming's
// optical fit.
//
// TreeConverter -h for more info.
//
// s.j.jenkins@liverpool.ac.uk
//********************************************



#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

#include "CalcGroupVelocity.h"

const std::string TAG = "[TreeConverter]: ";
const std::string ERR = "[ERROR]: ";
const std::string WAR = "[WARNING]: ";
bool verbose = false;

void ParseMacro(std::ifstream &file, double p[]);
double GetSolidAngle(double r, double R); //simple version, look into what Alie is using
void SendHelp();

int main(int argc, char *argv[]){

  //Command line arguments
  char *inFileName = NULL;
  char *outFileName = NULL;
  char *macroFileName = NULL;
  double wavelength = 385.;
  bool hybrid = true;
  bool isDiff = true; //currently only want to use diffuser
                      //should make this depend on opAng

  //Additional useful constants
  const int nPMTpermPMT = 19;
  const int nPMTTypes = 2;
  double pmtRadius[2];

  int opt;
  while ((opt = getopt(argc, argv, ":i:o:m:bvh")) != -1){
    switch(opt){
    case 'i':
      inFileName = optarg;
      break;
    case 'o':
      outFileName = optarg;
      break;
    case 'm':
      macroFileName = optarg;
      break;
    case 'b':
      hybrid = false;
      break;
    case 'v':
      verbose = true;
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

  //Open the input file, checking it exists
  if(inFileName == NULL){
    std::cout << ERR << "No input file found. Exiting.." << std::endl;
    return -1;
  }
  TFile *inFile = new TFile(inFileName, "READ");
  if(!inFile->IsOpen()){
    std::cout << ERR << "Cannot open input file: " << inFileName << std::endl;
    return -1;
  }

  //Open the macro file, checking it exists, and is a .mac file
  if(macroFileName == NULL){
    std::cout << ERR << "No macro file found. Required for source information. Exiting.." << std::endl;
    return -1;
  }
  std::string mString(macroFileName);
  if (mString.find(".mac") == std::string::npos){
    std::cout << ERR << "Supplied macro file is not .mac format. Please provide a .mac file. Exiting.." << std::endl;
    return -1;
  }
  std::ifstream macroFile;
  macroFile.open(macroFileName);
  if(!macroFile){
    std::cout << ERR << "Supplied macro file cannot be opened. Exiting.." << std::endl;
    return -1;
  }
  //Read the information required for source tree
  std::cout << TAG << "Parsing input macro file: " << macroFileName << std::endl;
  double params[3] = {-1., -1., -1}; //Nphotons, opening angle, wavelength
  ParseMacro(macroFile, params);
  int nPhotons = (int)params[0];
  double opAng = params[1];
  wavelength = params[2];


  //Open the output file
  if(outFileName == NULL) {
    std::cout << WAR << "No output filename provided. Falling back to default." << std::endl;
    outFileName = (char*)"output.root";
  }
  TFile *outFile = new TFile(outFileName, "RECREATE");
  std::cout << TAG << outFileName << " opened." << std::endl;
  

  //Access the geometry
  //This only needs one 'event'
  TTree *wcsimGeoT = (TTree*)inFile->Get("wcsimGeoT");
  WCSimRootGeom* geom = new WCSimRootGeom();
  wcsimGeoT->SetBranchAddress("wcsimrootgeom", &geom);
  wcsimGeoT->GetEntry(0);
  std::cout << TAG << "Geometry loaded." << std::endl;

  //Access the options
  //This only needs one 'event'
  TTree *wcsimOptT = (TTree*)inFile->Get("wcsimRootOptionsT");
  WCSimRootOptions *options = 0;
  wcsimOptT->SetBranchAddress("wcsimrootoptions", &options);
  wcsimOptT->GetEntry(0);

  double rayff = options->GetRayff();
  double bsrff = options->GetBsrff();
  double abwff = options->GetAbwff();
  double rgcff = options->GetRgcff();
  double mieff = options->GetMieff();
  double qeff  = options->GetQeff();

  std::cout << TAG << "Options loaded." << std::endl;
  if(verbose) options->Print();

  //Finally, access the event trees
  //Need two of these, one for mPMTs
  TTree *wcsimT = (TTree*)inFile->Get("wcsimT");
  WCSimRootEvent* wcsimrootsuperevent  = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  wcsimT->SetBranchAddress("wcsimrootevent", &wcsimrootsuperevent);
  wcsimT->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
  if(hybrid){
    wcsimT->SetBranchAddress("wcsimrootevent2", &wcsimrootsuperevent2);
    wcsimT->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }
  std::cout << TAG << "Trees loaded." << std::endl;

  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  //Save the tuning parameters
  TTree *tuning = new TTree("tuning", "tuning");
  tuning->Branch("rayff", &rayff);
  tuning->Branch("bsrff", &bsrff);
  tuning->Branch("abwff", &abwff);
  tuning->Branch("rgcff", &rgcff);
  tuning->Branch("mieff", &mieff);
  tuning->Branch("qeff",  &qeff);
  tuning->Fill();
  outFile->cd();
  tuning->Write();

  //Save the PMT geometry
  double dist, xpos, ypos, zpos, costh, cosths, phis, omega, rad, eff;
  int pmtID, mpmtID, isPMTOn;
  //B&L PMT tree
  TTree *pmt_geom = new TTree("pmt_geom", "pmt_geom");
  pmt_geom->Branch("R", &dist);          //Distance to PMT from source
  pmt_geom->Branch("x", &xpos);          //PMT x coordinate
  pmt_geom->Branch("y", &ypos);          //PMT y coordinate
  pmt_geom->Branch("z", &zpos);          //PMT z coordinate
  pmt_geom->Branch("costh", &costh);     //Photon angle relative to PMT
  pmt_geom->Branch("cosths", &cosths);   //PMT theta relative to source
  pmt_geom->Branch("phis", &phis);       //PMT phi relative to source
  pmt_geom->Branch("omega", &omega);     //Solid angle subtended
  pmt_geom->Branch("radius", &rad);      //PMT radius
  pmt_geom->Branch("eff", &eff);         //PMT efficiency
  pmt_geom->Branch("pmtID", &pmtID);     //PMT ID (enum. from 0)
  pmt_geom->Branch("isPMTOn", &isPMTOn);  //PMT status (0 = off, 1 = on) - can add extra status tags if needed
  

  //mPMT tree - may need to add additional info here
  TTree *mpmt_geom = new TTree("mpmt_geom", "mpmt_geom");
  mpmt_geom->Branch("R", &dist);          //Distance to PMT from source
  mpmt_geom->Branch("x", &xpos);          //PMT x coordinate
  mpmt_geom->Branch("y", &ypos);          //PMT y coordinate
  mpmt_geom->Branch("z", &zpos);          //PMT z coordinate
  mpmt_geom->Branch("costh", &costh);     //Photon angle relative to PMT
  mpmt_geom->Branch("cosths", &cosths);   //PMT theta relative to source
  mpmt_geom->Branch("phis", &phis);       //PMT phi relative to source
  mpmt_geom->Branch("omega", &omega);     //Solid angle subtended
  mpmt_geom->Branch("radius", &rad);      //PMT radius
  mpmt_geom->Branch("eff", &eff);         //PMT efficiency
  mpmt_geom->Branch("pmtID", &pmtID);     //PMT ID (enum. from 0)
  mpmt_geom->Branch("isPMTOn", &isPMTOn); //PMT status (0 = off, 1 = on) - can add extra status tags if needed
  mpmt_geom->Branch("mpmtID", &mpmtID);   //PMT ID within module


  //Injector position
  double vtx[3];
  wcsimT->GetEntry(0);
  wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

  for(int i=0; i<3; i++) vtx[i] = wcsimrootevent->GetVtx(i);

  //Calculate photon group velocity
  double vg = CalcGroupVelocity(wavelength);
  vg /= 1E9; //mns-1

  //Fill source tree
  TTree *source_tree = new TTree("source", "source");
  source_tree->Branch("vtx_x", &vtx[0]);       //Injector x position
  source_tree->Branch("vtx_y", &vtx[1]);       //Injector y position
  source_tree->Branch("vtx_z", &vtx[2]);       //Injector z position
  source_tree->Branch("lambda", &wavelength);  //Wavelength [nm]
  source_tree->Branch("vg", &vg);              //Photon group velocity [nms-1]
  source_tree->Branch("nPhotons", &nPhotons);  //Initial simulated photons
  source_tree->Branch("opAng", &opAng);        //Injector simulated opening angle

  source_tree->Fill();
  outFile->cd();
  source_tree->Write();

  
  //Set up injector local coordinate systems
  //Injection is inward perpendicular to wall
  //May switch these to TVector3s once I'm sure everything works correctly
  double vSourceDir[3];
  double vSourceXAxisLocal[3], vSourceYAxisLocal[3];
  double norm;
  double cut = geom->GetWCCylLength()/2*0.9;
  //Barrel injectors
  if(abs(vtx[2]) < cut){
    vSourceDir[0] = vtx[0];
    vSourceDir[1] = vtx[1];
    vSourceDir[2] = 0;
    norm = sqrt(vtx[0]*vtx[0] + vtx[1]+vtx[1]);
    //Want unit vectors - direct inwards
    vSourceDir[0] /= norm;
    vSourceDir[1] /= norm;
    //Construct local axes - injection is along
    //the z axis in local system
    vSourceXAxisLocal[0] = 0;
    vSourceXAxisLocal[1] = 0;
    vSourceXAxisLocal[2] = 1;
    vSourceYAxisLocal[0] = -vtx[1]/norm;
    vSourceYAxisLocal[1] = vtx[0]/norm;
    vSourceYAxisLocal[2] = 0;
  }
  //Endcaps
  else{
    vSourceDir[0] = 0;
    vSourceDir[1] = 0;
    //Top
    if(vtx[2] > cut){
      vSourceDir[2]        = -1;
      vSourceXAxisLocal[0] =  1;
      vSourceXAxisLocal[1] =  0;
      vSourceXAxisLocal[2] =  0;
      vSourceYAxisLocal[0] =  0;
      vSourceYAxisLocal[1] = -1;
      vSourceYAxisLocal[2] =  0;
    }
    else{
      vSourceDir[2]        = 1;
      vSourceXAxisLocal[0] = 1;
      vSourceXAxisLocal[1] = 0;
      vSourceXAxisLocal[2] = 0;
      vSourceYAxisLocal[0] = 0;
      vSourceYAxisLocal[1] = 1;
      vSourceYAxisLocal[2] = 0;
    }
  }

  std::cout << TAG << "Reading PMT information." << std::endl;

  int nPMTs_BL = geom->GetWCNumPMT();
  int nPMTs_m  = 0;
  if (hybrid) nPMTs_m = geom->GetWCNumPMT(true);

  if(verbose){
    printf("%s%i B&L PMTs.\n", TAG.c_str(), nPMTs_BL);
    printf("%s%i mPMTs.\n", TAG.c_str(), nPMTs_m);
  }

  pmtRadius[0] = geom->GetWCPMTRadius();
  pmtRadius[1] = geom->GetWCPMTRadius(true);
  
  double PMTpos[3], vOrient[3], vDir[3];
  double vecNorm, orientNorm, localx, localy;
  //Loop over pmts
  int nPMTs;
  for (int pmtType = 0; pmtType < nPMTTypes; pmtType++){

    nPMTs = pmtType == 0 ? nPMTs_BL : nPMTs_m;
    for (int i=0; i<nPMTs; i++){
      
      WCSimRootPMT pmt;
      pmt = geom->GetPMT(i, pmtType==1);
      pmtID = i;
      mpmtID = pmtType==0 ? 0 : pmtID%nPMTpermPMT;

      for(int j=0; j<3; j++){
	PMTpos[j]  = pmt.GetPosition(j);
	vOrient[j] = pmt.GetOrientation(j); //direction PMT faces
	vDir[j] = PMTpos[j] - vtx[j]; //vector from inj to PMT
      }
	
      xpos = PMTpos[0]; //For tree outputs
      ypos = PMTpos[1];
      zpos = PMTpos[2];

      //Calculate angles etc.
      vecNorm = sqrt(vDir[0]*vDir[0] + vDir[1]*vDir[1] + vDir[2]*vDir[2]);
      orientNorm = sqrt(vOrient[0]*vOrient[0] + vOrient[1]*vOrient[1] + vOrient[2]*vOrient[2]);
      for(int j=0; j<3; j++){
	vDir[j] /= vecNorm;
	vOrient[j] /= orientNorm;
      }
      dist = vecNorm;

      costh  = -1*(vDir[0]*vOrient[0] + vDir[1]*vOrient[1] + vDir[2]*vOrient[2]);
      cosths = vDir[0]*vSourceDir[0] + vDir[1]*vSourceDir[1] + vDir[2]*vSourceDir[2];
      localx = vDir[0]*vSourceXAxisLocal[0] + vDir[1]*vSourceXAxisLocal[1] + vDir[2]*vSourceXAxisLocal[2];
      localy = vDir[0]*vSourceYAxisLocal[1] + vDir[1]*vSourceYAxisLocal[1] + vDir[2]*vSourceYAxisLocal[2];
      phis   = atan2(localy, localx);
      
      rad = pmtType==0 ? pmtRadius[0] : pmtRadius[1];
      omega  = GetSolidAngle(rad, dist);

      //is PMT working - currently all marked as good
      isPMTOn = 1;

      //PMT efficiency - currently all 100% efficient
      eff = 1.0;
    
      if(pmtType==0) pmt_geom->Fill();
      if(pmtType==1) mpmt_geom->Fill();
    }//End of loop over individual PMTs
  }//End of loop over PMT types

  outFile->cd();
  pmt_geom->Write();
  if(hybrid) mpmt_geom->Write();

  
  //Now we can access the actual hit information
  //Will probably want to add additional info here,
  //but start with these for now.
  double nHits, nPE, time_tofcorr, time;
  TTree *pmt_hits = new TTree("pmt_hits", "pmt_hits");
  pmt_hits->Branch("nHits", &nHits);
  pmt_hits->Branch("nPE", &nPE);
  pmt_hits->Branch("time", &time);
  pmt_hits->Branch("time_tofcorr", &time_tofcorr);
  pmt_hits->Branch("pmtID", &pmtID);
  
  TTree *mpmt_hits = new TTree("mpmt_hits", "mpmt_hits");
  mpmt_hits->Branch("nHits", &nHits);
  mpmt_hits->Branch("nPE", &nPE);
  mpmt_hits->Branch("time", &time);
  mpmt_hits->Branch("time_tofcorr", &time_tofcorr);
  mpmt_hits->Branch("pmtID", &pmtID);

  std::cout << TAG << "Reading hit information." << std::endl;

  
  int nCDigiHits, nCDigiHits2, pmtNum;
  double q, photonDist, tof;
  double photonDir[3];
  //Loop over events
  for(int ev=0; ev<wcsimT->GetEntries(); ev++){

    wcsimT->GetEvent(ev);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);

    nCDigiHits  = wcsimrootevent->GetNcherenkovdigihits();
    nCDigiHits2 = hybrid ? wcsimrootevent2->GetNcherenkovdigihits() : 0;

    //Loop over PMT types
    for (int pmtType=0; pmtType < nPMTTypes; pmtType++){
      nHits = (pmtType==0) ? nCDigiHits : nCDigiHits2;

      for (int ihit=0; ihit<nHits; ihit++){
	WCSimRootCherenkovDigiHit *hit;
	if(pmtType==0) hit = (WCSimRootCherenkovDigiHit*)wcsimrootevent->GetCherenkovDigiHits()->At(ihit);
	else           hit = (WCSimRootCherenkovDigiHit*)wcsimrootevent2->GetCherenkovDigiHits()->At(ihit);

	q = hit->GetQ();
	time = hit->GetT();//digitised hit time
	pmtNum = hit->GetTubeId(); //1 higher
	pmtID = pmtNum-1;
	
	WCSimRootPMT pmt = geom->GetPMT(pmtNum-1, pmtType==1);
	for(int i=0; i<3; i++) 
	  photonDir[i] = pmt.GetPosition(i) - vtx[i];

	photonDist = sqrt(photonDir[0]*photonDir[0] + photonDir[1]*photonDir[1] + photonDir[2]*photonDir[2]);
	photonDist *= 1E-2; //Convert to m
	tof = photonDist/vg;
	time_tofcorr = time - tof;
	
	nPE = q;

	if(pmtType==0) pmt_hits->Fill();
	if(pmtType==1) mpmt_hits->Fill();
      }
    }
  }

  outFile->cd();
  pmt_hits->Write();
  if (hybrid) mpmt_hits->Write();

  outFile->Close();

  std::cout << TAG << "Tree conversion complete." << std::endl;
  std::cout << "\u3042\u308A\u304C\u3068\u3046\u3054\u3056\u3044\u307E\u3057\u305F!" << std::endl;
  return 0;
}

//Function to search the input macro file for source information
//There's a fair amount of hard-coding here, but as the commands
//for the injectors in the .mac file _should_ always be the same,
//it shouldn't cause issues.
void ParseMacro(std::ifstream &file, double p[]){
  
  std::string line;
  std::string nphot, opAng, wav;
  bool pFound, aFound, wFound;
  pFound = aFound = wFound = false;
  //Number of photons - initial source intensity
  while(std::getline(file, line)){
    if (line.find("/mygen/injector_nPhotons") == std::string::npos) continue;
    else{
      nphot = line.substr(line.find("nPhotons") + 9);
      pFound = true;
      if(verbose) std::cout << TAG << "NPhotons = " << nphot << std::endl;
      break;
    }
  }

  //Opening angle
  file.clear();
  file.seekg(0);
  while(std::getline(file, line)){
    if (line.find("/mygen/injector_opening_angle") == std::string::npos) continue;
    else{
      opAng = line.substr(line.find("opening_angle") + 14);
      aFound = true;
      if(verbose) std::cout << TAG << "Opening angle = " << opAng << std::endl;
      break;
    }
  }

  //Wavelength
  file.clear();
  file.seekg(0);
  while(std::getline(file, line)){
    if (line.find("/mygen/injector_wavelength") == std::string::npos) continue;
    else{
      wav = line.substr(line.find("wavelength") + 11);
      wFound = true;
      if(verbose) std::cout << TAG << "Wavelength = " << wav << std::endl;
      break;
    }
  }
  
  if(!pFound) std::cout << WAR << "Number of photons information missing from .mac file!" << std::endl;
  else p[0] = std::stod(nphot);
  if(!aFound) std::cout << WAR << "Opening angle information missing from .mac file!" << std::endl;
  else p[1] = std::stod(opAng);
  if(!wFound) std::cout << WAR << "Wavelength information missing from .mac file!" << std::endl;
  else p[2] = std::stod(wav);  

}

//This is a simple approximation
double GetSolidAngle(double r, double R){

  return TMath::Pi()*r*r/(R*R);

}

void SendHelp(){

  std::cout << TAG << "Script to read WCSim output and extract PMT hit info.\n"
	    << "USAGE: TreeConverter -i [input file] -m [macro file]\n"
	    << "OPTIONS:\n"
	    << "-i : Input file\n"
	    << "-o : Output file\n"
	    << "-m : Macro file from simulation\n"
	    << "-b : Use B&L PMTs only\n"
	    << "-v : Print debug info\n"
	    << "-h : Display this message\n";

}
