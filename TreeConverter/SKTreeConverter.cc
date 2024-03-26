//********************************************
// Modified script to extract information from
// SK root files and put them into library-
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
#include <sstream>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"


#include "CalcGroupVelocity.h"
#include "SKInjectorLocations.h"

const std::string TAG = "[TreeConverter]: ";
const std::string ERR = "[ERROR]: ";
const std::string WAR = "[WARNING]: ";
bool verbose = false;

double GetSolidAngle(double r, double R); //simple version, look into what Alie is using
void ReadPMTGeometry(std::ifstream &posFile, std::ifstream &oriFile, std::vector< std::vector<double> > &posVec, std::vector< std::vector<double> > &oriVec);
void SendHelp();

int main(int argc, char *argv[]){

  //Command line arguments
  char *inFileName = NULL;
  char *outFileName = NULL;
  double wavelength = 435.;
  double opAng = 0.;
  bool isDiff = false;
  bool isCol = false;
  bool isData = true;//Set to true as it should always be data for now
                     //If I end up having to integrate this with SKDetSim
                     //we can alter this.
  std::string injector = " ";

  double pmtRadius;


  int opt;
  while ((opt = getopt(argc, argv, ":i:o:b:dcrvh")) != -1){
    switch(opt){
    case 'i':
      inFileName = optarg;
      break;
    case 'o':
      outFileName = optarg;
      break;
    case 'b':
      injector = optarg;
      break;
    case 'd':
      isDiff = true;
      opAng = 40.;
      break;
    case 'c':
      isCol = true;
      opAng = 4.;
      break;
    case 'r':
      isData = true;
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
    std::cout << ERR << "No input file found. Pass -h for help. Exiting.." << std::endl;
    return -1;
  }
  TFile *inFile = new TFile(inFileName, "READ");
  if(!inFile->IsOpen()){
    std::cout << ERR << "Cannot open input file: " << inFileName << std::endl;
    return -1;
  }

  //Check an injector flag has been passed, this is required
  if(injector == " "){
    std::cout << ERR << "No injector position passed. This is required. Pass -h for help. Exiting.." << std::endl;
    return -1;
  }

  //Check optic flags
  if(isDiff == false && isCol == false){
    std::cout << ERR << "No optic information passed. Pass -h for help. Exiting.." << std::endl;
    return -1;
  }
  else if(isDiff == true && isCol == true){
    std::cout << ERR << "Cannot run with both diffuser and collimator flags, only one should be passed. Exiting.." << std::endl;
    return -1;
  }

  //Open the PMT position/orientation tables and read those in
  std::vector< std::vector<double> > posVec, oriVec;
  if(isData){
    std::ifstream posFile, oriFile;
    posFile.open("./tables/pmt_position_cyl.dat");
    oriFile.open("./tables/pmt_orientation_cyl.dat");

    ReadPMTGeometry(posFile, oriFile, posVec, oriVec);

    posFile.close();
    oriFile.close();
  }
  
  //Open the output file
  if(outFileName == NULL) {
    std::cout << WAR << "No output filename provided. Falling back to default." << std::endl;
    outFileName = (char*)"output.root";
  }
  TFile *outFile = new TFile(outFileName, "RECREATE");
  std::cout << TAG << outFileName << " opened." << std::endl;
  //Finished flag checks and file inputs


  

  //Finally, access the event tree
  TTree *inTree = (TTree*)inFile->Get("tqtree");
  if(verbose) inTree->Print();

  Int_t year = 0;
  Int_t month = 0;
  Int_t day = 0;
  Int_t hour = 0;
  Int_t minute = 0;
  Int_t second = 0;
  Int_t run = 0;
  Int_t subrun = 0;
  Int_t nev = 0;
  std::vector<int> *ihit_vec = 0;
  std::vector<int> *cable_vec = 0;
  std::vector<float> *charge_vec = 0;
  std::vector<double> *time_vec = 0;
  std::vector<double> *pmtx_vec = 0;
  std::vector<double> *pmty_vec = 0;
  std::vector<double> *pmtz_vec = 0;

  inTree->SetBranchAddress("year", &year);
  inTree->SetBranchAddress("month", &month);
  inTree->SetBranchAddress("day", &day);
  inTree->SetBranchAddress("hour", &hour);
  inTree->SetBranchAddress("minute", &minute);
  inTree->SetBranchAddress("second", &second);
  inTree->SetBranchAddress("run", &run);
  inTree->SetBranchAddress("subrun", &subrun);
  inTree->SetBranchAddress("nev", &nev);
  inTree->SetBranchAddress("ihit_vec", &ihit_vec);
  inTree->SetBranchAddress("cable_vec", &cable_vec);
  inTree->SetBranchAddress("charge_vec", &charge_vec);
  inTree->SetBranchAddress("time_vec", &time_vec);
  inTree->SetBranchAddress("pmtx_vec", &pmtx_vec);
  inTree->SetBranchAddress("pmty_vec", &pmty_vec);
  inTree->SetBranchAddress("pmtz_vec", &pmtz_vec);

  //Ignore monitor info for now

  std::cout << TAG << "Tree loaded." << std::endl;


  //Save the PMT geometry
  double dist, xpos, ypos, zpos, costh, cosths, phis, omega, rad, eff;
  int pmtID, isPMTOn;
  //PMT tree - this stores positional information for ALL PMTs, not just those with hits
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



  //Injector position
  std::vector<double> vtx(3);
  //Get the position via the injector tag
  vtx = GetInjectorLocation(injector);

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
  //source_tree->Branch("nPhotons", &nPhotons);  //Initial simulated photons
  source_tree->Branch("opAng", &opAng);        //Injector opening angle

  source_tree->Fill();
  outFile->cd();
  source_tree->Write();

  
  //Set up injector local coordinate systems
  //Injection is inward perpendicular to wall
  //May switch these to TVector3s once I'm sure everything works correctly
  double vSourceDir[3];
  double vSourceXAxisLocal[3], vSourceYAxisLocal[3];
  double norm;
  double cut = 1795;//cut for top and bottom caps
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

  int nPMTs = posVec.size();

  if(verbose){
    printf("%s%i B&L PMTs.\n", TAG.c_str(), nPMTs);
  }

  pmtRadius = 25.; //Hardcode, probably need a more precise value
  
  double PMTpos[3], vOrient[3], vDir[3];
  double vecNorm, orientNorm, localx, localy;
  //Loop over pmts
  
  for (int i=0; i<nPMTs; i++){
    
    pmtID = i+1; //PMT IDs start from 1
    
    for(int j=0; j<3; j++){
      PMTpos[j]  = posVec.at(i).at(j);
      vOrient[j] = oriVec.at(i).at(j); //direction PMT faces
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
      
    rad = pmtRadius;
    omega  = GetSolidAngle(rad, dist);

    //is PMT working - currently all marked as good
    isPMTOn = 1;

    //PMT efficiency - currently all 100% efficient
    eff = 1.0;
    
    pmt_geom->Fill();
    
  }//End of loop over individual PMTs

  outFile->cd();
  pmt_geom->Write();

  
  //Now we can access the actual hit information
  //Will probably want to add additional info here,
  //but start with these for now.
  int nHits;
  double nPE, time_tofcorr, time;
  TTree *pmt_hits = new TTree("pmt_hits", "pmt_hits");
  pmt_hits->Branch("nHits", &nHits);
  pmt_hits->Branch("nPE", &nPE);
  pmt_hits->Branch("time", &time);
  pmt_hits->Branch("time_tofcorr", &time_tofcorr);
  pmt_hits->Branch("pmtID", &pmtID);
  
  std::cout << TAG << "Reading hit information." << std::endl;

  
  int pmtNum;
  double q, photonDist, tof;
  double photonDir[3];
  //Loop over events
  for(int ev=0; ev<inTree->GetEntries(); ev++){

    inTree->GetEvent(ev);
    nHits = ihit_vec->size();

    for (int ihit=0; ihit<nHits; ihit++){

      q = charge_vec->at(ihit);
      time = time_vec->at(ihit);
      pmtNum = cable_vec->at(ihit);
      pmtID = pmtNum; //Don't need the offset here
	
      for(int i=0; i<3; i++) 
	photonDir[i] = posVec.at(pmtID-1).at(i) - vtx[i];

      photonDist = sqrt(photonDir[0]*photonDir[0] + photonDir[1]*photonDir[1] + photonDir[2]*photonDir[2]);
      photonDist *= 1E-2; //Convert to m
      tof = photonDist/vg;
      time_tofcorr = time - tof;
	
      nPE = q;

      pmt_hits->Fill();

    }
    
  }

  outFile->cd();
  pmt_hits->Write();

  outFile->Close();

  std::cout << TAG << "Tree conversion complete." << std::endl;
  std::cout << "\u3042\u308A\u304C\u3068\u3046\u3054\u3056\u3044\u307E\u3057\u305F!" << std::endl;
  return 0;
}

//Put the real PMT geometry into arrays to be accessed in main
void ReadPMTGeometry(std::ifstream &posFile, std::ifstream &oriFile, std::vector< std::vector<double> > &posVec, std::vector< std::vector<double> > &oriVec){

  std::string line;
  std::vector<double> tmp(3);
  while(std::getline(posFile, line)){

    std::stringstream line_stream(line);
    std::string entry;
    std::vector<std::string> entries;
    char delim = ' ';

    while(std::getline(line_stream, entry, delim))
      entries.push_back(entry);

    //Add the x,y,z coordinates to the vector
    for(int i=0; i<3; i++)
      tmp[i] = std::stod(entries[i+1]);

    posVec.push_back(tmp);
    
  }//end of loop over file lines

  //Do the same thing for orientation
  while(std::getline(oriFile, line)){

    std::stringstream line_stream(line);
    std::string entry;
    std::vector<std::string> entries;
    char delim = ' ';

    while(std::getline(line_stream, entry, delim))
      entries.push_back(entry);

    //Add the x,y,z coordinates to the vector
    for(int i=0; i<3; i++)
      tmp[i] = std::stod(entries[i+1]);

    oriVec.push_back(tmp);
    
  }//end of loop over file lines

  //Sanity check
  if(posVec.size() != oriVec.size())
    std::cout << WAR << "PMT position and orientation vectors are different sizes! Something has gone wrong." << std::endl;

  return;

}


//This is a simple approximation
double GetSolidAngle(double r, double R){

  return TMath::Pi()*r*r/(R*R);

}


void SendHelp(){

  std::cout << TAG << "Script to read WCSim output and extract PMT hit info.\n"
	    << "USAGE: TreeConverter -i [input file] -b [injector] -c/d\n"
	    << "OPTIONS:\n"
	    << "-i : Input file\n"
	    << "-o : Output file\n"
	    << "-b : Injector [pass as B{number}]\n"
	    << "-c : Collimator\n"
	    << "-d : Diffuser\n"
	    << "-d : Run on real data\n"
	    << "-v : Print debug info\n"
	    << "-h : Display this message\n";

}
