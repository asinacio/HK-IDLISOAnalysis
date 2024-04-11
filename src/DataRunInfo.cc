#include "DataRunInfo.hh"

#include <iostream>

using namespace std;

ClassImp( DataRunInfo )

//////////////////////////////////////
//////////////////////////////////////

void DataRunInfo::ClearRun(){ 

  // Set all the private member variables
  // to non-interpretive/physical values.

  SetSourcePos( -9999.9, -9999.9, -9999.9 );
  SetLambda( -10.0 );
  SetGroupVel( -10.0 );
  SetNPhotons( -10.0 );
  SetOpAng( -10.0 );
  
}

//////////////////////////////////////
//////////////////////////////////////

void DataRunInfo::FillRunInfo( const string runFileName ){

  // Open the input file and retrieve run info from tree
  TFile *inFile = new TFile( runFileName.c_str(), "read" );

  // If files does not exist, print a message and exit
  if( !inFile || inFile->IsZombie() ){
    std::cerr << "Error opening file!" << endl;
    exit(-1);
  }

  cout << "Reading run info from file " << runFileName << "..." << endl;

  // Loading the ttrees
  TTree *tsource = (TTree*)inFile->Get("source");

  // Loading information about the source and laser
  cout << "Loading calibration source information..." << endl;
  Double_t sourcePosx, sourcePosy, sourcePosz, sourceWl, sourceI, sourceAng, groupVel;
  tsource->SetBranchAddress("vtx_x", &sourcePosx);  // Source x coordinate
  tsource->SetBranchAddress("vtx_y", &sourcePosy);  // Source y coordinate
  tsource->SetBranchAddress("vtx_z", &sourcePosz);  // Source z coordinate
  tsource->SetBranchAddress("lambda", &sourceWl);   // Source wavelength (nm)
  tsource->SetBranchAddress("nPhotons", &sourceI);  // Source intensity, number of photons
  tsource->SetBranchAddress("opAng", &sourceAng);   // Source opening angle
  tsource->SetBranchAddress("vg", &groupVel);       // Photon group velocity (nm s-1)

  SetSourcePos( sourcePosx, sourcePosy, sourcePosz );
  SetLambda( sourceWl );
  SetGroupVel( groupVel );
  SetNPhotons( sourceI );
  SetOpAng( sourceAng );

  inFile->Close();
}
