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

  if ( !fDataPMTs.empty() ){ fDataPMTs.clear(); }
  if ( !fDatamPMTs.empty() ){ fDatamPMTs.clear(); }
  
}

//////////////////////////////////////
//////////////////////////////////////

void DataRunInfo::FillRunInfo( const string runFileName ){

  // Open the input file and retrieve run info from tree
  TFile *inFile = new TFile( runFileName.c_str() );

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
  double sourcePosx, sourcePosy, sourcePosz, sourceWl, sourceI, sourceAng, groupVel;
  tsource->SetBranchAddress("vtx_x", &sourcePosx);  // Source x coordinate
  tsource->SetBranchAddress("vtx_y", &sourcePosy);  // Source y coordinate
  tsource->SetBranchAddress("vtx_z", &sourcePosz);  // Source z coordinate
  tsource->SetBranchAddress("lambda", &sourceWl);   // Source wavelength (nm)
  tsource->SetBranchAddress("nPhotons", &sourceI);  // Source intensity, number of photons
  tsource->SetBranchAddress("opAng", &sourceAng);   // Source opening angle
  tsource->SetBranchAddress("vg", &groupVel);       // Photon group velocity (nm s-1)

  tsource->GetEntry(0);

  SetSourcePos( sourcePosx, sourcePosy, sourcePosz );
  SetLambda( sourceWl );
  SetGroupVel( groupVel );
  SetNPhotons( sourceI );
  SetOpAng( sourceAng );

  inFile->Close();
}

//////////////////////////////////////
//////////////////////////////////////

void DataRunInfo::FillPMTInfo( const string runFileName ){

  // Open the input file and retrieve run info from tree
  TFile *inFile = new TFile( runFileName.c_str(), "read" );

  // If files does not exist, print a message and exit
  if( !inFile || inFile->IsZombie() ){
    std::cerr << "Error opening file!" << endl;
    exit(-1);
  }

  cout << "Reading PMT info from file " << runFileName << "..." << endl;

  // Loading the ttrees
  TTree *tpmt_geom = (TTree*)inFile->Get("pmt_geom");
  TTree *tpmt_hits = (TTree*)inFile->Get("pmt_hits");

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

  for( Int_t i = 0; i < tpmt_geom->GetEntries(); i++ ){

    tpmt_geom->GetEntry(i);

    fDataPMTs[pmtID].SetPMTID( pmtID );
    fDataPMTs[pmtID].SetPMTType( 1 );
    fDataPMTs[pmtID].SetPMTStatus( isPMTOn );
    fDataPMTs[pmtID].SetPMTPos( xpos, ypos, zpos );
    fDataPMTs[pmtID].SetPMTSourceDist( dist );
    fDataPMTs[pmtID].SetPMTCosTh( cosths );
    fDataPMTs[pmtID].SetPMTPhi( phis );
    fDataPMTs[pmtID].SetPMTPhotonAng( costh );
    fDataPMTs[pmtID].SetPMTSolidAngle( omega );
    fDataPMTs[pmtID].SetPMTRadius( rad );
    fDataPMTs[pmtID].SetPMTEfficiency( eff );

  }

  double nHits, nPE, time_tofcorr, time;
  tpmt_hits->SetBranchAddress("nHits", &nHits);
  tpmt_hits->SetBranchAddress("nPE", &nPE);
  tpmt_hits->SetBranchAddress("time", &time);
  tpmt_hits->SetBranchAddress("time_tofcorr", &time_tofcorr);
  tpmt_hits->SetBranchAddress("pmtID", &pmtID);

  for( Int_t i = 0; i < tpmt_hits->GetEntries(); i++ ){

    tpmt_hits->GetEntry(i);

    fDataPMTs[pmtID].SetPMTnHits( nHits );
    fDataPMTs[pmtID].SetPMTnPE( nPE );
    fDataPMTs[pmtID].SetPMTTime( time );
    fDataPMTs[pmtID].SetPMTTime_ToFCorr( time_tofcorr );
    
  }

  inFile->Close();

}

//////////////////////////////////////
//////////////////////////////////////

void DataRunInfo::FillmPMTInfo( const string runFileName ){

  // Open the input file and retrieve run info from tree
  TFile *inFile = new TFile( runFileName.c_str(), "read" );

  // If files does not exist, print a message and exit
  if( !inFile || inFile->IsZombie() ){
    std::cerr << "Error opening file!" << endl;
    exit(-1);
  }

  cout << "Reading mPMT info from file " << runFileName << "..." << endl;

  // Loading the ttrees
  TTree *tmpmt_geom = (TTree*)inFile->Get("mpmt_geom");
  TTree *tmpmt_hits = (TTree*)inFile->Get("mpmt_hits");

  // Loading PMT information
  cout << "Loading mPMT data..." << endl;
  double dist, xpos, ypos, zpos, costh, cosths, phis, omega, rad, eff;
  int pmtID, mpmtID, isPMTOn;
  tmpmt_geom->SetBranchAddress("R", &dist);          // Distance to PMT from source
  tmpmt_geom->SetBranchAddress("x", &xpos);          // PMT x coordinate
  tmpmt_geom->SetBranchAddress("y", &ypos);          // PMT y coordinate
  tmpmt_geom->SetBranchAddress("z", &zpos);          // PMT z coordinate
  tmpmt_geom->SetBranchAddress("costh", &costh);     // Photon angle relative to PMT
  tmpmt_geom->SetBranchAddress("cosths", &cosths);   // PMT theta relative to source
  tmpmt_geom->SetBranchAddress("phis", &phis);       // PMT phi relative to source
  tmpmt_geom->SetBranchAddress("omega", &omega);     // Solid angle subtended
  tmpmt_geom->SetBranchAddress("radius", &rad);      // PMT radius
  tmpmt_geom->SetBranchAddress("eff", &eff);         // PMT efficiency
  tmpmt_geom->SetBranchAddress("pmtID", &pmtID);     // PMT ID (enum. from 0)
  tmpmt_geom->SetBranchAddress("isPMTOn", &isPMTOn); // PMT status (0 = off, 1 = on) - can add extra status tags if needed


  for( Int_t i = 0; i < tmpmt_geom->GetEntries(); i++ ){

    tmpmt_geom->GetEntry(i);

    fDatamPMTs[pmtID].SetPMTID( pmtID );
    fDatamPMTs[pmtID].SetPMTType( 2 );
    fDatamPMTs[pmtID].SetPMTStatus( isPMTOn );
    fDatamPMTs[pmtID].SetPMTPos( xpos, ypos, zpos );
    fDatamPMTs[pmtID].SetPMTSourceDist( dist );
    fDatamPMTs[pmtID].SetPMTCosTh( cosths );
    fDatamPMTs[pmtID].SetPMTPhi( phis );
    fDatamPMTs[pmtID].SetPMTPhotonAng( costh );
    fDatamPMTs[pmtID].SetPMTSolidAngle( omega );
    fDatamPMTs[pmtID].SetPMTRadius( rad );
    fDatamPMTs[pmtID].SetPMTEfficiency( eff );

  }
  
  double nHits, nPE, time_tofcorr, time;
  tmpmt_hits->SetBranchAddress("nHits", &nHits);
  tmpmt_hits->SetBranchAddress("nPE", &nPE);
  tmpmt_hits->SetBranchAddress("time", &time);
  tmpmt_hits->SetBranchAddress("time_tofcorr", &time_tofcorr);
  tmpmt_hits->SetBranchAddress("pmtID", &pmtID);

  for( Int_t i = 0; i < tmpmt_hits->GetEntries(); i++ ){

    tmpmt_hits->GetEntry(i);

    fDatamPMTs[pmtID].SetPMTnHits( nHits );
    fDatamPMTs[pmtID].SetPMTnPE( nPE );
    fDatamPMTs[pmtID].SetPMTTime( time );
    fDatamPMTs[pmtID].SetPMTTime_ToFCorr( time_tofcorr );
    
  }

  inFile->Close();

}

//////////////////////////////////////
//////////////////////////////////////

DataPMT& DataRunInfo::GetPMT( Int_t iPMT )
{

  // Return the PMT with ID 'iPMT'
  return fDataPMTs[ iPMT ];
  
}

//////////////////////////////////////
//////////////////////////////////////

DataPMT& DataRunInfo::GetmPMT( Int_t imPMT )
{

  // Return the mPMT with ID 'imPMT'
  return fDatamPMTs[ imPMT ];
  
}

