////////////////////////////////////////////////////////////////////
///
/// FILENAME: DataRunInfo.hh
///
/// CLASS: DataRunInfo
///
/// BRIEF: Run source level data structure for analysis
///          
/// AUTHOR: Ana Sofia Inacio <ana.carpinteiroinacio@physics.ox.ac.uk>
///
/// DETAIL: This data structure contains all the information
///         for a specific calibration run required as an input
///         into an optics fit/analysis.
///
///         Herein lies information specific to the run itself,
///         such as the injector position, the wavelength, photon
///         group velocity, initial number of photons and injector
///         opening angle. 
///
////////////////////////////////////////////////////////////////////

// NOTE by A.S.I: Will potentially want to include a map of PMT objects that correspond to this run. Look at the way it is done in the SNO+ OCA code

#ifndef _DataRunInfo_
#define _DataRunInfo_

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace std;

class DataRunInfo{
  
public:

  // The constructors and destructor for the DataRunInfo object.
  DataRunInfo(){
    ClearRun();
    //fRunPMTs.clear();
  }
  ~DataRunInfo(){ }

  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////
    
  // Initialise and clear all the private member variables
  // to non-physical/interpretive values before assigning them.
  void ClearRun();

  // Fill Run information from a TreeConverter file
  // into the run information stored here.
  void FillRunInfo( const string runFileName );

  /////////////////////////////////
  ////////     GETTERS     ////////
  /////////////////////////////////

  // Get the source position (mm).
  TVector3 GetSourcePos() const { return fSourcePos; } 

  // Get the wavelength of the laser (nm).
  Double_t GetLambda() const { return fLambda; }

  // Get the group velocity of the photons (nm s-1).
  Double_t GetGroupVel() const { return fGroupVel; } 
   
  // Get the number of photons.
  Double_t GetNPhotons() const { return fNPhotons; }

  // Get the injector opening angle.
  Double_t GetOpAng() const { return fOpAng; }
    
  /////////////////////////////////
  ////////     SETTERS     ////////
  /////////////////////////////////

  // Set the source position (mm).
  void SetSourcePos( const Double_t xPos, const Double_t yPos, const Double_t zPos ){ 
    fSourcePos.SetX( xPos );
    fSourcePos.SetY( yPos );
    fSourcePos.SetZ( zPos );
  }
  void SetSourcePos( const TVector3 xyzPos ){ 
    fSourcePos = xyzPos;
  }
    
  // Set the wavelength of the laser (nm).
  void SetLambda( const Double_t lambda ){ fLambda = lambda; }

  // Set the group velocity of the photons (nm s-1).
  void SetGroupVel( const Double_t vg ){ fGroupVel = vg; } 

  // Set the number of photons.
  void SetNPhotons( const Double_t nphotons ){ fNPhotons = nphotons; }
  
  // Set the injector opening angle.
  void SetOpAng( const Double_t opAng ){ fOpAng = opAng; } 
    
private:
    
  TVector3 fSourcePos;        // The source position
  Double_t fLambda;           // The wavelength of the laser
  Double_t fGroupVel;         // The photon group velocity
  Double_t fNPhotons;         // The number of initial photons
  Double_t fOpAng;            // The source opening angle
  
  //RunPMT Object

  ClassDef(DataRunInfo,1)
    
};
#endif
