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

#ifndef _DataRunInfo_
#define _DataRunInfo_

#include "DataPMT.hh"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <map>

using namespace std;

class DataRunInfo : public TObject{
  
public:

  // The constructors and destructor for the DataRunInfo object.
  DataRunInfo(){
    ClearRun();
    fDataPMTs.clear();
    fDatamPMTs.clear();
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

  // Fill PMT information for this run object.
  void FillPMTInfo( const string runFileName );

  // Fill mPMT information for this run object.
  void FillmPMTInfo( const string runFileName );

  // Get a PMT by PMT ID.
  DataPMT& GetPMT( const Int_t iPMT );

  // Get a PMT by PMT ID.
  DataPMT& GetmPMT( const Int_t imPMT );

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

  // Return the iterators corresponding to the beginning 
  // and end of the PMT map.
  map<Int_t, DataPMT>::iterator GetPMTIterBegin() { return fDataPMTs.begin(); }
  map<Int_t, DataPMT>::iterator GetPMTIterEnd() { return fDataPMTs.end(); }

  // Return the iterators corresponding to the beginning 
  // and end of the mPMT map.
  map<Int_t, DataPMT>::iterator GetmPMTIterBegin() { return fDatamPMTs.begin(); }
  map<Int_t, DataPMT>::iterator GetmPMTIterEnd() { return fDatamPMTs.end(); }
    
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
  map<Int_t, DataPMT> fDataPMTs; // Map of PMTs in this DataRun object
  map<Int_t, DataPMT> fDatamPMTs; // Map of mPMTs in this DataRun object

  ClassDef(DataRunInfo,1)
    
};
#endif
