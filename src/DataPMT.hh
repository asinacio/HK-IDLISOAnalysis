////////////////////////////////////////////////////////////////////
///
/// FILENAME: DataPMT.hh
///
/// CLASS: DataPMT
///
/// BRIEF: PMT level data structure for the optics analysis
///          
/// AUTHOR: Ana Sofia Inacio <ana.carpinteiroinacio@physics.ox.ac.uk>
///
/// DETAIL: This data structure contains all the information
///         for a specific PMT in a run, required as an input
///         to a DataPMT object for use in an optics 
///         fit/analysis.
///
///         Herein lies information specific to the PMT itself,
///         such as distance from the source, PMT ID, solid angle etc.
///
////////////////////////////////////////////////////////////////////

#ifndef _DataPMT_
#define _DataPMT_

#include "TVector3.h"
#include "TObject.h"

class DataPMT : public TObject{
public:

  // The constructors and destructors for the DataPMT object.
  DataPMT(){ ClearPMT(); }
  DataPMT( Int_t pmtID ){ 
    fPMTID = pmtID; 
    ClearPMT(); 
  }
  ~DataPMT(){ }

  /////////////////////////////////
  ////////     METHODS     ////////
  /////////////////////////////////
    
  // Initialise and clear all the private member variables
  // to non-physical values before assigning them.
  void ClearPMT();

  /////////////////////////////////
  ////////     GETTERS     ////////
  /////////////////////////////////
    
  // Get the ID for the PMT.
  Int_t GetPMTID() const { return fPMTID; }

  // Set the PMT type (1 = normal, 2 = mPMT).
  Int_t GetPMTType() const { return fPMTType; }

  // Set the PMT status (0 = off, 1 = on).
  Int_t GetPMTStatus() const { return fPMTOn; }

  // Set the PMT position.
  TVector3 GetPMTPos() const { return fPMTPos; }

  // Set the distance from source to PMT.
  Double_t GetPMTSourceDist() const { return fPMTSourceDist; }

  // Set the PMT angle theta relative to the source.
  Double_t GetPMTCosTh() const { return fPMTCosTh; }

  // Set the PMT angle phi relative to the source.
  Double_t GetPMTPhi() const { return fPMTPhi; }

  // Set the photon angle relative to the PMT normal.
  Double_t GetPMTPhotonAng() const { return fPMTPhotonAng; }

  // Set the PMT subtended solid angle.
  Double_t GetPMTSolidAngle() const { return fPMTSolidAngle; }

  // Set the PMT radius.
  Double_t GetPMTRadius() const { return fPMTRadius; }

  // Set the PMT efficiency.
  Double_t GetPMTEfficiency() const { return fPMTEff; }

  // Set the nHits recorded by the PMT.
  Double_t GetPMTnHits() const { return fPMTnHits; }

  // Set the number of photoelectrons (PE) recorded by the PMT.
  Double_t GetPMTnPE() const { return fPMTnPE; }

  // Set the time recorded by the PMT.
  Double_t GetPMTTime() const { return fPMTTime; }

  // Set the time recorded by the PMT, corrected by the time-of-flight.
  Double_t GetPMTTime_ToFCorr() const { return fPMTTime_ToFCorr; }

  /////////////////////////////////
  ////////     SETTERS     ////////
  /////////////////////////////////
    
  // Set the ID for the PMT.
  void SetPMTID( const Int_t pmtID ){ fPMTID = pmtID; }

  // Set the PMT type (1 = normal, 2 = mPMT).
  void SetPMTType( const Int_t pmtType ){ fPMTType = pmtType; }

  // Set the PMT status (0 = off, 1 = on).
  void SetPMTStatus( const Int_t pmtStatus ){ fPMTOn = pmtStatus; }

  // Set the PMT position.
  void SetPMTPos( const Double_t xPos, const Double_t yPos, const Double_t zPos ){ 
    fPMTPos.SetX( xPos );
    fPMTPos.SetY( yPos );
    fPMTPos.SetZ( zPos );
  }
  void SetPMTPos( const TVector3 xyzPos ){ 
    fPMTPos = xyzPos;
  }

  // Set the distance from source to PMT.
  void SetPMTSourceDist( const Double_t dist ){ fPMTSourceDist = dist; }

  // Set the PMT angle theta relative to the source.
  void SetPMTCosTh( const Double_t costh ){ fPMTCosTh = costh; }

  // Set the PMT angle phi relative to the source.
  void SetPMTPhi( const Double_t phi ){ fPMTPhi = phi; }

  // Set the photon angle relative to the PMT normal.
  void SetPMTPhotonAng( const Double_t photonAng ){ fPMTPhotonAng = photonAng; }

  // Set the PMT subtended solid angle.
  void SetPMTSolidAngle( const Double_t solidAng ){ fPMTSolidAngle = solidAng; }

  // Set the PMT radius.
  void SetPMTRadius( const Double_t radius ){ fPMTRadius = radius; }

  // Set the PMT efficiency.
  void SetPMTEfficiency( const Double_t eff ){ fPMTEff = eff; }

  // Set the nHits recorded by the PMT.
  void SetPMTnHits( const Double_t nhit ){ fPMTnHits = nhit; }

  // Set the number of photoelectrons (PE) recorded by the PMT.
  void SetPMTnPE( const Double_t npe ){ fPMTnPE = npe; }

  // Set the time recorded by the PMT.
  void SetPMTTime( const Double_t time ){ fPMTTime = time; }

  // Set the time recorded by the PMT, corrected by the time-of-flight.
  void SetPMTTime_ToFCorr( const Double_t time_corr ){ fPMTTime_ToFCorr = time_corr; }

private:
    
  Int_t fPMTID;                // PMT ID
  Int_t fPMTType;              // The PMT type (1 = normal, 2 = mPMT)
  Int_t fPMTOn;                // PMT status (0 = off, 1 = on) - can add extra status tags if needed
  TVector3 fPMTPos;            // PMT Position
  Double_t fPMTSourceDist;     // Distance to PMT from source
  Double_t fPMTCosTh;          // PMT theta relative to source
  Double_t fPMTPhi;            // PMT phi relative to source
  Double_t fPMTPhotonAng;      // Photon angle relative to PMT
  Double_t fPMTSolidAngle;     // Solid angle subtended
  Double_t fPMTRadius;         // PMT Radius
  Double_t fPMTEff;            // PMT efficiency
  Double_t fPMTnHits;          // nHits recorded by PMT
  Double_t fPMTnPE;            // Number of PE recorded by PMT
  Double_t fPMTTime;           // Time recorded by PMT
  Double_t fPMTTime_ToFCorr;   // Time recorded by PMT, corrected by the time-of-flight

  ClassDef( DataPMT, 1 );
    
};
#endif
