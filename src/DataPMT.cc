#include "DataPMT.hh"

#include <iostream>

using namespace std;

ClassImp( DataPMT )

//////////////////////////////////////
//////////////////////////////////////

void DataPMT::ClearPMT(){

  // Set all the private member variables
  // to non-interpretive/physical values.

  SetPMTID( -1 );
  SetPMTType( -1 );
  SetPMTStatus( -1 );

  SetPMTPos( -99999.9, -99999.9, -99999.9 );
  SetPMTSourceDist( -99999.9 );
  SetPMTCosTh( -99999.9 );
  SetPMTPhi( -99999.9 );
  SetPMTPhotonAng( -99999.9 );
  SetPMTSolidAngle( -99999.9 );
  SetPMTRadius( -99999.9 );
  SetPMTEfficiency( -99999.9 );
  SetPMTnHits( -99999.9 );
  SetPMTnPE( -99999.9 );
  SetPMTTime( -99999.9 );
  SetPMTTime_ToFCorr( -99999.9 );

}
