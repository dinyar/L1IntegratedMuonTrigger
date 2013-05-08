#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "FWCore/Utilities/interface/Exception.h"

L1ITMu::PrimitiveCombiner::PrimitiveCombiner( const DTChamberId & dtId, 
					      const L1ITMu::PrimitiveCombiner::resolutions & res )
  : _dtId(dtId), _resol(res), _bx(0),
    _radialAngle(0), _bendingAngle(0), _bendingResol(0), _isValid(0)
{}



void L1ITMu::PrimitiveCombiner::addDt( const L1ITMu::TriggerPrimitive & dt )
{
  int qualityCode = dt.getDTData().qualityCode;
  switch ( qualityCode ) {
  case 2 : addDtHO( dt ); break;
  case 3 : addDtHI( dt ); break;
  default :
    throw cms::Exception("Invalid DT Quality") 
      << "This module can combine only HI to HO (quality2/3), provided : "
      << qualityCode << std::endl;
  }

}


/// FIXME: check new exposed variables and methods from L1ITMu::TriggerPrimitive
void L1ITMu::PrimitiveCombiner::addDtHI( const L1ITMu::TriggerPrimitive & prim )
{
  if ( _dtHI )
    throw cms::Exception("Invalid DT Quality") 
      << "DT primitive with quality HI already provided"
      << std::endl;

  float x = 0; // xPos( dtId.wheel(), dtId.sector(), dtId.station() );
  float z = 0; // zPos( dtId.wheel(), dtId.sector(), prim.getCMSGlobalPhi() );
  _dtHI = L1ITMu::PrimitiveCombiner::primitive( &prim, x, z );
}
/// FIXME END


/// FIXME: check new exposed variables and methods from L1ITMu::TriggerPrimitive
void L1ITMu::PrimitiveCombiner::addDtHO( const L1ITMu::TriggerPrimitive & prim )
{
  if ( _dtHO )
    throw cms::Exception("Invalid DT Quality") 
      << "DT primitive with quality HO already provided"
      << std::endl;

  float x = 0; // xPos( dtId.wheel(), dtId.sector(), dtId.station() );
  float z = 0; // zPos( dtId.wheel(), dtId.sector(), prim.getCMSGlobalPhi() );
  _dtHO = L1ITMu::PrimitiveCombiner::primitive( &prim, x, z );
}
/// FIXME END


/// FIXME: check new exposed variables and methods from L1ITMu::TriggerPrimitive
void L1ITMu::PrimitiveCombiner::addRpcIn( const L1ITMu::TriggerPrimitive & prim )
{
  float x = 0; // xPos( dtId.wheel(), dtId.sector(), dtId.station() );
  float z = 0; // zPos( dtId.wheel(), dtId.sector(), prim.getCMSGlobalPhi() );
  _rpcIn = L1ITMu::PrimitiveCombiner::primitive( &prim, x, z );
}
/// FIXME END


/// FIXME: check new exposed variables and methods from L1ITMu::TriggerPrimitive
void L1ITMu::PrimitiveCombiner::addRpcOut( const L1ITMu::TriggerPrimitive & prim )
{
  float x = 0; // xPos( dtId.wheel(), dtId.sector(), dtId.station() );
  float z = 0; // zPos( dtId.wheel(), dtId.sector(), prim.getCMSGlobalPhi() );
  _rpcOut = L1ITMu::PrimitiveCombiner::primitive( &prim, x, z );
}
/// FIXME END




void L1ITMu::PrimitiveCombiner::combine()
{
  typedef L1ITMu::PrimitiveCombiner::results localResult;
  std::vector<localResult> localResults;

  if ( _dtHI && _dtHO ) {
    localResults.push_back( combineDt( _dtHI, _dtHO ) );
  }

  if ( _dtHI ) {
    if ( _rpcIn ) {
      localResults.push_back( combineDtRpc( _dtHI, _rpcIn ) );
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHI, _rpcOut ) );
    }
  }


  if ( _dtHO ) {
    if ( _rpcIn ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcIn ) );
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcOut ) );
    }
  }

  double weightSum = 0;
  _bendingResol = 0;
  _bendingAngle = 0;

  std::vector<localResult>::const_iterator it = localResults.begin();
  std::vector<localResult>::const_iterator itend = localResults.end();
  for ( ; it != itend; ++it ) {
    weightSum += it->weight;
    _bendingAngle += it->phiBComb * it->weight;
    _bendingResol += it->phiBCombResol * it->weight;
  }

  _bendingAngle /= weightSum;
  _bendingResol /= weightSum;

  /// FIXME: what goes here?
  _isValid = true;
  _bendingAngle = 0;
  /// FIXME END

}



L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDt( const primitive & dt1,
				      const primitive & dt2 )
{

  results localResult;

  localResult.phiBComb = phiBCombined( dt1.x, dt2.x, dt1.z, dt2.z );
  localResult.phiBCombResol = phiBCombinedResol( _resol.xDt, _resol.xDt );
  localResult.radialAngle = 0;
  localResult.bendingAngle = 0;
  localResult.weight = 0;

  return localResult;

}

L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDtRpc( const primitive & dt,
					 const primitive & rpc )
{

  results localResult;

  localResult.phiBComb = phiBCombined( dt.x, rpc.x, dt.z, rpc.z );
  localResult.phiBCombResol = phiBCombinedResol( _resol.xDt, _resol.xRpc );
  localResult.radialAngle = 0;
  localResult.bendingAngle = 0;
  localResult.weight = 0;

  return localResult;

}
