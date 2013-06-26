#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Utilities/interface/Exception.h"

L1ITMu::PrimitiveCombiner::PrimitiveCombiner( const DTChamberId & dtId, 
					      const L1ITMu::PrimitiveCombiner::resolutions & res )
  : _dtId(dtId), _resol(res), _bx(0),
    _radialAngle(0), _bendingAngle(0), _bendingResol(0), _isValid(0),
    _dtHI(0), _dtHO(0), _rpcIn(0), _rpcOut(0)
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


void L1ITMu::PrimitiveCombiner::addDtHI( const L1ITMu::TriggerPrimitive & prim )
{
  if ( _dtHI )
    throw cms::Exception("Invalid DT Quality") 
      << "DT primitive with quality HI already provided"
      << std::endl;

  _dtHI = &prim;
}


void L1ITMu::PrimitiveCombiner::addDtHO( const L1ITMu::TriggerPrimitive & prim )
{
  if ( _dtHO )
    throw cms::Exception("Invalid DT Quality") 
      << "DT primitive with quality HO already provided"
      << std::endl;

  _dtHO = &prim;
}


void L1ITMu::PrimitiveCombiner::addRpcIn( const L1ITMu::TriggerPrimitive & prim )
{
  _rpcIn = &prim;
}


void L1ITMu::PrimitiveCombiner::addRpcOut( const L1ITMu::TriggerPrimitive & prim )
{
  _rpcOut = &prim;
}



void L1ITMu::PrimitiveCombiner::combine()
{
  typedef L1ITMu::PrimitiveCombiner::results localResult;
  std::vector<localResult> localResults;

  _radialAngle = 0;
  if ( _dtHI && _dtHO ) {
    localResults.push_back( combineDt( _dtHI, _dtHO ) );
    _radialAngle = localResults.back().radialAngle;
  }

  if ( _dtHI ) {
    if ( _rpcIn ) {
      localResults.push_back( combineDtRpc( _dtHI, _rpcIn ) );
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHI, _rpcOut ) );
    }
    if ( ! _radialAngle ) _radialAngle = localResults.back().radialAngle;
  }


  if ( _dtHO ) {
    if ( _rpcIn ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcIn ) );
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcOut ) );
    }
    if ( ! _radialAngle ) _radialAngle = localResults.back().radialAngle;
  }

  double weightSum = 0;
  _bendingResol = 0;
  _bendingAngle = 0;

  std::vector<localResult>::const_iterator it = localResults.begin();
  std::vector<localResult>::const_iterator itend = localResults.end();
  for ( ; it != itend; ++it ) {
    weightSum += it->bendingResol;
    _bendingAngle += it->bendingAngle * it->bendingResol;
    _bendingResol += it->bendingResol;
  }

  _bendingAngle /= weightSum;
  _bendingResol /= weightSum;

  _isValid = true;

}



L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDt( const L1ITMu::TriggerPrimitive * dt1,
				      const L1ITMu::TriggerPrimitive * dt2 )
{

  GlobalPoint point1 = dt1->getCMSGlobalPoint();
  GlobalPoint point2 = dt2->getCMSGlobalPoint();

  results localResult;
  localResult.bendingAngle = phiBCombined( point1.x(), point2.x(), point1.z(), point2.z() );
  localResult.bendingResol = phiBCombinedResol( _resol.xDt, _resol.xDt );
  localResult.radialAngle = (dt1->getCMSGlobalPhi() + dt2->getCMSGlobalPhi()) * 0.5;
  return localResult;

}

L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDtRpc( const L1ITMu::TriggerPrimitive * dt,
					 const L1ITMu::TriggerPrimitive * rpc )
{

  GlobalPoint point1 = dt->getCMSGlobalPoint();
  GlobalPoint point2 = rpc->getCMSGlobalPoint();

  results localResult;
  localResult.bendingAngle = phiBCombined( point1.x(), point2.x(), point1.z(), point2.z() );
  localResult.bendingResol = phiBCombinedResol( _resol.xDt, _resol.xRpc );
  localResult.radialAngle = dt->getCMSGlobalPhi();
  return localResult;

}


// dt+rpc solo bending, phi dt
// dt+dt bending, media della posizione, direzione in base alla differenza dei due punti
// cancellare seconda traccia (bx diverso)
