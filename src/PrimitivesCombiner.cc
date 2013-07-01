#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "FWCore/Utilities/interface/Exception.h"


L1ITMu::PrimitiveCombiner::PrimitiveCombiner( const L1ITMu::PrimitiveCombiner::resolutions & res,
					      edm::ESHandle<DTGeometry> & muonGeom )
  : _resol(res), _muonGeom(muonGeom),
    _bx(0), _radialAngle(0), _bendingAngle(0), _bendingResol(0),
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
  if ( !isValid() ) return;

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
      _radialAngle = _radialAngle ? _radialAngle :_dtHI->getDTData().radialAngle;
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHI, _rpcOut ) );
      _radialAngle = _radialAngle ? _radialAngle :_dtHI->getDTData().radialAngle;
    }
  }


  if ( _dtHO ) {
    if ( _rpcIn ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcIn ) );
      _radialAngle = _radialAngle ? _radialAngle :_dtHO->getDTData().radialAngle;
    }
    if ( _rpcOut ) {
      localResults.push_back( combineDtRpc( _dtHO, _rpcOut ) );
      _radialAngle = _radialAngle ? _radialAngle :_dtHO->getDTData().radialAngle;
    }
  }

  double weightSum = 0;
  _bendingResol = 0;
  _bendingAngle = 0;

  std::vector<localResult>::const_iterator it = localResults.begin();
  std::vector<localResult>::const_iterator itend = localResults.end();
  for ( ; it != itend; ++it ) {
    weightSum += it->bendingResol;
    _bendingAngle += it->bendingAngle * it->bendingResol;
    _bendingResol += it->bendingResol * it->bendingResol;
  }

  _bendingAngle /= weightSum;
  _bendingResol = sqrt( _bendingResol );

}



L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDt( const L1ITMu::TriggerPrimitive * dt1,
				      const L1ITMu::TriggerPrimitive * dt2 )
{

  const DTChamber* chamb1 = _muonGeom->chamber( dt1->detId<DTChamberId>() );
  LocalPoint point1 = chamb1->toLocal( dt1->getCMSGlobalPoint() );

  const DTChamber* chamb2 = _muonGeom->chamber( dt2->detId<DTChamberId>() );
  LocalPoint point2 = chamb2->toLocal( dt2->getCMSGlobalPoint() );

  results localResult;
  localResult.bendingAngle = phiBCombined( point1.x(), point1.z(), point2.x(), point2.z() ) * 512;
  localResult.bendingResol = phiBCombinedResol( _resol.xDt, _resol.xDt );
  localResult.radialAngle = 0.5 * ( dt1->getDTData().radialAngle + dt2->getDTData().radialAngle );

  //  std::cout << "dt-dt radial : " << dt1->getDTData().radialAngle << " * " << dt2->getDTData().radialAngle << " = " << localResult.radialAngle << '\n';
//   std::cout << '\t' << point1.x() << '\t' << point1.z() << '\t' << dt1->getDTData().qualityCode << '\n';
//     std::cout << '\t' << point2.x() << '\t' << point2.z() << '\t' << dt2->getDTData().qualityCode << '\n';
//   std::cout << "dt-dt bending : " << dt1->getDTData().bendingAngle << " * " << dt2->getDTData().bendingAngle << " = "
// 	    << localResult.bendingAngle << '\n';

  return localResult;

}

L1ITMu::PrimitiveCombiner::results
L1ITMu::PrimitiveCombiner::combineDtRpc( const L1ITMu::TriggerPrimitive * dt,
					 const L1ITMu::TriggerPrimitive * rpc )
{

  const DTChamber* chamb1 = _muonGeom->chamber( dt->detId<DTChamberId>() );
  LocalPoint point1 = chamb1->toLocal( dt->getCMSGlobalPoint() );

  int station = rpc->detId<RPCDetId>().station();
  int sector  = rpc->detId<RPCDetId>().sector();
  int wheel = rpc->detId<RPCDetId>().ring();
  const DTChamber* chamb2 = _muonGeom->chamber( DTChamberId( wheel, station, sector ) );
  LocalPoint point2 = chamb2->toLocal( rpc->getCMSGlobalPoint() );

  results localResult;
  localResult.bendingAngle = phiBCombined( point1.x(), point1.z(), point2.x(), point2.z() ) * 512;
  localResult.bendingResol = phiBCombinedResol( _resol.xDt, _resol.xRpc );
  localResult.radialAngle = dt->getDTData().radialAngle;

  //std::cout << "dt-rpc bending : " << dt->getDTData().bendingAngle << " -> " << localResult.bendingAngle << '\n';

  return localResult;

}


// dt+rpc solo bending, phi dt
// dt+dt bending, media della posizione, direzione in base alla differenza dei due punti
// cancellare seconda traccia (bx diverso)
