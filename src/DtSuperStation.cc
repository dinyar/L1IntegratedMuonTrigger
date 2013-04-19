#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStation.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/Math/interface/deltaPhi.h"

L1ITMu::DtSuperStation::DtSuperStation( const DTChamberId & dtId )
{
  _wheel = dtId.wheel();
  _sector = dtId.sector();
  _station = dtId.station();
}

void L1ITMu::DtSuperStation::addStub(const TriggerPrimitiveRef& stub)
{

  TriggerPrimitive::subsystem_type type = stub->subsystem();

  switch ( type ) {
  case TriggerPrimitive::kDT :
    _dtAssociatedStubs.push_back( stub );
    break;
  case TriggerPrimitive::kRPC : {
    const RPCDetId & rpcId = stub->detId<RPCDetId>();

    if ( rpcId.region() ) { // endcap
      throw cms::Exception("Invalid Subsytem")
	<< "The RPC stub is not in a barrel layer" << std::endl;
    }

    if ( rpcId.layer() == 1 ) _rpcInAssociatedStubs.push_back( stub );
    else if ( rpcId.layer() == 2 ) _rpcOutAssociatedStubs.push_back( stub );
    else throw cms::Exception("Invalid Subsytem")
      << "The RPC layer is not a barrel layer" << std::endl;
    break;
  }
  default :
    throw cms::Exception("Invalid Subsytem") 
      << "The specified subsystem for this track stub is out of range"
      << std::endl;
  }
}



void L1ITMu::DtSuperStation::associate( double minRpcPhi )
{

  size_t dtSize = _dtAssociatedStubs.size();
  size_t rpcInSize = _rpcInAssociatedStubs.size();
  size_t rpcOutSize = _rpcOutAssociatedStubs.size();
  _dtMapAss.resize( dtSize );

  for ( size_t iDt = 0; iDt < dtSize; ++iDt ) {

    double phi = _dtAssociatedStubs.at(iDt)->getCMSGlobalPhi();
    std::map< double, size_t > rpcInIdx;
    std::map< double, size_t > rpcOutIdx;

    for ( size_t iIn = 0; iIn < rpcInSize; ++iIn ) {
      double phiIn = _rpcInAssociatedStubs.at( iIn )->getCMSGlobalPhi();
      double deltaPhiIn = fabs( reco::deltaPhi( phi, phiIn ) );
      if ( deltaPhiIn < minRpcPhi ) rpcInIdx[ deltaPhiIn ] = iIn;
    }

    for ( size_t iOut = 0; iOut < rpcOutSize; ++iOut ) {
      double phiOut = _rpcOutAssociatedStubs.at( iOut )->getCMSGlobalPhi();
      double deltaPhiOut = fabs( reco::deltaPhi( phi, phiOut ) );
      if ( deltaPhiOut < minRpcPhi ) rpcOutIdx[ deltaPhiOut ] = iOut;
    }


    DtSuperStation::primitiveAssociation & dtAss = _dtMapAss.at(iDt);

    /// fill up index for In associations
    std::map< double, size_t >::const_iterator it = rpcInIdx.begin();
    std::map< double, size_t >::const_iterator itend = rpcInIdx.end();
    dtAss.rpcIn.reserve( rpcInIdx.size() );
    for ( ; it != itend; ++it ) dtAss.rpcIn.push_back( it->second );

    /// fill up index for Out associations
    it = rpcOutIdx.begin();
    itend = rpcOutIdx.end();
    dtAss.rpcOut.reserve( rpcOutIdx.size() );
    for ( ; it != itend; ++it ) dtAss.rpcOut.push_back( it->second );

  }

}


L1ITMu::TriggerPrimitiveList L1ITMu::DtSuperStation::getRpcInAssociatedStubs( size_t dtIndex ) const
{

  L1ITMu::TriggerPrimitiveList returnList;

  try {
    const primitiveAssociation & prim = _dtMapAss.at(dtIndex);
    std::vector<size_t>::const_iterator it = prim.rpcIn.begin();
    std::vector<size_t>::const_iterator itend = prim.rpcIn.end();

    for ( ; it != itend; ++it ) returnList.push_back( _rpcInAssociatedStubs.at( *it ) );

  } catch ( const std::out_of_range & e ) {
    throw cms::Exception("DT Chamber Out of Range") 
      << "Requested DT primitive in position " << dtIndex << " out of " << _dtMapAss.size() << " total primitives"
      << std::endl;
  }

  return returnList;

}


L1ITMu::TriggerPrimitiveList L1ITMu::DtSuperStation::getRpcOutAssociatedStubs( size_t dtIndex ) const
{

  L1ITMu::TriggerPrimitiveList returnList;

  try {
    const primitiveAssociation & prim = _dtMapAss.at(dtIndex);

    std::vector<size_t>::const_iterator it = prim.rpcOut.begin();
    std::vector<size_t>::const_iterator itend = prim.rpcOut.end();

    for ( ; it != itend; ++it ) returnList.push_back( _rpcOutAssociatedStubs.at( *it ) );
  } catch ( const std::out_of_range & e ) {
    throw cms::Exception("DT Chamber Out of Range") 
      << "The number of dt primitives in sector are " << _dtMapAss.size()
      << std::endl;
  }

  return returnList;


}


