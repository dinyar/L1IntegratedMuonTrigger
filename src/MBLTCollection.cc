#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/Math/interface/deltaPhi.h"

L1ITMu::MBLTCollection::MBLTCollection( const DTChamberId & dtId )
{
  _wheel = dtId.wheel();
  _sector = dtId.sector();
  _station = dtId.station();
}

void L1ITMu::MBLTCollection::addStub(const TriggerPrimitiveRef& stub)
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



void L1ITMu::MBLTCollection::associate( double minRpcPhi )
{

  size_t dtSize = _dtAssociatedStubs.size();
  size_t rpcInSize = _rpcInAssociatedStubs.size();
  size_t rpcOutSize = _rpcOutAssociatedStubs.size();
  _dtMapAss.resize( dtSize );


//   std::vector< std::map<double, size_t> > dtIdxIn;
//   dtIdxIn.resize(rpcInSize);
//   std::vector< std::map<double, size_t> > dtIdxOut;
//   dtIdxOut.resize(rpcOutSize);

  for ( size_t iDt = 0; iDt < dtSize; ++iDt ) {

    double phi = _dtAssociatedStubs.at(iDt)->getCMSGlobalPhi();
    std::map< double, size_t > rpcInIdx;
    std::map< double, size_t > rpcOutIdx;

    for ( size_t iIn = 0; iIn < rpcInSize; ++iIn ) {
      double phiIn = _rpcInAssociatedStubs.at( iIn )->getCMSGlobalPhi();
      double deltaPhiIn = fabs( reco::deltaPhi( phi, phiIn ) );
      if ( deltaPhiIn < minRpcPhi ) {
	rpcInIdx[ deltaPhiIn ] = iIn;
// 	dtIdxIn[iIn][ deltaPhiIn ] = iDt;
      }
    }

    for ( size_t iOut = 0; iOut < rpcOutSize; ++iOut ) {
      double phiOut = _rpcOutAssociatedStubs.at( iOut )->getCMSGlobalPhi();
      double deltaPhiOut = fabs( reco::deltaPhi( phi, phiOut ) );
      if ( deltaPhiOut < minRpcPhi ) {
	rpcOutIdx[ deltaPhiOut ] = iOut;
// 	dtIdxOut[iOut][ deltaPhiOut ] = iDt;
      }
    }


    MBLTCollection::primitiveAssociation & dtAss = _dtMapAss.at(iDt);

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


L1ITMu::TriggerPrimitiveList L1ITMu::MBLTCollection::getRpcInAssociatedStubs( size_t dtIndex ) const
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


L1ITMu::TriggerPrimitiveList L1ITMu::MBLTCollection::getRpcOutAssociatedStubs( size_t dtIndex ) const
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






L1ITMu::MBLTCollection::bxMatch L1ITMu::MBLTCollection::haveCommonRpc( size_t dt1, size_t dt2 ) const
{
  
  L1ITMu::MBLTCollection::bxMatch ret_val = NOMATCH;

  if ( dt1 == dt2 ) {
    throw cms::Exception("DT primitive compared to itself") 
      << "The two id passed refer to the same primitive"
      << std::endl;
  }

  try {
    const primitiveAssociation & prim1 = _dtMapAss.at(dt1);
    const primitiveAssociation & prim2 = _dtMapAss.at(dt2);

//     bool in_match = false;
//     bool out_match = false;

//     size_t rpcInSize1 = prim1.rpcIn.size();
//     size_t rpcInSize2 = prim2.rpcIn.size();
//     for ( size_t i = 0; i < rpcInSize1; ++i ) 
//       for ( size_t j = 0; j < rpcInSize2; ++j )
// 	if ( prim1.rpcIn[i] == prim1.rpcIn[j] ) {
// 	  in_match = true;
// 	  i = rpcInSize1;
// 	  break;
// 	}

//     size_t rpcOutSize1 = prim1.rpcOut.size();
//     size_t rpcOutSize2 = prim2.rpcOut.size();
//     for ( size_t i = 0; i < rpcOutSize1; ++i ) 
//       for ( size_t j = 0; j < rpcOutSize2; ++j )
// 	if ( prim1.rpcOut[i] == prim1.rpcOut[j] ) {
// 	  out_match = true;
// 	  i = rpcOutSize1;
// 	  break;
// 	}
//     if ( in_match && out_match ) return FULLMATCH;
//     else if ( in_match ) return INMATCH;
//     else if ( out_match ) return OUTMATCH;
//     return NOMATCH;


    if ( !prim1.rpcIn.empty() && !prim2.rpcIn.empty() ) {
      if ( prim1.rpcIn.front() == prim2.rpcIn.front() ) {
    	ret_val = INMATCH;
      }
    }

    if ( !prim1.rpcOut.empty() && !prim2.rpcOut.empty() ) {
      if ( prim1.rpcOut.front() == prim2.rpcOut.front() ) {
    	ret_val = ( ret_val == INMATCH ) ? FULLMATCH : OUTMATCH;
      }
    }
    return ret_val;

  } catch ( const std::out_of_range & e ) {
    throw cms::Exception("DT Chamber Out of Range") 
      << "The number of dt primitives in sector are " << _dtMapAss.size()
      << std::endl;
  }

  return ret_val;

}


