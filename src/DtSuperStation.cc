#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStation.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"


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
