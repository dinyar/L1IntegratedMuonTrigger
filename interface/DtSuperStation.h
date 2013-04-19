#ifndef L1Trigger_L1IntegratedMuonTrigger_DtSuperStation_h_
#define L1Trigger_L1IntegratedMuonTrigger_DtSuperStation_h_
// 
// Class: L1ITMu:: 
//
// Info: This track represents a DT(1 station) plus eventual RPC station(s)
//       based track seed, from which a full multi-station track can be
//       built.
//
// Author: 
//

#include <iostream>

#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"


namespace L1ITMu {
  
  class DtSuperStation {

  public:

    /// internal enum for subdetector stub identification
    enum subsystem_offset{ kDT, kRPCb, kCSC, kRPCf };

    /// structure for internal indexing
    struct primitiveAssociation {
      std::vector< size_t > rpcIn;
      std::vector< size_t > rpcOut;
    };
    
    /// default constructor
    DtSuperStation() :_wheel(0),_sector(0),_station(0) {};

    /// construction out of DTChamberId: automatically extracts info
    DtSuperStation( const DTChamberId & dtId );
    ~DtSuperStation() {};

    /// selectively add Trigger Primitives to the SuperSegment
    /// dt, rpc up layer and rpc down layer are stored in separated collections
    void addStub( const TriggerPrimitiveRef& stub );


    /// return a reference to the DT only segments
    const TriggerPrimitiveList & getDtSegments() const {
      return _dtAssociatedStubs;
    }

    /// rpc inner layer hits only
    const TriggerPrimitiveList & getRpcInner() const {
      return _rpcInAssociatedStubs;
    }

    /// rpc outer layer hits only
    const TriggerPrimitiveList & getRpcOuter() const {
      return _rpcOutAssociatedStubs;
    }

    /// returns wheel
    inline int wheel() const { return _wheel; }

    /// returns sector
    inline int sector() const { return _sector; }

    /// returns station
    inline int station() const { return _station; }

    /// rpc inner layer hits associated to a given dt station
    TriggerPrimitiveList getRpcInAssociatedStubs( size_t dtIndex ) const;

    /// rpc outer layer hits associated to a given dt station
    TriggerPrimitiveList getRpcOutAssociatedStubs( size_t dtIndex ) const;

    // build association map among dt and rpc primitives
    void associate( double );

  private :

    /// dt segments
    TriggerPrimitiveList _dtAssociatedStubs;

    /// rpc inner layer hits
    TriggerPrimitiveList _rpcInAssociatedStubs;

    /// rpc outer layer hits
    TriggerPrimitiveList _rpcOutAssociatedStubs;

    /// space coordinates
    int _wheel, _sector, _station;

    // association map among dt and rpc primitives
    std::vector< primitiveAssociation > _dtMapAss;

  };
}

#endif
