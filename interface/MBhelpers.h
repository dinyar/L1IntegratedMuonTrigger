#ifndef __L1ITMU_MBHELPERS_H__
#define __L1ITMU_MBHELPERS_H__
// 
// Info: This is a collection of helpful free functions to making dealing
//       trigger primitives a bit more straight forward.
//
// Author: L. Gray (FNAL)
//

// #include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"
// #include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Ref.h"

namespace L1ITMu {
  namespace MBhelpers {

    MBLTVectorRef
      getPrimitivesByMBTriggerInfo(const int wheel,
                                   const int sp_wheel,
                                   const int sector,
                                   const int bx,
                                   const edm::Handle<MBLTContainer>& mbs,
                                   const unsigned mode,
                                   const std::vector<unsigned>& addrs);
  }
}

#endif
