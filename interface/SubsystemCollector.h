#ifndef __L1ITMU_SUBSYSTEMCOLLECTOR_H__
#define __L1ITMU_SUBSYSTEMCOLLECTOR_H__
// 
// Class: L1ITMu::SubsystemCollector
//
// Info: This is the base class for a object that eats a specified subsystem
//       and turns those digis into L1ITMu::TriggerPrimitives
//
// Author: L. Gray (FNAL)
//
#include <vector>
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace L1ITMu {
  
  class SubsystemCollector {
  public:
    SubsystemCollector(const edm::ParameterSet&);
    virtual ~SubsystemCollector() {}

    virtual void extractPrimitives(const edm::Event&, const edm::EventSetup&, 
				   std::vector<TriggerPrimitive>&) const = 0;
  protected:
    edm::InputTag _src;
    edm::InputTag _src2;
  };
}

#endif
