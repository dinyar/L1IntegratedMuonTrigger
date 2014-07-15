#include "L1Trigger/L1IntegratedMuonTrigger/interface/SubsystemCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace L1ITMu;

SubsystemCollector::SubsystemCollector(const edm::ParameterSet& ps):
  _src(ps.getParameter<edm::InputTag>("src")) ,
  _src2(ps.getParameter<edm::InputTag>("src2")) {
}
