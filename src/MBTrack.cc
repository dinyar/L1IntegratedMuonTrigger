#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"



///  construct starting froma DTTF candidate
L1ITMu::MBTrack::MBTrack(const L1MuDTTrackCand& dttrack) :
  L1MuRegionalCand(dttrack) {
  _wheel = dttrack.whNum();
  _sector = dttrack.scNum();
}

void L1ITMu::MBTrack::addStub(const MBLTCollection& stub)
{ 
  _associatedStubs[stub.detId()] = stub;
}

void L1ITMu::MBTrack::addStub(const std::pair<const DTChamberId, L1ITMu::MBLTCollection>& stub)
{ 
  _associatedStubs[stub.first] = stub.second;
}
