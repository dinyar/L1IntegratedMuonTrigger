#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"



///  construct starting from a DTTF candidate
L1ITMu::MBTrack::MBTrack(const L1MuDTTrackCand& dttrack) :
  L1MuRegionalCand(dttrack) {
  _wheel = dttrack.whNum();
  _sector = dttrack.scNum();
}

void L1ITMu::MBTrack::addStub(const MBLTContainerRefPair& stub)
{ 
  _associatedStubs.push_back( stub );
}

// void L1ITMu::MBTrack::addStub(const std::pair<const DTChamberId, L1ITMu::MBLTCollection>& stub)
// { 
//   _associatedStubs.push_back( stub.second );
// }

/// associate a GMTout 
void L1ITMu::MBTrack::associateGMTout(const L1MuGMTExtendedCand& gmtout)
{
  _gmtout.push_back(gmtout);
}
   
/// associate a GMTin 
void L1ITMu::MBTrack::associateGMTin(const L1MuRegionalCand& gmtin)
{
  _gmtin.push_back(gmtin);  
}
