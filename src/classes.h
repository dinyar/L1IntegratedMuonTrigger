#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/InternalTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/InternalTrackFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/CandidateTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/CandidateTrackFwd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/RegionalTracksFwd.h"

#include "DataFormats/RPCDigi/interface/RPCDigiL1Link.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackContainer.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStation.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStationFwd.h"

namespace {
  using namespace L1ITMu;
  struct dictionary {
    // trigger primitive information
    TriggerPrimitive::RPCData rpcD;
    TriggerPrimitive::CSCData cscD;
    TriggerPrimitive::DTData  dtD;
    
    edm::Wrapper<TriggerPrimitive::RPCData> wrpcD;
    edm::Wrapper<TriggerPrimitive::CSCData> wcscD;
    edm::Wrapper<TriggerPrimitive::DTData>  wdtD;

    TriggerPrimitive trigPrim;
    edm::Wrapper<TriggerPrimitive> wtrigPrim;
    TriggerPrimitiveCollection vtrigPrim;
    edm::Wrapper<TriggerPrimitiveCollection> wvtrigPrim;

    TriggerPrimitiveRef rtrigPrim;
    TriggerPrimitiveList rvtrigPrim;
    TriggerPrimitiveStationMap smtrigPrim;
    
    TriggerPrimitivePtr ptrigPrim;   

    // internal tracks
    InternalTrack tk;
    edm::Wrapper<InternalTrack> wtk;

    InternalTrackCollection ctk;
    edm::Wrapper<InternalTrackCollection> wctk;

    InternalTrackRef rtk;
    InternalTrackPtr ptk;

    // subdetector tracks
    DtSuperStation dss;
    edm::Wrapper<DtSuperStation> wdss;

    DtSuperStationPair pdss;
    edm::Wrapper<DtSuperStationPair> wpdss;

    DtSuperStationMap mdss;
    edm::Wrapper<DtSuperStationMap> wmdss;

    DtSuperStation::primitiveAssociation sspa;
    edm::Wrapper<DtSuperStation::primitiveAssociation> wsspa;

    std::vector<DtSuperStation::primitiveAssociation> vsspa;
    edm::Wrapper<std::vector<DtSuperStation::primitiveAssociation> > wvsspa;

    // Candidate tracks
    CandidateTrack cctk;
    edm::Wrapper<CandidateTrack> wcctk;

    CandidateTrackCollection ccctk;
    edm::Wrapper<CandidateTrackCollection> wccctk;

    CandidateTrackRef rcctk;
    CandidateTrackPtr pcctk;
    
    // regional muon system tracks removed from their containers
    // so we can save references to them
    RegionalCandBaseRef rcR2B;
    RegionalCandPtr     rcPtr;
    RegionalCandRef     rfRef;

    DTTrackCollection dtTrkColl;
    edm::Wrapper<DTTrackCollection> wdtTrkColl;
    DTTrackPtr dtTrkPtr;
    DTTrackRef dtTrackRef;
    
    CSCTrackCollection cscTrkColl;
    edm::Wrapper<CSCTrackCollection> wcscTrkColl;
    CSCTrackPtr cscTrkPtr;
    CSCTrackRef cscTrkRef;

    RPCL1LinkPtr prpcL1link;
    RPCL1LinkRef rrpcL1link;

    edm::reftobase::Holder<L1MuRegionalCand,RegionalCandRef>  r2rholder;
    edm::reftobase::Holder<L1MuRegionalCand,DTTrackRef>  r2dtholder;
    edm::reftobase::Holder<L1MuRegionalCand,CSCTrackRef>  r2cscholder;
    edm::reftobase::Holder<L1MuRegionalCand,InternalTrackRef>  r2itholder;
    edm::reftobase::RefHolder<RegionalCandRef>  r2rrefholder;
    edm::reftobase::RefHolder<DTTrackRef>  r2dtrefholder;
    edm::reftobase::RefHolder<CSCTrackRef>  r2cscrefholder;
    edm::reftobase::RefHolder<InternalTrackRef>  r2itrefholder;
  };
}
