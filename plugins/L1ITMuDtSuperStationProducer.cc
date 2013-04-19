
// framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/MakerMacros.h"

// user include files
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStation.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStationFwd.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

class L1ITMuDtSuperStationProducer : public edm::EDProducer {
public:
  L1ITMuDtSuperStationProducer(const edm::ParameterSet&);
  ~L1ITMuDtSuperStationProducer() {}

  void produce(edm::Event&, const edm::EventSetup&);  
private:
  edm::InputTag _trigPrimSrc;
  double _maxDeltaPhi;
};


L1ITMuDtSuperStationProducer::L1ITMuDtSuperStationProducer( const edm::ParameterSet& iConfig )
  : _trigPrimSrc(iConfig.getParameter<edm::InputTag>("TriggerPrimitiveSrc")),
    _maxDeltaPhi(iConfig.getParameter<double>("MaxDeltaPhi"))
{
  produces<L1ITMu::DtSuperStationMap>().setBranchAlias("DtSuperStationMap");
}

void 
L1ITMuDtSuperStationProducer::produce( edm::Event& iEvent, 
				       const edm::EventSetup& iSetup )
{

  std::auto_ptr<L1ITMu::DtSuperStationMap>
    out ( new L1ITMu::DtSuperStationMap );
  L1ITMu::DtSuperStationMap & tracksMap = *out;

  edm::Handle<L1ITMu::TriggerPrimitiveCollection> tps;
  iEvent.getByLabel(_trigPrimSrc, tps);

  auto tp = tps->cbegin();
  auto tpbeg = tps->cbegin();
  auto tpend = tps->cend();
  for( ; tp != tpend; ++tp ) {

    DTChamberId key;

    L1ITMu::TriggerPrimitive::subsystem_type type = tp->subsystem();
    switch ( type ) {

    case L1ITMu::TriggerPrimitive::kDT :
      key = tp->detId<DTChamberId>();
      break;

    case L1ITMu::TriggerPrimitive::kRPC : {
      if ( tp->detId<RPCDetId>().region() ) continue; // endcap
      int station = tp->detId<RPCDetId>().station();
      int sector  = tp->detId<RPCDetId>().sector();
      int wheel = tp->detId<RPCDetId>().ring();
      key = DTChamberId( wheel, station, sector );
      break;
    }

    default : continue;
    }

    if ( out->find( key ) == out->end() ) {
      tracksMap[key] = L1ITMu::DtSuperStation( key );
    }

    L1ITMu::TriggerPrimitiveRef tpref(tps, tp - tpbeg);
    tracksMap[key].addStub( tpref );
  }


  L1ITMu::DtSuperStationMap::iterator st = out->begin();
  L1ITMu::DtSuperStationMap::iterator stend = out->end();
  for ( ; st != stend; ++st ) st->second.associate( _maxDeltaPhi );

  iEvent.put(out);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuDtSuperStationProducer);
