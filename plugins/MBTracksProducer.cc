// 
// Class: L1DTTFTrackConverter
//
// Info: This producer eats DTTF tracks (pre GMT) and matches them to 
//       L1ITMu::TriggerPrimitives. In the process of doing so it
//       converts the DTTF tracks into a collection L1ITMu::InternalTrack
//       each of which contains the track stubs it was matched to.
//
// Author: 
//

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrackFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/RegionalTracksFwd.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackContainer.h"

#include "DataFormats/Common/interface/RefProd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBhelpers.h"

// this magic file contains a DT TrackClass -> mode LUT
#include "L1Trigger/DTTrackFinder/src/L1MuDTTrackAssParam.h"

// GMT
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"


// using namespace L1ITMu;

typedef edm::ParameterSet PSet;

class MBTracksProducer : public edm::EDProducer {
public:
  MBTracksProducer(const PSet&);
  ~MBTracksProducer() {}

  void produce(edm::Event&, const edm::EventSetup&);  
private:
  edm::InputTag _dtTrackSrc, _gmtSrc, _trigPrimSrc;
  const int min_bx, max_bx;
};

MBTracksProducer::MBTracksProducer(const PSet& ps):
  _dtTrackSrc(ps.getParameter<edm::InputTag>("DTTrackSrc")),
  _gmtSrc(ps.getParameter<edm::InputTag>("gtDigis")),
  _trigPrimSrc(ps.getParameter<edm::InputTag>("TriggerPrimitiveSrc")),
  min_bx(ps.getParameter<int>("BX_min")),
  max_bx(ps.getParameter<int>("BX_max")) {
  produces<L1ITMu::DTTrackCollection>("input");
  produces<L1ITMu::MBTrackCollection>();
}

void MBTracksProducer::produce(edm::Event& ev, 
				    const edm::EventSetup& es) {
  std::auto_ptr<L1ITMu::MBTrackCollection> 
    convertedTracks (new L1ITMu::MBTrackCollection());
  std::auto_ptr<L1ITMu::DTTrackCollection> inputTracks(new L1ITMu::DTTrackCollection);

  // get the RefProd so we can make persistable references to
  // the track the converted track was made from
  edm::RefProd<L1ITMu::DTTrackCollection> dttfProd = 
    ev.getRefBeforePut<L1ITMu::DTTrackCollection>("input");

  edm::Handle<L1MuDTTrackContainer> dtTracks;
  ev.getByLabel(_dtTrackSrc,dtTracks);
  
  edm::Handle<L1ITMu::MBLTContainer> trigPrims;
  ev.getByLabel(_trigPrimSrc,trigPrims);

  // get GMT readout collection
  edm::Handle<L1MuGMTReadoutCollection> pCollection;
  ev.getByLabel(_gmtSrc, pCollection);            
  L1MuGMTReadoutCollection const* gmtrc = pCollection.product();
  
  // get GMT readout record
  std::vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
  std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr = gmt_records.begin();
  std::vector<L1MuGMTReadoutRecord>::const_iterator RREnd = gmt_records.end();
  
  int wheel;
  // DT sector processors have wheels [-3,-2,-1,1,2,3] since 
  // wheel zero needs two SPs
  for( int sp_wheel = -3 ; sp_wheel <= 3; ++sp_wheel ) {
    if( sp_wheel == 0 ) continue;
    wheel = std::abs(sp_wheel)-1;    
    wheel = sp_wheel < 0 ? -wheel : wheel;
    for( int sector = 0; sector <= 11; ++sector ) {
      for( int bx = min_bx; bx <= max_bx; ++bx ) {
	for( int itrk = 1; itrk <=2; ++itrk ) {
	  std::unique_ptr<L1MuDTTrackCand> dttrk;
	  if( itrk == 1 ) 
	    dttrk.reset(dtTracks->dtTrackCand1(sp_wheel,sector,bx));
	  else            
	    dttrk.reset(dtTracks->dtTrackCand2(sp_wheel,sector,bx));
	  
	  if( dttrk ) {
	    inputTracks->push_back(*dttrk);


	    L1ITMu::MBTrack trk(*dttrk);
	    L1ITMu::DTTrackRef parentRef(dttfProd,inputTracks->size() - 1);
	    L1ITMu::RegionalCandBaseRef parentBaseRef(parentRef);
	    trk.setParent(parentBaseRef);
            
            
            /// JP: GMT-DTTF matching            
            for ( ; RRItr != RREnd; ++RRItr ) {
              
              // loop over GMT input collection
              std::vector<L1MuRegionalCand> dttfCands = RRItr->getDTBXCands();
              std::vector<L1MuRegionalCand>::iterator dttfCand = dttfCands.begin();
              std::vector<L1MuRegionalCand>::iterator dttfCandEnd = dttfCands.end();
              
              for( ; dttfCand != dttfCandEnd; ++dttfCand ) {
                
                if(dttfCand->empty()) continue;                 
                if( (dttfCand->quality() == dttrk->quality()) &&
                    (fabs(dttfCand->phiValue() - dttrk->phiValue()) < 1e-5) ) {
                
                  trk.associateGMTin(*dttfCand);
                    
                }                                
              }

              // loop over GMT output collection              
              std::vector<L1MuGMTExtendedCand> gmtCands = RRItr->getGMTBrlCands();
              std::vector<L1MuGMTExtendedCand>::iterator gmtCand = gmtCands.begin();
              std::vector<L1MuGMTExtendedCand>::iterator gmtCandEnd = gmtCands.end();
              
              for( ; gmtCand != gmtCandEnd; ++gmtCand ) {
                
                if(gmtCand->empty()) continue;        
                
                float dphi = gmtCand->phiValue() - dttrk->phiValue();
                float deta = gmtCand->etaValue() - dttrk->etaValue();
                float dr   = sqrt(dphi*dphi + deta*deta);
                if( dr < 0.3 ) {

                  trk.associateGMTout(*gmtCand);               
                    
                }
              } 
              
            }

            /// JP: MB-DTTF matching                                  
            std::vector<unsigned> addrs;
            addrs.reserve(4);	     	   
	    for( int station = 1; station <= 4; ++ station ) {
	      addrs.push_back(dttrk->stNum(station));
	    }	

	    // in DTs the mode is encoded by the track class
	    // mode is a 4 bit word , the bit position indicates the station
	    // if the bit is 1 then the station was used in track building
	    const unsigned mode = tc2bitmap((TrackClass)dttrk->TCNum());

            /// JP : implementation through the MBhelpers class
            ///      the same code could be placed somewhere else
            ///      to be more consistent
	    L1ITMu::MBLTContainer tplist =  L1ITMu::MBhelpers::getPrimitivesByMBTriggerInfo(wheel,
						    sp_wheel,sector+1,
						    trigPrims,mode,
						    addrs);
	                            
	    auto stub = tplist.cbegin();
	    auto stend = tplist.cend();
	    for( ; stub != stend; ++stub ) {
	      trk.addStub(*stub);      
	    }

	    convertedTracks->push_back(trk);
	  }
	  dttrk.release();
	}
      }
    }
  }
  ev.put(inputTracks,"input");
  ev.put(convertedTracks);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MBTracksProducer);
