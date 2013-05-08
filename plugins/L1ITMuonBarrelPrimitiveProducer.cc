// framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// L1IT include files
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/PrimitiveCombiner.h"

// user include files
#include "DataFormats/Math/interface/deltaPhi.h"

class L1ITMuonBarrelPrimitiveProducer : public edm::EDProducer {

public:
  L1ITMuonBarrelPrimitiveProducer(const edm::ParameterSet&);
  ~L1ITMuonBarrelPrimitiveProducer();
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  edm::InputTag _mbltCollectionInput;
   const L1ITMu::PrimitiveCombiner::resolutions _resol;

};


L1ITMuonBarrelPrimitiveProducer::~L1ITMuonBarrelPrimitiveProducer()
{
}

L1ITMuonBarrelPrimitiveProducer::L1ITMuonBarrelPrimitiveProducer( const edm::ParameterSet& iConfig )
  : _mbltCollectionInput( iConfig.getParameter<edm::InputTag>("MBLTCollection") ),
    _resol( iConfig.getParameter<double>("xDtResol"), 
	    iConfig.getParameter<double>("xRpcResol"),
	    iConfig.getParameter<double>("phibDtCorrResol"),
	    iConfig.getParameter<double>("phibDtUnCorrResol") )
{
  produces<L1MuDTChambPhContainer>();
  // produces<std::vector<L1MuDTChambPhDigi> >();
}



void L1ITMuonBarrelPrimitiveProducer::produce( edm::Event& iEvent, 
					      const edm::EventSetup& iSetup )
{

//   std::auto_ptr<std::vector<L1MuDTChambPhDigi> >
//     out ( new std::vector<L1MuDTChambPhDigi> );
  std::auto_ptr<L1MuDTChambPhContainer> out(new L1MuDTChambPhContainer);
  std::vector<L1MuDTChambPhDigi> phiChambVector;
  
  edm::Handle<L1ITMu::MBLTContainer> mbltContainer;
  iEvent.getByLabel( _mbltCollectionInput, mbltContainer );

  L1ITMu::MBLTContainer::const_iterator st = mbltContainer->begin();
  L1ITMu::MBLTContainer::const_iterator stend = mbltContainer->end();

  L1MuDTChambPhContainer phiContainer;
  std::vector<L1MuDTChambPhDigi> phiVector;

  for ( ; st != stend; ++st ) {

    const L1ITMu::MBLTCollection & mbltStation = st->second;

    /// useful index
    int station = mbltStation.station();
    int wheel = mbltStation.wheel();
    int sector = mbltStation.sector();
    ///

    /// get dt to rpc associations
    size_t dtListSize = mbltStation.getDtSegments().size();
    std::vector<size_t> uncorrelated;
    for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {

      const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);
      int dtquality = dt->getDTData().qualityCode;
      
      /// define new set of qualities
      /// skip for the moment uncorrelated
      int qualityCode = -2;
      switch ( dtquality ) {
      case -1 : continue;/// -1 are theta
      case 0 : qualityCode = -2; break;
      case 1 : qualityCode = -2; break;
      case 2 : uncorrelated.push_back( iDt ); continue;
      case 3 : uncorrelated.push_back( iDt ); continue;
      case 4 : qualityCode = 5; break;
      case 5 : qualityCode = 5; break;
      case 6 : qualityCode = 6; break;
      default : qualityCode = dtquality; break;
      }

      L1MuDTChambPhDigi chamb( dt->getBX(), wheel, sector, station, dt->getDTData().radialAngle,
			       dt->getDTData().bendingAngle, qualityCode,
			       dt->getDTData().Ts2TagCode, dt->getDTData().BxCntCode );
      phiChambVector.push_back( chamb );
    }


    /// now deal with uncorrelated
    /// basic idea:
    /// - use reduced list of array index for cpu saving
    /// - try to match only with closest DT segment (deltaPhi)
    /// - try to match only HI with HO if at different bx
    /// - match valid if they share the closest rpc hit (in/out/in+out)
    /// - remove used dt primitives ( set -1 the index at the corresponding position)

    size_t uncSize = uncorrelated.size(); 
    for ( size_t idxDt = 0; idxDt < uncSize; ++idxDt ) {

      int iDt = uncorrelated.at(idxDt);
      if ( iDt < 0 ) continue;
      const L1ITMu::TriggerPrimitive & dt = *mbltStation.getDtSegments().at(iDt);

      /// check if there is a pair of HI+HO at different bx
      int closest = -1;
      int closestIdx = -1;
      double minDeltaPhiDt = 9999999999;
      for ( size_t jdxDt = idxDt+1; jdxDt < uncSize; ++jdxDt ) {

	int jDt = uncorrelated.at(jdxDt);
	if ( jDt < 0 ) continue;

	const L1ITMu::TriggerPrimitiveRef & dtM = mbltStation.getDtSegments().at(jDt);
	if ( dt.getBX() == dtM->getBX() || dt.getDTData().qualityCode == dtM->getDTData().qualityCode )
	  continue;

	double deltaPhiDt = fabs( reco::deltaPhi( dt.getCMSGlobalPhi(), dtM->getCMSGlobalPhi() ) );
	if ( deltaPhiDt < minDeltaPhiDt ) {
	  closest = jDt;
	  closestIdx = jdxDt;
	}
      }

      /// check if the pair shares the closest rpc hit
      L1ITMu::MBLTCollection::bxMatch match = L1ITMu::MBLTCollection::NOMATCH;
      // if ( closest > 0 && minDeltaPhiDt < 0.05 ) {
      if ( closest > 0 ) {
	match = mbltStation.haveCommonRpc( iDt, closest );
      }

      /// this is just a seto of output variables for building L1ITMuDTChambPhDigi
      int qualityCode = dt.getDTData().qualityCode;
      int bx = -2;
      int radialAngle = 0;
      int bendingAngle = 0;
      L1ITMu::PrimitiveCombiner combiner( st->first, _resol );
      /// association HI/HO provided by the tool
      combiner.addDt( dt );

      /// there is a pair HI+HO with a shared inner RPC hit
      if ( match != L1ITMu::MBLTCollection::NOMATCH ) {
	uncorrelated[closestIdx] = -1;

	/// association HI/HO provided by the tool
	combiner.addDt( *mbltStation.getDtSegments().at(closest) );

	/// redefine quality
	qualityCode = 4;
	L1ITMu::TriggerPrimitiveList rpcInMatch = mbltStation.getRpcInAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );

	/// there is a pair HI+HO with a shared inner RPC hit
	if ( match == L1ITMu::MBLTCollection::INMATCH ) {

	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  combiner.addRpcIn( rpcIn );
	  bx = rpcIn.getBX();

	  /// there is a pair HI+HO with a shared outer RPC hit
	} else if ( match == L1ITMu::MBLTCollection::OUTMATCH ) {

	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  combiner.addRpcOut( rpcOut );
	  bx = rpcOut.getBX();

	  /// there is a pair HI+HO with both shared inner and outer RPC hit
	} else if ( match == L1ITMu::MBLTCollection::FULLMATCH ) {

	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  combiner.addRpcIn( rpcIn );
	  combiner.addRpcOut( rpcOut );
	  bx = rpcIn.getBX();
	}


      } else { /// there is no match

	bool hasRpcMatch = false;

	L1ITMu::TriggerPrimitiveList rpcInMatch = mbltStation.getRpcInAssociatedStubs( iDt );
	L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );
	size_t rpcInMatchSize = rpcInMatch.size();
	size_t rpcOutMatchSize = rpcOutMatch.size();

	/// the uncorrelated has possibly inner and outer confirmation
	if ( rpcInMatchSize && rpcOutMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  /// only the first is real...
	  if ( dt.getBX() == rpcIn.getBX() && dt.getBX() == rpcOut.getBX() ) {
	    bx = rpcIn.getBX();
	    hasRpcMatch = true;
	    combiner.addRpcIn( rpcIn );
	    combiner.addRpcOut( rpcOut );
	  } else if ( dt.getBX() == rpcIn.getBX() ) {
	    bx = rpcIn.getBX();
	    hasRpcMatch = true;
	    combiner.addRpcIn( rpcIn );
	  } else if ( dt.getBX() == rpcOut.getBX() ) {
	    bx = rpcOut.getBX();
	    hasRpcMatch = true;
	    combiner.addRpcOut( rpcOut );
	  }

	/// the uncorrelated has a possible inner confirmation
	} else if ( rpcInMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcIn = *rpcInMatch.front();
	  if ( dt.getBX() == rpcIn.getBX() ) {
	    bx = rpcIn.getBX();
	    hasRpcMatch = true;
	    combiner.addRpcIn( rpcIn );
	  }

	/// the uncorrelated has a possible outer confirmation
	} else if ( rpcOutMatchSize ) {
	  const L1ITMu::TriggerPrimitive & rpcOut = *rpcOutMatch.front();
	  if ( dt.getBX() == rpcOut.getBX() ) {
	    bx = rpcOut.getBX();
	    hasRpcMatch = true;
	    combiner.addRpcOut( rpcOut );
	  }

	}

	// match found, PrimitiveCombiner has the needed variables already calculated
	if ( hasRpcMatch ) {
	  combiner.combine();
	  radialAngle = combiner.radialAngle();
	  bendingAngle = combiner.bendingAngle();
	} else {
	  // no match found, keep the primitive as it is
	  bx = dt.getBX();
	  radialAngle = dt.getDTData().radialAngle;
	  bendingAngle = dt.getDTData().bendingAngle;
	  qualityCode = ( qualityCode == 2 ) ? 0 : 1;
	}

      }
      L1MuDTChambPhDigi chamb( bx, wheel, sector, station, radialAngle,
			       bendingAngle, qualityCode,
			       dt.getDTData().Ts2TagCode, dt.getDTData().BxCntCode );
      phiChambVector.push_back( chamb );

    } /// end of the Uncorrelated loop



  }

  out->setContainer( phiChambVector );
  /// fill event
  iEvent.put(out);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelPrimitiveProducer);

