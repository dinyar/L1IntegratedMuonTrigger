// 
// Class: L1ITMuDtSuperStationPlots
//
// Info: Processes a track into histograms of delta-phis and such
//
// Author: L. Gray (FNAL)
//

#include <memory>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStation.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/DtSuperStationFwd.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

using namespace L1ITMu;

typedef edm::ParameterSet PSet;

class L1ITMuDtSuperStationPlots : public edm::EDAnalyzer {  
  
public:
  L1ITMuDtSuperStationPlots( const PSet& );
  ~L1ITMuDtSuperStationPlots();

  void analyze( const edm::Event&, const edm::EventSetup& );  
private:
  edm::InputTag _dtSuperStInput;
  edm::Service<TFileService> _fs;
  TH1F * confirmed[4];
  TH1F * timing[4];

  TH2F * dtDist[4];
  TH2F * rpcInDist[4];
  TH2F * rpcOutDist[4];

};


L1ITMuDtSuperStationPlots::L1ITMuDtSuperStationPlots(const PSet& p)
  : _dtSuperStInput( p.getParameter<edm::InputTag>("L1ITMuDtSuperStation") )
{

  for ( int i = 0; i < 4; ++i ) {
    TString stname = Form( "st%d", i+1 );
    confirmed[i] = _fs->make<TH1F>( "confirmed_"+stname, "confirmed in "+stname, 4, 0, 4 );

    confirmed[i]->GetXaxis()->SetBinLabel( 1, "Unconfirmed");
    confirmed[i]->GetXaxis()->SetBinLabel( 2, "In+Out");
    confirmed[i]->GetXaxis()->SetBinLabel( 3, "In");
    confirmed[i]->GetXaxis()->SetBinLabel( 4, "Out");

    timing[i] = _fs->make<TH1F>( "timing_"+stname, "out of time dt segments in "+stname, 4, 0, 4 );


    dtDist[i] = _fs->make<TH2F>( "dtDist_"+stname, "dt primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcInDist[i] = _fs->make<TH2F>( "rpcInDist_"+stname, "rpc inner primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcOutDist[i] = _fs->make<TH2F>( "rpcOutDist_"+stname, "rpc outer primitives in "+stname, 5, -2, 3, 12, 1, 13 );


  }


}
  
L1ITMuDtSuperStationPlots::~L1ITMuDtSuperStationPlots()
{
}


void L1ITMuDtSuperStationPlots::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::DtSuperStationMap> dtSuperStMap;
  iEvent.getByLabel( _dtSuperStInput, dtSuperStMap );

  L1ITMu::DtSuperStationMap::const_iterator st = dtSuperStMap->begin();
  L1ITMu::DtSuperStationMap::const_iterator stend = dtSuperStMap->end();

  for ( ; st != stend; ++st ) {

    const L1ITMu::DtSuperStation & superSt = st->second;

    size_t dtseg = superSt.getDtSegments().size();
    size_t rpcIn = superSt.getRpcInner().size();
    size_t rpcOut = superSt.getRpcOuter().size();

    int station = superSt.station();
    int index = station - 1;
    int wheel = superSt.wheel();
    int sector = superSt.sector();

    if ( index < 0 || index > 3 )
      throw cms::Exception("Invalid Station")
	<< "wrong station number " << station << std::endl;

    if ( dtseg ) dtDist[index]->Fill( wheel, sector );
    if ( rpcIn ) rpcInDist[index]->Fill( wheel, sector );
    if ( rpcOut ) rpcOutDist[index]->Fill( wheel, sector );

    if ( dtseg ) {

      if ( rpcIn && rpcOut ) confirmed[index]->Fill( 1 );
      else if ( rpcIn )      confirmed[index]->Fill( 2 );
      else if ( rpcOut )     confirmed[index]->Fill( 3 );
      else confirmed[index]->Fill( 0 );


      if ( rpcIn || rpcOut ) {

	TriggerPrimitiveList::const_iterator it = superSt.getDtSegments().begin();
	TriggerPrimitiveList::const_iterator itend = superSt.getDtSegments().end();
	int bx = 0;

	for ( ; it != itend; ++it ) {
	  bx += ( (*it)->getBX() ? 1 : 0);
	}
	timing[index]->Fill( bx );
      }

    }


  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuDtSuperStationPlots);
