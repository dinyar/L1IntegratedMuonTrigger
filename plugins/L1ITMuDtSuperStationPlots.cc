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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
  TH1F * timingConf[4];
  TH1F * timingConfIn[4];
  TH1F * timingConfOut[4];

  TH2F * dtDist[4];
  TH2F * rpcInDist[4];
  TH2F * rpcOutDist[4];

  TH1F * deltaPhi;
  TH1F * deltaEta;
  TH1F * deltaR;
  TH2F * deltaPhiR;

  TH2F * dtQuality;
  TH1F * rpcInHitsPerDtseg;
  TH1F * rpcOutHitsPerDtseg;

};


L1ITMuDtSuperStationPlots::L1ITMuDtSuperStationPlots(const PSet& p)
  : _dtSuperStInput( p.getParameter<edm::InputTag>("L1ITMuDtSuperStation") )
{

  for ( int i = 0; i < 4; ++i ) {
    TString stname = Form( "st%d", i+1 );
    confirmed[i] = _fs->make<TH1F>( "confirmed_"+stname, "confirmed in "+stname, 5, -1, 4 );

    confirmed[i]->GetXaxis()->SetBinLabel( 1, "RPC Only");
    confirmed[i]->GetXaxis()->SetBinLabel( 2, "Unconfirmed");
    confirmed[i]->GetXaxis()->SetBinLabel( 3, "In+Out");
    confirmed[i]->GetXaxis()->SetBinLabel( 4, "In");
    confirmed[i]->GetXaxis()->SetBinLabel( 5, "Out");

    timingConf[i] = _fs->make<TH1F>( "timing_"+stname, "dt vs rpc timing in "+stname, 4, -2, 2 );
    timingConfIn[i] = _fs->make<TH1F>( "timingIn_"+stname, "dt vs rpcIn timing in "+stname, 3, -1, 2 );
    timingConfOut[i] = _fs->make<TH1F>( "timingOut_"+stname, "dt vs rpcOut timing in "+stname, 3, -1, 2 );


    dtDist[i] = _fs->make<TH2F>( "dtDist_"+stname, "dt primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcInDist[i] = _fs->make<TH2F>( "rpcInDist_"+stname, "rpc inner primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcOutDist[i] = _fs->make<TH2F>( "rpcOutDist_"+stname, "rpc outer primitives in "+stname, 5, -2, 3, 12, 1, 13 );

  }

  deltaPhi = _fs->make<TH1F>( "deltaPhi", "deltaPhi", 400, -0.2, 0.2 );
  deltaEta = _fs->make<TH1F>( "deltaEta", "deltaEta", 40, 0, 0.4 );
  deltaR = _fs->make<TH1F>( "deltaR", "deltaR", 100, 0, 1 );
  rpcInHitsPerDtseg = _fs->make<TH1F>( "rpcInHitsPerDtseg", "number of inner rpc hits matching a dt sector", 20, 0, 20 );
  rpcOutHitsPerDtseg = _fs->make<TH1F>( "rpcOutHitsPerDtseg", "number of outer rpc hits matching a dt sector", 20, 0, 20 );

  deltaPhiR = _fs->make<TH2F>( "deltaPhiR", "deltaR vs deltaPhi", 400, 0, 0.4, 100, 0, 0.5 );
  dtQuality = _fs->make<TH2F>( "dtQuality", "dt quality vs rpc match", 4, 0, 4, 8, 1, 7);
  dtQuality->GetXaxis()->SetBinLabel( 1, "Unconfirmed");
  dtQuality->GetXaxis()->SetBinLabel( 2, "In+Out");
  dtQuality->GetXaxis()->SetBinLabel( 3, "In");
  dtQuality->GetXaxis()->SetBinLabel( 4, "Out");

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

    /// useful index
    int station = superSt.station();
    int index = station - 1;
    int wheel = superSt.wheel();
    int sector = superSt.sector();

    if ( index < 0 || index > 3 )
      throw cms::Exception("Invalid Station")
	<< "wrong station number " << station << std::endl;

    /// size for primitives vectors
    size_t dtListSize = superSt.getDtSegments().size();
    size_t rpcInListSize = superSt.getRpcInner().size();
    size_t rpcOutListSize = superSt.getRpcOuter().size();

    /// fill general distribution plots
    if ( dtListSize ) dtDist[index]->Fill( wheel, sector );
    if ( rpcInListSize ) rpcInDist[index]->Fill( wheel, sector );
    if ( rpcOutListSize ) rpcOutDist[index]->Fill( wheel, sector );

    if ( !dtListSize ) {
      confirmed[index]->Fill( -1 );
      continue;
    }

    /// get dt to rpc associations
    for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {

      const L1ITMu::TriggerPrimitiveRef & dt = superSt.getDtSegments().at(iDt);
      // const L1ITMu::TriggerPrimitive & dt = *(superSt.getDtSegments().at( iDt ) );
      double eta = dt->getCMSGlobalEta();
      double phi = dt->getCMSGlobalPhi();
      int dtbx = dt->getBX();

      L1ITMu::TriggerPrimitiveList rpcInMatch = superSt.getRpcInAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList rpcOutMatch = superSt.getRpcOutAssociatedStubs( iDt );
      size_t rpcInMatchSize = rpcInMatch.size();
      size_t rpcOutMatchSize = rpcOutMatch.size();


      /// count matching
      rpcInHitsPerDtseg->Fill( rpcInMatchSize );
      rpcOutHitsPerDtseg->Fill( rpcOutMatchSize );

      if ( rpcInMatchSize && rpcOutMatchSize ) {
	confirmed[index]->Fill( 1 );
	dtQuality->Fill( 1, dt->getDTData().qualityCode );
      } else if ( rpcInMatchSize ) {
	confirmed[index]->Fill( 2 );
	dtQuality->Fill( 2, dt->getDTData().qualityCode );
      } else if ( rpcOutMatchSize ) {
	confirmed[index]->Fill( 3 );
	dtQuality->Fill( 3, dt->getDTData().qualityCode );
      } else {
	confirmed[index]->Fill( 0 );
	dtQuality->Fill( 0., dt->getDTData().qualityCode );
      }


      /// angular differences
      TriggerPrimitiveList::const_iterator rpcInEnd = rpcInMatch.end();
      TriggerPrimitiveList::const_iterator rpcOutEnd = rpcOutMatch.end();

      TriggerPrimitiveList::const_iterator rpcInIt = rpcInMatch.begin();
      for ( ; rpcInIt != rpcInEnd; ++rpcInIt ) {
	double rpcPhi = (*rpcInIt)->getCMSGlobalPhi();
	double rpcEta = (*rpcInIt)->getCMSGlobalEta();
	double dR = reco::deltaR( eta, phi, rpcEta, rpcPhi);
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcPhi );
	deltaEta->Fill( fabs( eta - rpcEta ) );
	deltaPhi->Fill( rpcDeltaPhi );
	deltaR->Fill( dR );
	deltaPhiR->Fill( fabs( rpcDeltaPhi ), dR );

      }

      TriggerPrimitiveList::const_iterator rpcOutIt = rpcOutMatch.begin();
      for ( ; rpcOutIt != rpcOutEnd; ++rpcOutIt ) {
	double rpcPhi = (*rpcOutIt)->getCMSGlobalPhi();
	double rpcEta = (*rpcOutIt)->getCMSGlobalEta();
	double dR = reco::deltaR( eta, phi, rpcEta, rpcPhi);
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcPhi );
	deltaEta->Fill( fabs( eta - rpcEta ) );
	deltaPhi->Fill( rpcDeltaPhi );
	deltaR->Fill( dR );
	deltaPhiR->Fill( fabs( rpcDeltaPhi ), dR );

      }


      /// timing differences
      if ( rpcInMatchSize && rpcOutMatchSize ) {
	int rpcbxIn = rpcInMatch.front()->getBX();
	int rpcbxOut = rpcOutMatch.front()->getBX();
	int bx = 0;
	if ( rpcbxOut != rpcbxIn ) bx = -2;
	else if ( dtbx == rpcbxIn ) bx = 0;
	else if ( dtbx > rpcbxIn ) bx = 1;
	else if ( dtbx < rpcbxIn ) bx = -1;
	timingConf[index]->Fill( bx );
      } else if ( rpcInMatchSize ) {
	int rpcbx = rpcInMatch.front()->getBX();
	int bx = 0;
	if ( dtbx == rpcbx ) bx = 0;
	else if ( dtbx > rpcbx ) bx = 1;
	else if ( dtbx < rpcbx ) bx = -1;
	timingConfIn[index]->Fill( bx );
      } else if ( rpcOutMatchSize ) {
	int rpcbx = rpcOutMatch.front()->getBX();
	int bx = 0;
	if ( dtbx == rpcbx ) bx = 0;
	else if ( dtbx > rpcbx ) bx = 1;
	else if ( dtbx < rpcbx ) bx = -1;
	timingConfOut[index]->Fill( bx );
      }

    }

  }


}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuDtSuperStationPlots);
