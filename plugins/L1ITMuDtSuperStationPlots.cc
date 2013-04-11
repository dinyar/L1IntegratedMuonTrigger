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

    timingConf[i] = _fs->make<TH1F>( "timing_"+stname, "dt vs rpc timing in "+stname, 4, -2, 2 );
    timingConfIn[i] = _fs->make<TH1F>( "timingIn_"+stname, "dt vs rpcIn timing in "+stname, 3, -1, 2 );
    timingConfOut[i] = _fs->make<TH1F>( "timingOut_"+stname, "dt vs rpcOut timing in "+stname, 3, -1, 2 );


    dtDist[i] = _fs->make<TH2F>( "dtDist_"+stname, "dt primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcInDist[i] = _fs->make<TH2F>( "rpcInDist_"+stname, "rpc inner primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    rpcOutDist[i] = _fs->make<TH2F>( "rpcOutDist_"+stname, "rpc outer primitives in "+stname, 5, -2, 3, 12, 1, 13 );

  }

  deltaPhi = _fs->make<TH1F>( "deltaPhi", "deltaPhi", 40, 0, 0.4 );
  deltaEta = _fs->make<TH1F>( "deltaEta", "deltaEta", 40, 0, 0.4 );
  deltaR = _fs->make<TH1F>( "deltaR", "deltaR", 100, 0, 1 );
  deltaPhiR = _fs->make<TH2F>( "deltaPhiR", "deltaR vs deltaPhi", 40, 0, 0.4, 100, 0, 0.5 );

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

      if ( rpcIn || rpcOut ) {

	TriggerPrimitiveList::const_iterator dt = superSt.getDtSegments().begin();
	TriggerPrimitiveList::const_iterator dtend = superSt.getDtSegments().end();

	TriggerPrimitiveList::const_iterator rpcInEnd = superSt.getRpcInner().end();
	TriggerPrimitiveList::const_iterator rpcOutEnd = superSt.getRpcOuter().end();
	TriggerPrimitiveList::const_iterator rpcIn = rpcInEnd;
	TriggerPrimitiveList::const_iterator rpcOut = rpcOutEnd;

	for ( ; dt != dtend; ++dt ) {

	  double eta = (*dt)->getCMSGlobalEta();
	  double phi = (*dt)->getCMSGlobalPhi();

	  double minRpcPhi = 99;
	  TriggerPrimitiveList::const_iterator rpcInIt = superSt.getRpcInner().begin();
	  for ( ; rpcInIt != rpcInEnd; ++rpcInIt ) {

	    double rpcPhi = (*rpcInIt)->getCMSGlobalPhi();
	    if ( fabs( phi - rpcPhi ) < minRpcPhi ) {
	      rpcIn = rpcInIt;
	      minRpcPhi = rpcPhi;
	    }

	  }

	  minRpcPhi = 99;
	  TriggerPrimitiveList::const_iterator rpcOutIt = superSt.getRpcOuter().begin();
	  for ( ; rpcOutIt != rpcOutEnd; ++rpcOutIt ) {

	    double rpcPhi = (*rpcOutIt)->getCMSGlobalPhi();
	    if ( fabs( phi - rpcPhi ) < minRpcPhi ) {
	      rpcOut = rpcOutIt;
	      minRpcPhi = rpcPhi;
	    }
	  }

	  if ( rpcIn != rpcInEnd ) {
	    double rpcPhi = (*rpcIn)->getCMSGlobalPhi();
	    double rpcEta = (*rpcIn)->getCMSGlobalEta();
	    double dR = reco::deltaR( eta, phi, rpcEta, rpcPhi);
	    deltaEta->Fill( fabs( eta - rpcEta ) );
	    deltaPhi->Fill( fabs( phi - rpcPhi ) );
	    deltaR->Fill( dR );
	    deltaPhiR->Fill( fabs( phi - rpcPhi ), dR );
	  }

	  if ( rpcOut != rpcOutEnd ) {
	    double rpcPhi = (*rpcOut)->getCMSGlobalPhi();
	    double rpcEta = (*rpcOut)->getCMSGlobalEta();
	    double dR = reco::deltaR( eta, phi, rpcEta, rpcPhi);
	    deltaEta->Fill( fabs( eta - rpcEta ) );
	    deltaPhi->Fill( fabs( phi - rpcPhi ) );
	    deltaR->Fill( dR );
	    deltaPhiR->Fill( fabs( phi - rpcPhi ), dR );
	  }


	  if ( rpcIn != rpcInEnd && rpcOut != rpcOutEnd ) {
	    confirmed[index]->Fill( 1 );

	    int dtbx = (*dt)->getBX();
	    int rpcbx = (*rpcIn)->getBX();
	    int rpcbxOut = (*rpcOut)->getBX();
	    int bx = 0;
	    if ( rpcbxOut != rpcbx ) bx = -2;
	    else if ( dtbx == rpcbx ) bx = 0;
	    else if ( dtbx > rpcbx ) bx = 1;
	    else if ( dtbx < rpcbx ) bx = -1;
	    timingConf[index]->Fill( bx );

	  } else if ( rpcIn != rpcInEnd ) {
	    confirmed[index]->Fill( 2 );

	    int dtbx = (*dt)->getBX();
	    int rpcbx = (*rpcIn)->getBX();
	    int bx = 0;
	    if ( dtbx == rpcbx ) bx = 0;
	    else if ( dtbx > rpcbx ) bx = 1;
	    else if ( dtbx < rpcbx ) bx = -1;
	    timingConfIn[index]->Fill( bx );

	  } else if ( rpcOut != rpcOutEnd ) {
	    confirmed[index]->Fill( 3 );

	    int dtbx = (*dt)->getBX();
	    int rpcbx = (*rpcOut)->getBX();
	    int bx = 0;
	    if ( dtbx == rpcbx ) bx = 0;
	    else if ( dtbx > rpcbx ) bx = 1;
	    else if ( dtbx < rpcbx ) bx = -1;
	    timingConfOut[index]->Fill( bx );

	  }  else confirmed[index]->Fill( 0 );

	}
      }
      else confirmed[index]->Fill( 0 );

    }


  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuDtSuperStationPlots);
