// 
// Class: L1ITMuonBarrelPlots
//
// Info: Processes a track into histograms of delta-phis and such
//
// Author: 
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

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrackFwd.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

using namespace L1ITMu;

typedef edm::ParameterSet PSet;

class L1ITMuonBarrelPlots : public edm::EDAnalyzer {  
  
public:
  L1ITMuonBarrelPlots( const PSet& );
  ~L1ITMuonBarrelPlots();

  void analyze( const edm::Event&, const edm::EventSetup& );  
private:
  void endJob();


  edm::InputTag _mbltCollectionInput;
  edm::InputTag _mbTracksCollectionInput;
  edm::InputTag _l1itDtPhiChInput;
  edm::Service<TFileService> _fs;
  TH1F * confirmed[4];
  TH1F * timingConf[4];
  TH1F * timingConfIn[4];
  TH1F * timingConfOut[4];

  TH2F * dtDist[4];
  TH2F * rpcInDist[4];
  TH2F * rpcOutDist[4];
  TH2F * dtQuality[4];

  TH1F * deltaPhi;
  TH1F * deltaEta;
  TH1F * deltaR;
  TH2F * deltaPhiR;
  TH1F * deltaPhiBin;

  TH1F * deltaPhiDt;

  TH1F * rpcInHitsPerDtseg;
  TH1F * rpcOutHitsPerDtseg;
  TH2F * dtQualityNew;
  
};


L1ITMuonBarrelPlots::L1ITMuonBarrelPlots(const PSet& p)
  : _mbltCollectionInput( p.getParameter<edm::InputTag>("MBLTCollection") ),
    _mbTracksCollectionInput( p.getParameter<edm::InputTag>("MBTracksCollection") ),
    _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") )
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

    /// DT quality stuff
    dtQuality[i] = _fs->make<TH2F>( "dtQuality_"+stname, "dt quality vs rpc match in "+stname, 4, 0, 4, 7, 0, 7);
    dtQuality[i]->GetXaxis()->SetBinLabel( 1, "Unconfirmed");
    dtQuality[i]->GetXaxis()->SetBinLabel( 2, "In+Out");
    dtQuality[i]->GetXaxis()->SetBinLabel( 3, "In");
    dtQuality[i]->GetXaxis()->SetBinLabel( 4, "Out");

    dtQuality[i]->GetYaxis()->SetBinLabel( 1, "LI" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 2, "LO" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 3, "HI" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 4, "HO" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 5, "LL" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 6, "HL" );
    dtQuality[i]->GetYaxis()->SetBinLabel( 7, "HH" );

  }

  deltaPhiDt = _fs->make<TH1F>( "deltaPhiDt", "deltaPhiDt", 400, -0.2, 0.2 );

  deltaPhi = _fs->make<TH1F>( "deltaPhi", "deltaPhi", 400, -0.2, 0.2 );
  deltaPhiBin = _fs->make<TH1F>( "deltaPhiBin", "deltaPhi bins", 3, -0.15, 0.15 );

  deltaEta = _fs->make<TH1F>( "deltaEta", "deltaEta", 40, 0, 0.4 );
  deltaR = _fs->make<TH1F>( "deltaR", "deltaR", 100, 0, 1 );
  rpcInHitsPerDtseg = _fs->make<TH1F>( "rpcInHitsPerDtseg", "number of inner rpc hits matching a dt sector", 20, 0, 20 );
  rpcOutHitsPerDtseg = _fs->make<TH1F>( "rpcOutHitsPerDtseg", "number of outer rpc hits matching a dt sector", 20, 0, 20 );

  deltaPhiR = _fs->make<TH2F>( "deltaPhiR", "deltaR vs deltaPhi", 400, 0, 0.4, 100, 0, 0.5 );

  /// DT quality stuff according to new definition
  dtQualityNew = _fs->make<TH2F>( "dtQualityNew", "dt quality vs station", 4, 1, 5, 7, 0, 7);
  dtQualityNew->GetXaxis()->SetBinLabel( 1, "MB1");
  dtQualityNew->GetXaxis()->SetBinLabel( 2, "MB2");
  dtQualityNew->GetXaxis()->SetBinLabel( 3, "MB3");
  dtQualityNew->GetXaxis()->SetBinLabel( 4, "MB4");

  dtQualityNew->GetYaxis()->SetBinLabel( 1, "HI" );
  dtQualityNew->GetYaxis()->SetBinLabel( 2, "HO" );
  dtQualityNew->GetYaxis()->SetBinLabel( 3, "HI+RPC" );
  dtQualityNew->GetYaxis()->SetBinLabel( 4, "HO+RPC" );
  dtQualityNew->GetYaxis()->SetBinLabel( 5, "#splitline{(HI+HO)}{+RPC@bx0}" );
  dtQualityNew->GetYaxis()->SetBinLabel( 6, "(LL || HL)" );
  dtQualityNew->GetYaxis()->SetBinLabel( 7, "HH" );

}
  
L1ITMuonBarrelPlots::~L1ITMuonBarrelPlots()
{
}




void L1ITMuonBarrelPlots::endJob()
{

  for ( int i = 0; i < 4; ++i ) {

    /// per station quality plots
    for ( int q = 1; q < 8; ++q ) {
      double quality = 0;
      for ( int m = 1; m < 5; ++m ) {
	quality += dtQuality[i]->GetBinContent( m, q );
      }
      for ( int m = 1; m < 5; ++m ) {
	if ( quality )
	  dtQuality[i]->SetBinContent( m, q, dtQuality[i]->GetBinContent( m, q ) / quality );
	else
	  dtQuality[i]->SetBinContent( m, q, 0 );
      }
    }

    /// new quality plot
    int st = i + 1;
    double qualityNew = 0;
    for ( int q = 1; q < 8; ++q ) {
      qualityNew += dtQualityNew->GetBinContent( st, q );
    }
    for ( int q = 1; q < 8; ++q ) {
      if ( qualityNew )
	dtQualityNew->SetBinContent( st, q, dtQualityNew->GetBinContent( st, q ) / qualityNew );
      else
	dtQualityNew->SetBinContent( st, q, 0 );
    }
  }

}


// double phiBending()
// {
// phib_DT-RPC = (x_RPC - x_DT) / (y_RPC - y_DT) and resol_phib-RPC = sqrt(resol_x_RPC^2 + resol_x_DT^2) 


// }


void L1ITMuonBarrelPlots::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBLTContainer> mbltContainer;
  iEvent.getByLabel( _mbltCollectionInput, mbltContainer );

  L1ITMu::MBLTContainer::const_iterator st = mbltContainer->begin();
  L1ITMu::MBLTContainer::const_iterator stend = mbltContainer->end();

  for ( ; st != stend; ++st ) {

    const L1ITMu::MBLTCollection & mbltStation = st->second;

    /// useful index
    int station = mbltStation.station();
    int index = station - 1;
    int wheel = mbltStation.wheel();
    int sector = mbltStation.sector();

    if ( index < 0 || index > 3 )
      throw cms::Exception("Invalid Station")
	<< "wrong station number " << station << std::endl;

    /// size for primitives vectors
    size_t dtListSize = mbltStation.getDtSegments().size();
    size_t rpcInListSize = mbltStation.getRpcInner().size();
    size_t rpcOutListSize = mbltStation.getRpcOuter().size();

    /// fill general distribution plots
    if ( dtListSize ) dtDist[index]->Fill( wheel, sector );
    if ( rpcInListSize ) rpcInDist[index]->Fill( wheel, sector );
    if ( rpcOutListSize ) rpcOutDist[index]->Fill( wheel, sector );

    /// get dt to rpc associations
    size_t dtInTime = 0;
    for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {

      const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);
      int dtquality = dt->getDTData().qualityCode;

      /// skip theta segments
      if ( dtquality < 0 ) continue;

      double eta = dt->getCMSGlobalEta();
      double phi = dt->getCMSGlobalPhi();
      int dtbx = dt->getBX();

      /// rpc associated hits collections
      L1ITMu::TriggerPrimitiveList rpcInMatch = mbltStation.getRpcInAssociatedStubs( iDt );
      L1ITMu::TriggerPrimitiveList rpcOutMatch = mbltStation.getRpcOutAssociatedStubs( iDt );
      size_t rpcInMatchSize = rpcInMatch.size();
      size_t rpcOutMatchSize = rpcOutMatch.size();


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

      /// let's keep only bx=0
      if ( dtbx ) continue;
      ++dtInTime;

      // delta phi
      for ( size_t jDt = iDt+1; jDt < dtListSize; ++jDt ) {
	double phi2 = mbltStation.getDtSegments().at(jDt)->getCMSGlobalPhi();
	double deltaPhiSt = reco::deltaPhi( phi, phi2 );
	deltaPhiDt->Fill( deltaPhiSt );
// 	if ( deltaPhiSt < 0.05 ) {
// 	}
      }


      /// count matching
      rpcInHitsPerDtseg->Fill( rpcInMatchSize );
      rpcOutHitsPerDtseg->Fill( rpcOutMatchSize );

      if ( rpcInMatchSize && rpcOutMatchSize ) {
	confirmed[index]->Fill( 1 );
	dtQuality[index]->Fill( 1, dtquality );
      } else if ( rpcInMatchSize ) {
	confirmed[index]->Fill( 2 );
	dtQuality[index]->Fill( 2, dtquality );
      } else if ( rpcOutMatchSize ) {
	confirmed[index]->Fill( 3 );
	dtQuality[index]->Fill( 3, dtquality );
      } else {
	confirmed[index]->Fill( 0 );
	dtQuality[index]->Fill( 0., dtquality );
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

      if ( rpcInMatchSize ) {
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcInMatch.front()->getCMSGlobalPhi() );
	if ( rpcDeltaPhi < -0.05 ) deltaPhiBin->Fill( -0.1 );
	else if ( rpcDeltaPhi > 0.05 ) deltaPhiBin->Fill( 0.1 );
	else deltaPhiBin->Fill( rpcDeltaPhi );
      }

      if ( rpcOutMatchSize ) {
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcOutMatch.front()->getCMSGlobalPhi() );
	if ( rpcDeltaPhi < -0.05 ) deltaPhiBin->Fill( -0.1 );
	else if ( rpcDeltaPhi > 0.05 ) deltaPhiBin->Fill( 0.1 );
	else deltaPhiBin->Fill( rpcDeltaPhi );
      }


    }

    if ( !dtInTime && ( rpcInListSize || rpcOutListSize ) ) {
      confirmed[index]->Fill( -1 );
      continue;
    }

  }


  /// New primitives loop
  edm::Handle<L1MuDTChambPhContainer> phiChambContainer;
  iEvent.getByLabel( _l1itDtPhiChInput, phiChambContainer );
  std::vector<L1MuDTChambPhDigi>* phTrigs = phiChambContainer->getContainer();

  std::vector<L1MuDTChambPhDigi>::const_iterator iph  = phTrigs->begin();
  std::vector<L1MuDTChambPhDigi>::const_iterator iphe = phTrigs->end();

  for(; iph !=iphe ; ++iph) {
    if ( !iph->bxNum() ) dtQualityNew->Fill( iph->stNum(), iph->code() );
  }
//     iph->phi();
//     iph->phiB();
//     iph->code();

  /// JP : loop over MBtracks
  edm::Handle<L1ITMu::MBTrackCollection> mbtrContainer;
  iEvent.getByLabel( _mbTracksCollectionInput, mbtrContainer);

  L1ITMu::MBTrackCollection::const_iterator tr = mbtrContainer->begin();
  L1ITMu::MBTrackCollection::const_iterator trend = mbtrContainer->end();

  std::cout << "Getting the MBTrackCollection (size " << mbtrContainer->size() << ")\n";
  
  for ( ; tr != trend; ++tr ) {
    
    const L1ITMu::MBTrack & mbtrack = *tr;
    
    std::cout << " - MBTrack :\n";
    
    /// loop over MBCollection
    const L1ITMu::MBLTContainer & mbltContainer = mbtrack.getStubs();

    L1ITMu::MBLTContainer::const_iterator st = mbltContainer.begin();
    L1ITMu::MBLTContainer::const_iterator stend = mbltContainer.end();

    std::cout << "  - Getting the MBLTCollection (size " << mbltContainer.size() << ") from the getStubs() method\n";
    
    for ( ; st != stend; ++st ) {
      
      const L1ITMu::MBLTCollection & mbltStation = st->second;
    }
    
    /// loop over GMTout (L1MuGMTExtendedCand)
    const std::vector<L1MuGMTExtendedCand> l1gmtextcand = mbtrack.getAssociatedGMTout(); 

    std::vector<L1MuGMTExtendedCand>::const_iterator igmtout  = l1gmtextcand.begin();
    std::vector<L1MuGMTExtendedCand>::const_iterator igmtoutend = l1gmtextcand.end();

    std::cout << "  - Getting the vector<L1MuGMTExtendedCand> (size " << l1gmtextcand.size() << ") from the getAssociatedGMTout() method\n";
    
    for ( ; igmtout != igmtoutend; ++igmtout ) {
      
      const L1MuGMTExtendedCand & gmtout = *igmtout;

      std::cout << "   - gmtOutput:\n";
      
      std::cout << "    - Quality = " << gmtout.quality() << std::endl;
      std::cout << "    - Phi     = " << gmtout.phiValue() << std::endl;
      std::cout << "    - Eta     = " << gmtout.etaValue() << std::endl;
      std::cout << "    - Pt      = " << gmtout.ptValue() << std::endl;
      std::cout << "    - Bx      = " << gmtout.bx() << std::endl;
      
    }

    /// loop over GMTin (L1MuRegionalCand)
    const std::vector<L1MuRegionalCand> l1muregcand = mbtrack.getAssociatedGMTin(); 
    
    std::vector<L1MuRegionalCand>::const_iterator igmtin  = l1muregcand.begin();
    std::vector<L1MuRegionalCand>::const_iterator igmtinend = l1muregcand.end();

    std::cout << "  - Getting the vector<L1MuRegionalCand> (size " << l1muregcand.size() << ") from the getAssociatedGMTin() method\n";

    for ( ; igmtin != igmtinend; ++igmtin ) {
      
      const L1MuRegionalCand & gmtin = *igmtin;
      
      std::cout << "   - gmtInput:\n";
      
      std::cout << "    - Quality = " << gmtin.quality() << std::endl;
      std::cout << "    - Phi     = " << gmtin.phiValue() << std::endl;
      std::cout << "    - Bx      = " << gmtin.bx() << std::endl;
      
    }

  }
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelPlots);
