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

#include "DataFormats/Common/interface/RefToBase.h"

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
  TH1F * _confirmed[4];
  TH1F * _timingConf[4];
  TH1F * _timingConfIn[4];
  TH1F * _timingConfOut[4];

  TH2F * _dtDist[4];
  TH2F * _rpcInDist[4];
  TH2F * _rpcOutDist[4];
  TH2F * _dtQuality[4];

  TH1F * _deltaPhi;
  TH1F * _deltaEta;
  TH1F * _deltaR;
  TH2F * _deltaPhiR;
  TH1F * _deltaPhiBin;

  TH1F * _deltaPhiDt;

  TH1F * _rpcInHitsPerDtseg;
  TH1F * _rpcOutHitsPerDtseg;
  TH2F * _dtQualityNew;

  TH1F * _deltaPhiDttfGmtOut;
  TH1F * _deltaPhiDttfGmtOutSameBx;
//   TH2F * _dtQualityNew;
  
};


L1ITMuonBarrelPlots::L1ITMuonBarrelPlots(const PSet& p)
  : _mbltCollectionInput( p.getParameter<edm::InputTag>("MBLTCollection") ),
    _mbTracksCollectionInput( p.getParameter<edm::InputTag>("MBTracksCollection") ),
    _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") )
{

  for ( int i = 0; i < 4; ++i ) {
    TString stname = Form( "st%d", i+1 );
    _confirmed[i] = _fs->make<TH1F>( "confirmed_"+stname, "confirmed in "+stname, 5, -1, 4 );

    _confirmed[i]->GetXaxis()->SetBinLabel( 1, "RPC Only");
    _confirmed[i]->GetXaxis()->SetBinLabel( 2, "Un_confirmed");
    _confirmed[i]->GetXaxis()->SetBinLabel( 3, "In+Out");
    _confirmed[i]->GetXaxis()->SetBinLabel( 4, "In");
    _confirmed[i]->GetXaxis()->SetBinLabel( 5, "Out");

    _timingConf[i] = _fs->make<TH1F>( "timing_"+stname, "dt vs rpc timing in "+stname, 4, -2, 2 );
    _timingConfIn[i] = _fs->make<TH1F>( "timingIn_"+stname, "dt vs rpcIn timing in "+stname, 3, -1, 2 );
    _timingConfOut[i] = _fs->make<TH1F>( "timingOut_"+stname, "dt vs rpcOut timing in "+stname, 3, -1, 2 );

    _dtDist[i] = _fs->make<TH2F>( "dtDist_"+stname, "dt primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    _rpcInDist[i] = _fs->make<TH2F>( "rpcInDist_"+stname, "rpc inner primitives in "+stname, 5, -2, 3, 12, 1, 13 );
    _rpcOutDist[i] = _fs->make<TH2F>( "rpcOutDist_"+stname, "rpc outer primitives in "+stname, 5, -2, 3, 12, 1, 13 );

    /// DT quality stuff
    _dtQuality[i] = _fs->make<TH2F>( "dtQuality_"+stname, "dt quality vs rpc match in "+stname, 4, 0, 4, 7, 0, 7);
    _dtQuality[i]->GetXaxis()->SetBinLabel( 1, "Un_confirmed");
    _dtQuality[i]->GetXaxis()->SetBinLabel( 2, "In+Out");
    _dtQuality[i]->GetXaxis()->SetBinLabel( 3, "In");
    _dtQuality[i]->GetXaxis()->SetBinLabel( 4, "Out");

    _dtQuality[i]->GetYaxis()->SetBinLabel( 1, "LI" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 2, "LO" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 3, "HI" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 4, "HO" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 5, "LL" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 6, "HL" );
    _dtQuality[i]->GetYaxis()->SetBinLabel( 7, "HH" );

  }

  _deltaPhiDt = _fs->make<TH1F>( "deltaPhiDt", "deltaPhiDt", 400, -0.2, 0.2 );

  _deltaPhi = _fs->make<TH1F>( "deltaPhi", "deltaPhi", 400, -0.2, 0.2 );
  _deltaPhiBin = _fs->make<TH1F>( "deltaPhiBin", "deltaPhi bins", 3, -0.15, 0.15 );

  _deltaEta = _fs->make<TH1F>( "deltaEta", "deltaEta", 40, 0, 0.4 );
  _deltaR = _fs->make<TH1F>( "deltaR", "deltaR", 100, 0, 1 );
  _rpcInHitsPerDtseg = _fs->make<TH1F>( "rpcInHitsPerDtseg", "number of inner rpc hits matching a dt sector", 20, 0, 20 );
  _rpcOutHitsPerDtseg = _fs->make<TH1F>( "rpcOutHitsPerDtseg", "number of outer rpc hits matching a dt sector", 20, 0, 20 );

  _deltaPhiR = _fs->make<TH2F>( "deltaPhiR", "deltaR vs deltaPhi", 400, 0, 0.4, 100, 0, 0.5 );

  /// DT quality stuff according to new definition
  _dtQualityNew = _fs->make<TH2F>( "dtQualityNew", "dt quality vs station", 4, 1, 5, 7, 0, 7);
  _dtQualityNew->GetXaxis()->SetBinLabel( 1, "MB1");
  _dtQualityNew->GetXaxis()->SetBinLabel( 2, "MB2");
  _dtQualityNew->GetXaxis()->SetBinLabel( 3, "MB3");
  _dtQualityNew->GetXaxis()->SetBinLabel( 4, "MB4");
  _dtQualityNew->GetYaxis()->SetBinLabel( 1, "HI" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 2, "HO" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 3, "HI+RPC" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 4, "HO+RPC" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 5, "#splitline{(HI+HO)}{+RPC@bx0}" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 6, "(LL || HL)" );
  _dtQualityNew->GetYaxis()->SetBinLabel( 7, "HH" );

  _deltaPhiDttfGmtOut       = _fs->make<TH1F>( "deltaPhiDttfGmtOut", "deltaPhiDttfGmtOut", 400, -0.2, 0.2 );
  _deltaPhiDttfGmtOutSameBx = _fs->make<TH1F>( "deltaPhiDttfGmtOutSameBx", "deltaPhiDttfGmtOutSameBx", 400, -0.2, 0.2 );
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
	quality += _dtQuality[i]->GetBinContent( m, q );
      }
      for ( int m = 1; m < 5; ++m ) {
	if ( quality )
	  _dtQuality[i]->SetBinContent( m, q, _dtQuality[i]->GetBinContent( m, q ) / quality );
	else
	  _dtQuality[i]->SetBinContent( m, q, 0 );
      }
    }

    /// new quality plot
    int st = i + 1;
    double qualityNew = 0;
    for ( int q = 1; q < 8; ++q ) {
      qualityNew += _dtQualityNew->GetBinContent( st, q );
    }
    for ( int q = 1; q < 8; ++q ) {
      if ( qualityNew )
	_dtQualityNew->SetBinContent( st, q, _dtQualityNew->GetBinContent( st, q ) / qualityNew );
      else
	_dtQualityNew->SetBinContent( st, q, 0 );
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
    if ( dtListSize ) _dtDist[index]->Fill( wheel, sector );
    if ( rpcInListSize ) _rpcInDist[index]->Fill( wheel, sector );
    if ( rpcOutListSize ) _rpcOutDist[index]->Fill( wheel, sector );

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
	_timingConf[index]->Fill( bx );
      } else if ( rpcInMatchSize ) {
	int rpcbx = rpcInMatch.front()->getBX();
	int bx = 0;
	if ( dtbx == rpcbx ) bx = 0;
	else if ( dtbx > rpcbx ) bx = 1;
	else if ( dtbx < rpcbx ) bx = -1;
	_timingConfIn[index]->Fill( bx );
      } else if ( rpcOutMatchSize ) {
	int rpcbx = rpcOutMatch.front()->getBX();
	int bx = 0;
	if ( dtbx == rpcbx ) bx = 0;
	else if ( dtbx > rpcbx ) bx = 1;
	else if ( dtbx < rpcbx ) bx = -1;
	_timingConfOut[index]->Fill( bx );
      }

      /// let's keep only bx=0
      if ( dtbx ) continue;
      ++dtInTime;

      // delta phi
      for ( size_t jDt = iDt+1; jDt < dtListSize; ++jDt ) {
	double phi2 = mbltStation.getDtSegments().at(jDt)->getCMSGlobalPhi();
	double deltaPhiSt = reco::deltaPhi( phi, phi2 );
	_deltaPhiDt->Fill( deltaPhiSt );
// 	if ( deltaPhiSt < 0.05 ) {
// 	}
      }


      /// count matching
      _rpcInHitsPerDtseg->Fill( rpcInMatchSize );
      _rpcOutHitsPerDtseg->Fill( rpcOutMatchSize );

      if ( rpcInMatchSize && rpcOutMatchSize ) {
	_confirmed[index]->Fill( 1 );
	_dtQuality[index]->Fill( 1, dtquality );
      } else if ( rpcInMatchSize ) {
	_confirmed[index]->Fill( 2 );
	_dtQuality[index]->Fill( 2, dtquality );
      } else if ( rpcOutMatchSize ) {
	_confirmed[index]->Fill( 3 );
	_dtQuality[index]->Fill( 3, dtquality );
      } else {
	_confirmed[index]->Fill( 0 );
	_dtQuality[index]->Fill( 0., dtquality );
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
	_deltaEta->Fill( fabs( eta - rpcEta ) );
	_deltaPhi->Fill( rpcDeltaPhi );
	_deltaR->Fill( dR );
	_deltaPhiR->Fill( fabs( rpcDeltaPhi ), dR );
      }

      TriggerPrimitiveList::const_iterator rpcOutIt = rpcOutMatch.begin();
      for ( ; rpcOutIt != rpcOutEnd; ++rpcOutIt ) {
	double rpcPhi = (*rpcOutIt)->getCMSGlobalPhi();
	double rpcEta = (*rpcOutIt)->getCMSGlobalEta();
	double dR = reco::deltaR( eta, phi, rpcEta, rpcPhi);
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcPhi );
	_deltaEta->Fill( fabs( eta - rpcEta ) );
	_deltaPhi->Fill( rpcDeltaPhi );
	_deltaR->Fill( dR );
	_deltaPhiR->Fill( fabs( rpcDeltaPhi ), dR );
      }

      if ( rpcInMatchSize ) {
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcInMatch.front()->getCMSGlobalPhi() );
	if ( rpcDeltaPhi < -0.05 ) _deltaPhiBin->Fill( -0.1 );
	else if ( rpcDeltaPhi > 0.05 ) _deltaPhiBin->Fill( 0.1 );
	else _deltaPhiBin->Fill( rpcDeltaPhi );
      }

      if ( rpcOutMatchSize ) {
	double rpcDeltaPhi = reco::deltaPhi( phi, rpcOutMatch.front()->getCMSGlobalPhi() );
	if ( rpcDeltaPhi < -0.05 ) _deltaPhiBin->Fill( -0.1 );
	else if ( rpcDeltaPhi > 0.05 ) _deltaPhiBin->Fill( 0.1 );
	else _deltaPhiBin->Fill( rpcDeltaPhi );
      }


    }

    if ( !dtInTime && ( rpcInListSize || rpcOutListSize ) ) {
      _confirmed[index]->Fill( -1 );
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
    if ( !iph->bxNum() ) _dtQualityNew->Fill( iph->stNum(), iph->code() );
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

    std::cout << "    - Station = " << mbltStation.station() << std::endl;
    std::cout << "    - Wheel   = " << mbltStation.wheel() << std::endl;
    std::cout << "    - Sector  = " << mbltStation.sector() - 1 << std::endl;
    
    /// useful index
    int station = mbltStation.station();
    int index = station - 1;
    int wheel = mbltStation.wheel();
    int sector = mbltStation.sector() - 1 ;

    if ( index < 0 || index > 3 )
      throw cms::Exception("Invalid Station")
        << "wrong station number " << station << std::endl;

    /// size for primitives vectors
    size_t dtListSize = mbltStation.getDtSegments().size();

    std::cout << "     - dtListSize = " << mbltStation.getDtSegments().size() << std::endl;

  }


    /// get dttf    
    const L1MuDTTrackCand& dttf = mbtrack.parent();
    std::cout << "  - Getting the L1MuDTTrackCand DTTF from the parent() method\n";
    
    int dttf_bx      =  dttf.bx();
    int dttf_qual    =  dttf.quality();
    int phi_local = dttf.phi_packed();
    if ( phi_local > 15 ) phi_local -= 32;
    double dttf_phi_global = (phi_local*(M_PI/72.))+((M_PI/6.)*dttf.scNum());
    if(dttf_phi_global < 0) dttf_phi_global+=2*M_PI;
    if(dttf_phi_global > 2*M_PI) dttf_phi_global-=2*M_PI;
    
    std::cout << "    - Quality = " << dttf.quality() << std::endl;
    std::cout << "    - Bx      = " << dttf.bx() << std::endl;
    std::cout << "    - Phi     = " << dttf_phi_global << std::endl;
    std::cout << "    - WhNum   = " << dttf.whNum() << std::endl;
    std::cout << "    - ScNum   = " << dttf.scNum() << std::endl;

    
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
      
      double dttfGmtDeltaPhi = reco::deltaPhi( dttf_phi_global, gmtout.phiValue() );
      
      _deltaPhiDttfGmtOut->Fill(dttfGmtDeltaPhi);
      
      if ( dttf_bx == gmtout.bx() ) { 
        _deltaPhiDttfGmtOutSameBx->Fill(dttfGmtDeltaPhi);
      }
      
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
