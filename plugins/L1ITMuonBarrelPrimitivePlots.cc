// 
// Class: L1ITMuonBarrelPrimitivePlots
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

class L1ITMuonBarrelPrimitivePlots : public edm::EDAnalyzer {  
  
public:
  L1ITMuonBarrelPrimitivePlots( const PSet& );
  ~L1ITMuonBarrelPrimitivePlots();

  void analyze( const edm::Event&, const edm::EventSetup& );  
private:
  void endJob();
  void beginJob();
  
  edm::InputTag _l1itDtPhiChInput;
  edm::InputTag _l1itDtPhiChInputNew;
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
};


L1ITMuonBarrelPrimitivePlots::L1ITMuonBarrelPrimitivePlots(const PSet& p)
  : _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") ),
    _l1itDtPhiChInputNew( p.getParameter<edm::InputTag>("L1ITDTChambPhContainerNew") )
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

}
  
L1ITMuonBarrelPrimitivePlots::~L1ITMuonBarrelPrimitivePlots()
{
}




void L1ITMuonBarrelPrimitivePlots::beginJob()
{
  
}



void L1ITMuonBarrelPrimitivePlots::endJob()
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


void L1ITMuonBarrelPrimitivePlots::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{

  /// Old primitives loop


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
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelPrimitivePlots);
