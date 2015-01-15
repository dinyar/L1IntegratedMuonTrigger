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

  TH1F * _uncorrelated;
  TH1F * _uncorrelatedNorm;
  // TH1F * _deltaPhi;

  TH2F * _bendingAngle;
  TH2F * _radialAngle;

  TH1F * _bendingAngleDiff;
  TH1F * _radialAngleDiff;

  TH2F * _bendingAnglePhys;
  TH2F * _radialAnglePhys;

  TH1F * _bendingAnglePhysDiff;
  TH1F * _radialAnglePhysDiff;

  TH2F * _dtQuality;
  TH2F * _dtQualityNew;
};


L1ITMuonBarrelPrimitivePlots::L1ITMuonBarrelPrimitivePlots(const PSet& p)
  : _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") ),
    _l1itDtPhiChInputNew( p.getParameter<edm::InputTag>("L1ITDTChambPhContainerNew") )
{

  _uncorrelated = _fs->make<TH1F>( "uncorrelated", "uncorrelated", 4, 1, 5 );
  _uncorrelated->GetXaxis()->SetBinLabel( 1, "total");
  _uncorrelated->GetXaxis()->SetBinLabel( 2, "final");
  _uncorrelated->GetXaxis()->SetBinLabel( 3, "final+rpc");
  _uncorrelated->GetXaxis()->SetBinLabel( 4, "rematched");
  _uncorrelatedNorm = _fs->make<TH1F>( "uncorrelatedNorm", "uncorrelatedNorm", 4, 1, 5 );
  _uncorrelatedNorm->GetXaxis()->SetBinLabel( 1, "total");
  _uncorrelatedNorm->GetXaxis()->SetBinLabel( 2, "final");
  _uncorrelatedNorm->GetXaxis()->SetBinLabel( 3, "final+rpc");
  _uncorrelatedNorm->GetXaxis()->SetBinLabel( 4, "rematched");


  _bendingAngle = _fs->make<TH2F>( "bendingAngle", "bendingAngle", 104, -520, 520, 104, -520, 520 );
  _radialAngle = _fs->make<TH2F>( "radialAngle", "radialAngle", 360, -1800, 1800, 360, -1800, 1800 );

  _bendingAngleDiff = _fs->make<TH1F>( "bendingAngleDiff", "bendingAngleDiff", 500, -500, 500 );
  _radialAngleDiff = _fs->make<TH1F>( "radialAngleDiff", "radialAngleDiff", 200, -100, 100 );


  _bendingAnglePhys = _fs->make<TH2F>( "bendingAnglePhys", "bendingAnglePhys", 63, -3.15, 3.15, 63, -3.15, 3.15 );
  _radialAnglePhys = _fs->make<TH2F>( "radialAnglePhys", "radialAnglePhys", 63, -3.15, 3.15, 63, -3.15, 3.15 );

  _bendingAnglePhysDiff = _fs->make<TH1F>( "bendingAnglePhysDiff", "bendingAnglePhysDiff", 100, -1, 1 );
  _radialAnglePhysDiff = _fs->make<TH1F>( "radialAnglePhysDiff", "radialAnglePhysDiff", 200, -2, 2 );


  /// DT quality stuff
  _dtQuality = _fs->make<TH2F>( "dtQuality", "dt quality vs station", 4, 1, 5, 7, 0, 7);
  _dtQuality->GetXaxis()->SetBinLabel( 1, "MB1");
  _dtQuality->GetXaxis()->SetBinLabel( 2, "MB2");
  _dtQuality->GetXaxis()->SetBinLabel( 3, "MB3");
  _dtQuality->GetXaxis()->SetBinLabel( 4, "MB4");
  _dtQuality->GetYaxis()->SetBinLabel( 1, "LI" );
  _dtQuality->GetYaxis()->SetBinLabel( 2, "LO" );
  _dtQuality->GetYaxis()->SetBinLabel( 3, "HI" );
  _dtQuality->GetYaxis()->SetBinLabel( 4, "HO" );
  _dtQuality->GetYaxis()->SetBinLabel( 5, "LL" );
  _dtQuality->GetYaxis()->SetBinLabel( 6, "HL" );
  _dtQuality->GetYaxis()->SetBinLabel( 7, "HH" );

  /// DT quality stuff according to new definition
  _dtQualityNew = _fs->make<TH2F>( "dtQualityNew", "New dt quality vs station", 4, 1, 5, 7, 0, 7);
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

  for ( int st = 1; st < 5; ++st ) {

    /// new quality plot
    double quality = 0;
    double qualityNew = 0;
    for ( int q = 1; q < 8; ++q ) {
      quality += _dtQuality->GetBinContent( st, q );
      qualityNew += _dtQualityNew->GetBinContent( st, q );
    }
    for ( int q = 1; q < 8; ++q ) {
      if ( quality )
	_dtQuality->SetBinContent( st, q, _dtQuality->GetBinContent( st, q ) / quality );
      else
	_dtQuality->SetBinContent( st, q, 0 );

      if ( qualityNew )
	_dtQualityNew->SetBinContent( st, q, _dtQualityNew->GetBinContent( st, q ) / qualityNew );
      else
	_dtQualityNew->SetBinContent( st, q, 0 );
    }
  }

  double totUnc = _uncorrelated->GetBinContent( 1 );
  for ( int i = 1; i < 5; ++i ) {
    _uncorrelatedNorm->SetBinContent( i, _uncorrelated->GetBinContent( i ) / totUnc );
  }

}



void L1ITMuonBarrelPrimitivePlots::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{

  std::vector<double> oldBendingAngle;
  std::vector<double> oldRadialAngle;
  std::vector<double> newBendingAngle;
  std::vector<double> newRadialAngle;

  /// Old primitives loop
  edm::Handle<L1MuDTChambPhContainer> phiChambContainer;
  iEvent.getByLabel( _l1itDtPhiChInput, phiChambContainer );

  const std::vector<L1MuDTChambPhDigi>& phTrigs
    = *( phiChambContainer->getContainer() );

  for ( const auto & iph : phTrigs ) {

    if ( iph.code() == 2 || iph.code() == 3 ) {
      _uncorrelated->Fill( 1 );
      oldBendingAngle.push_back( iph.phiB() );
      oldRadialAngle.push_back( iph.phi() );
    }

    if ( iph.bxNum() ) continue;
    _dtQuality->Fill( iph.stNum(), iph.code() );
  }

  /// New primitives loop
  edm::Handle<L1MuDTChambPhContainer> phiChambContainerNew;
  iEvent.getByLabel( _l1itDtPhiChInputNew, phiChambContainerNew );

  const std::vector<L1MuDTChambPhDigi> & phTrigsNew =
    *( phiChambContainerNew->getContainer() );

  for ( const auto & iphNew : phTrigsNew ) {

    if ( iphNew.code() < 4  ) {

      // any uncorrelated is bin 2
      _uncorrelated->Fill( 2 );

      // uncorrelated with rpc match is bin 3 and has bending recalculation
      if ( iphNew.code() > 1  ) {
	_uncorrelated->Fill( 3 );
	newBendingAngle.push_back( iphNew.phiB() );
      }

      /// recombined uncorrelated is bin 4 and has phi and phiB recalculated
    } else if ( iphNew.code() == 4 ) {
      _uncorrelated->Fill( 4 );
      newRadialAngle.push_back( iphNew.phi() );
      newBendingAngle.push_back( iphNew.phiB() );
    }

    if ( iphNew.bxNum() ) continue;
    _dtQualityNew->Fill( iphNew.stNum(), iphNew.code() );

  }

  /// the following loops make sense only if there is a new dt+dt+rpc or dt+rpc match
  for ( size_t j = 0; j < newBendingAngle.size(); ++j ) {
    for ( size_t i = 0; i < oldBendingAngle.size(); ++i ) {

      double diff = oldBendingAngle.at(i) - newBendingAngle.at(j);
      _bendingAngle->Fill( oldBendingAngle.at(i), newBendingAngle.at(j) );
      _bendingAngleDiff->Fill( diff);

      _bendingAnglePhys->Fill( oldBendingAngle.at(i)/515., newBendingAngle.at(j)/515. );
      _bendingAnglePhysDiff->Fill( diff/515. );
    }
  }

  for ( size_t j = 0; j < newRadialAngle.size(); ++j ) {
    for ( size_t i = 0; i < oldRadialAngle.size(); ++i ) {

      double diff = oldRadialAngle.at(i) - newRadialAngle.at(j);
      _radialAngle->Fill( oldRadialAngle.at(i), newRadialAngle.at(j) );
      _radialAngleDiff->Fill( oldRadialAngle.at(i) - newRadialAngle.at(j) );

      _radialAnglePhys->Fill( oldRadialAngle.at(i)/4096., newRadialAngle.at(j)/4096. );
      _radialAnglePhysDiff->Fill( diff/4096. );
    }
  }


  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelPrimitivePlots);
