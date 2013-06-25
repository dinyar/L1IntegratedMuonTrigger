// 
// Class: L1ITMuonBarrelTrackPlots
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrackFwd.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"

#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

using namespace L1ITMu;

typedef edm::ParameterSet PSet;

class L1ITMuonBarrelTrackPlots : public edm::EDAnalyzer {  
  
public:
  L1ITMuonBarrelTrackPlots( const PSet& );
  ~L1ITMuonBarrelTrackPlots();

  void analyze( const edm::Event&, const edm::EventSetup& );  
private:
  void endJob();
  void beginJob();
  
  edm::InputTag _mbltCollectionInput;
  edm::InputTag _mbTracksCollectionInput;
  edm::InputTag _l1itDtPhiChInput;
  edm::Service<TFileService> _fs;
  
  TH1F * _deltaPhiDttfGmtInAllBx;
  TH1F * _deltaPhiDttfGmtIn;
  TH1F * _deltaPhiDttfGmtOutAllBx;
  TH1F * _deltaPhiDttfGmtOut;
  
//   TH1F * _deltaRDttfGmtOutSameBx;
//   TH1F * _deltaRDttfGmtOut;
  
  TH2F * _dttfQualityVsL1Path;
  
  TH1F * _nGmtIn;
  TH1F * _nGmtOut;
  TH2F * _nGmtInVsGmtOutPerDTTF;
  
  TH2F * _nGmtInVsDttfQual;
  TH2F * _nGmtOutVsDttfQual; 
  
  TH2F * _dtQualityNewVsL1Path;
  TH2F * _dtQualityNewVsL1PathAllBx;
  
  TH2F * _nDttfVsGmtIn; // effect of the sorting
  TH2F * _nDttfVsGmtInAllBx; // effect of the sorting

  TH2F * _nDttfVsGmtOut; // effect of the sorting
  TH2F * _nDttfVsGmtOutAllBx; // effect of the sorting
  
  TH2F * _nGmtInVsGmtOut; // effect of the sorting
  TH2F * _nGmtInVsGmtOutAllBx; // effect of the sorting
  
  TH2F * _positionGmtOutNoGmtIn; // position (eta/phi) of GMT out matched witha DTTF but with NO-GMTin (ARE THEY CSC-like?)
  TH2F * _positionGmtOutNoGmtInAllBx; // position (eta/phi) of GMT out matched witha DTTF but with NO-GMTin (ARE THEY CSC-like?)

  TH2F * _GmtQualityVsL1Path;
  TH2F * _GmtQualityVsEta;
  TH2F * _GmtQualityVsEtaNoGmtIn;

  
//   int nTr,nMoreGMTin,nNoGMTinMoreGMTout,nNoGMTin,xxx,yyy,zzz,aaa;
};


L1ITMuonBarrelTrackPlots::L1ITMuonBarrelTrackPlots(const PSet& p)
  : _mbltCollectionInput( p.getParameter<edm::InputTag>("MBLTCollection") ),
    _mbTracksCollectionInput( p.getParameter<edm::InputTag>("MBTracksCollection") ),
    _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") )
{

  _deltaPhiDttfGmtInAllBx  = _fs->make<TH1F>( "deltaPhiDttfGmtInAllBx", "deltaPhiDttfGmtInAllBx", 200, -0.00001, 0.00001 );
  _deltaPhiDttfGmtIn       = _fs->make<TH1F>( "deltaPhiDttfGmtIn", "deltaPhiDttfGmtIn", 200, -0.00001, 0.00001 );
  _deltaPhiDttfGmtOutAllBx = _fs->make<TH1F>( "deltaPhiDttfGmtOutAllBx", "deltaPhiDttfGmtOutAllBx", 600, -0.3, 0.3 );
  _deltaPhiDttfGmtOut      = _fs->make<TH1F>( "deltaPhiDttfGmtOut", "deltaPhiDttfGmtOut", 600, -0.3, 0.3 );

//   _deltaRDttfGmtOutAllBx = _fs->make<TH1F>( "deltaRDttfGmtOutAllBx", "deltaRDttfGmtOutAllBx", 600, -0.3, 0.3 );
//   _deltaRDttfGmtOut      = _fs->make<TH1F>( "deltaRDttfGmtOut", "deltaRDttfGmtOut", 600, -0.3, 0.3 );

  _dttfQualityVsL1Path = _fs->make<TH2F>( "dttfQual", "DTTF quality", 7, 0, 7, 7, 1, 8);
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 5, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 6, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _dttfQualityVsL1Path->GetXaxis()->SetBinLabel( 7, "#splitline{DTTF + > 1 GMTin}{or > 2 GMTout}");
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 1, "T34" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 4, "T234" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 5, "T134" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _dttfQualityVsL1Path->GetYaxis()->SetBinLabel( 7, "T1234" );
 
  _nGmtIn       = _fs->make<TH1F>( "nGmtIn",  "number of GMTin per DTTF", 4, 0, 4 );
  _nGmtOut      = _fs->make<TH1F>( "nGmtOut", "number of GMTout per DTTF", 4, 0, 4 );
  _nGmtInVsGmtOutPerDTTF      = _fs->make<TH2F>( "nGmtInVsGmtOutPerDTTF", "number of GMTin Vs number of GMTout per DTTF", 4, 0, 4 , 4, 0, 4 );
  _nGmtInVsGmtOutPerDTTF->SetXTitle("n GMT input per DTTF");
  _nGmtInVsGmtOutPerDTTF->SetYTitle("n GMT output per DTTF");
  
  _nGmtInVsDttfQual = _fs->make<TH2F>( "nGmtInVsDttfQual", "number of GMTin per DTTF vs DTTF quality", 4, 0, 4, 7, 1, 8);
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtInVsDttfQual->SetXTitle("n GMT input per DTTF");

  _nGmtOutVsDttfQual = _fs->make<TH2F>( "nGmtOutVsDttfQual", "number of GMTout per DTTF vs DTTF quality", 4, 0, 4, 7, 1, 8);
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtOutVsDttfQual->SetXTitle("n GMT output per DTTF");
  
  _dtQualityNewVsL1Path = _fs->make<TH2F>( "dtQualityNewVsL1Path", "dt quality", 7, 0, 7, 7, 0, 7);
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 5, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 6, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _dtQualityNewVsL1Path->GetXaxis()->SetBinLabel( 7, "#splitline{DTTF + > 1 GMTin}{or > 2 GMTout}");
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 1, "HI" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 2, "HO" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 3, "HI+RPC" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 4, "HO+RPC" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 5, "#splitline{(HI+HO)}{+RPC@bx0}" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 6, "(LL || HL)" );
  _dtQualityNewVsL1Path->GetYaxis()->SetBinLabel( 7, "HH" );
  
  _dtQualityNewVsL1PathAllBx = _fs->make<TH2F>( "dtQualityNewVsL1PathAllBx", "dt quality", 7, 0, 7, 7, 0, 7);
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 5, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 6, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _dtQualityNewVsL1PathAllBx->GetXaxis()->SetBinLabel( 7, "#splitline{DTTF + > 1 GMTin}{or > 2 GMTout}");
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 1, "HI" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 2, "HO" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 3, "HI+RPC" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 4, "HO+RPC" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 5, "#splitline{(HI+HO)}{+RPC@bx0}" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 6, "(LL || HL)" );
  _dtQualityNewVsL1PathAllBx->GetYaxis()->SetBinLabel( 7, "HH" );

  _nDttfVsGmtIn = _fs->make<TH2F>( "nDttfVsGmtIn", "number of DTTF vs number of GMTin", 10, 0, 10, 4, 0, 4);
  _nDttfVsGmtIn->SetXTitle("n DTTF");
  _nDttfVsGmtIn->SetYTitle("n GMTin");

  _nDttfVsGmtInAllBx = _fs->make<TH2F>( "nDttfVsGmtInAllBx", "number of DTTF vs number of GMTin", 10, 0, 10, 4, 0, 4);
  _nDttfVsGmtInAllBx->SetXTitle("n DTTF");
  _nDttfVsGmtInAllBx->SetYTitle("n GMTin");

  _nDttfVsGmtOut = _fs->make<TH2F>( "nDttfVsGmtOut", "number of DTTF vs number of GMTout", 10, 0, 10, 4, 0, 4);
  _nDttfVsGmtOut->SetXTitle("n DTTF");
  _nDttfVsGmtOut->SetYTitle("n GMTout");

  _nDttfVsGmtOutAllBx = _fs->make<TH2F>( "nDttfVsGmtOutAllBx", "number of DTTF vs number of GMTout", 10, 0, 10, 4, 0, 4);
  _nDttfVsGmtOutAllBx->SetXTitle("n DTTF");
  _nDttfVsGmtOutAllBx->SetYTitle("n GMTout");

  _nGmtInVsGmtOut = _fs->make<TH2F>( "nGmtInVsGmtOut", "number of GMTin vs number of GMTout", 10, 0, 10, 4, 0, 4);
  _nGmtInVsGmtOut->SetXTitle("n GMTin");
  _nGmtInVsGmtOut->SetYTitle("n GMTout");

  _nGmtInVsGmtOutAllBx = _fs->make<TH2F>( "nGmtInVsGmtOutAllBx", "number of GMTin vs number of GMTout", 10, 0, 10, 4, 0, 4);
  _nGmtInVsGmtOutAllBx->SetXTitle("n GMTin");
  _nGmtInVsGmtOutAllBx->SetYTitle("n GMTout");
  
  
  _positionGmtOutNoGmtIn = _fs->make<TH2F>( "positionGmtOutNoGmtIn", "position(#eta-#phi) of DTTF-matched GMTout with no GMTin", 30, -1.5, 1.5, 64, 0, 6.4);
  _positionGmtOutNoGmtIn->SetXTitle("GMTout #eta");
  _positionGmtOutNoGmtIn->SetYTitle("GMTout #phi");
  
  _positionGmtOutNoGmtInAllBx = _fs->make<TH2F>( "positionGmtOutNoGmtInAllBx", "position(#eta-#phi) of DTTF-matched GMTout with no GMTin", 30, -1.5, 1.5, 64, 0, 6.4);
  _positionGmtOutNoGmtInAllBx->SetXTitle("GMTout #eta");
  _positionGmtOutNoGmtInAllBx->SetYTitle("GMTout #phi");

  _GmtQualityVsL1Path = _fs->make<TH2F>( "GmtQualityVsL1Path", "GMTout quality", 4, 0, 4, 64, 0, 6.4);
  _GmtQualityVsL1Path->GetXaxis()->SetBinLabel( 1, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _GmtQualityVsL1Path->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _GmtQualityVsL1Path->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _GmtQualityVsL1Path->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtQualityVsL1Path->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );

  _GmtQualityVsEta = _fs->make<TH2F>( "GmtQualityVsEta", "GMTout quality", 30, -1.5, 1.5 ,7, 1, 8);
  _GmtQualityVsEta->SetXTitle("GMTout #eta");
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtQualityVsEta->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );

  _GmtQualityVsEtaNoGmtIn = _fs->make<TH2F>( "GmtQualityVsEtaNoGmtIn", "GMTout quality", 30, -1.5, 1.5 ,7, 1, 8);
  _GmtQualityVsEtaNoGmtIn->SetXTitle("GMTout #eta");
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );
  
}
  
L1ITMuonBarrelTrackPlots::~L1ITMuonBarrelTrackPlots()
{
}




void L1ITMuonBarrelTrackPlots::beginJob()
{
  
}



void L1ITMuonBarrelTrackPlots::endJob()
{


}


void L1ITMuonBarrelTrackPlots::analyze( const edm::Event& iEvent, 
					 const edm::EventSetup& iSetup )
{
  


std::vector<L1MuRegionalCand> vGmtIn;
std::vector<L1MuGMTExtendedCand> vGmtOut;

  /// JP : loop over MBtracks
  edm::Handle<L1ITMu::MBTrackCollection> mbtrContainer;
  iEvent.getByLabel( _mbTracksCollectionInput, mbtrContainer);

  L1ITMu::MBTrackCollection::const_iterator tr = mbtrContainer->begin();
  L1ITMu::MBTrackCollection::const_iterator trend = mbtrContainer->end();

  int nTotGmtIn=0, // total n of GMTin in an event (FIXME -> should be checked against double counting due to the matching procedure...)
      nTotGmtInSameBx=0, // same as before, but for the same bx as the DTTF
      nTotGmtOut=0, // total n of GMTout in an event (FIXME -> should be checked against double counting due to the matching procedure...)
      nTotGmtOutSameBx=0, // same as before, but for the same bx as the DTTF
      nMBTr=0;  // same as n of DTTF!!!

  nMBTr = mbtrContainer->size();
  
  for ( ; tr != trend; ++tr ) {
    
    const L1ITMu::MBTrack & mbtrack = *tr;
    
    int nGmtIn=0,
        nGmtInSameBx=0,
        nGmtOut=0,
        nGmtOutSameBx=0;
    
    /// get dttf    
    const L1MuDTTrackCand& dttf = mbtrack.parent();
    
    int dttf_bx      =  dttf.bx();
    int dttf_qual    =  dttf.quality();
    int phi_local = dttf.phi_packed();
    if ( phi_local > 15 ) phi_local -= 32;
    double dttf_phi_global = (phi_local*(M_PI/72.))+((M_PI/6.)*dttf.scNum());
    if(dttf_phi_global < 0) dttf_phi_global+=2*M_PI;
    if(dttf_phi_global > 2*M_PI) dttf_phi_global-=2*M_PI;
    
    /// loop over GMTin (L1MuRegionalCand)
    const std::vector<L1MuRegionalCand> l1muregcand = mbtrack.getAssociatedGMTin(); 
    std::vector<L1MuRegionalCand>::const_iterator igmtin  = l1muregcand.begin();
    std::vector<L1MuRegionalCand>::const_iterator igmtinend = l1muregcand.end();
    nGmtIn = l1muregcand.size();
    
    for ( ; igmtin != igmtinend; ++igmtin ) {
      
      const L1MuRegionalCand & gmtin = *igmtin;
            
//       ////////////////////////////////
//       std::vector<L1MuRegionalCand>::const_iterator it  = vGmtIn.begin();
//       std::vector<L1MuRegionalCand>::const_iterator itend = vGmtIn.end();
//       for( ; it != itend; ++it){
//         if (*igmtin == *it) std::cout << "HO TROVATO UN GMT input DOPPIONE NELL'EVENTO\n";
//       }
//       ////////////////////////////////
      
      double dph = reco::deltaPhi( dttf_phi_global, gmtin.phiValue() );      
      _deltaPhiDttfGmtInAllBx->Fill(dph);

      if ( dttf_bx == gmtin.bx() ) { 
        
//         vGmtIn.push_back(gmtin);
        
        _deltaPhiDttfGmtIn->Fill(dph);
        nGmtInSameBx++;
      }   
      
    } // end loop over GMTin
    
    /// loop over GMTout (L1MuGMTExtendedCand)
    const std::vector<L1MuGMTExtendedCand> l1gmtextcand = mbtrack.getAssociatedGMTout(); 
    std::vector<L1MuGMTExtendedCand>::const_iterator igmtout  = l1gmtextcand.begin();
    std::vector<L1MuGMTExtendedCand>::const_iterator igmtoutend = l1gmtextcand.end();
    nGmtOut = l1gmtextcand.size();
    
    for ( ; igmtout != igmtoutend; ++igmtout ) {
      
      const L1MuGMTExtendedCand & gmtout = *igmtout;

      ////////////////////////////////
      std::vector<L1MuGMTExtendedCand>::const_iterator it  = vGmtOut.begin();
      std::vector<L1MuGMTExtendedCand>::const_iterator itend = vGmtOut.end();
      for( ; it != itend; ++it){
        if (*igmtout == *it) std::cout << "HO TROVATO UN GMT output DOPPIONE NELL'EVENTO\n";
      }
      ////////////////////////////////
      
      double dph = reco::deltaPhi( dttf_phi_global, gmtout.phiValue() );      
//       double dr  = reco::deltaR( dttf_phi_global, dttf_eta_global, gmtout.phiValue() , gmtout.etaValue() );   
      
      _deltaPhiDttfGmtOutAllBx->Fill(dph);
//       _deltaRDttfGmtOutAllBx->Fill(dr);
      if (!nGmtIn) _positionGmtOutNoGmtInAllBx->Fill(gmtout.etaValue(),gmtout.phiValue());
      
      if ( dttf_bx == gmtout.bx() ) { 
        nGmtOutSameBx++;
        _deltaPhiDttfGmtOut->Fill(dph);
//         _deltaRDttfGmtOut->Fill(dr);
        _GmtQualityVsEta->Fill(gmtout.etaValue(),gmtout.quality());
        if (!nGmtInSameBx) {
          _positionGmtOutNoGmtIn->Fill(gmtout.etaValue(),gmtout.phiValue());
          _GmtQualityVsEtaNoGmtIn->Fill(gmtout.etaValue(),gmtout.quality());          
        }
      }     
      
    } // end loop over GMTout
    
/// --- ///    
    // 0 -> DTTF only
    // 1 -> DTTF + 1 GMTin + 0 GMTout
    // 2 -> DTTF + 1 GMTin + 1 GMTout
    // 3 -> DTTF + 1 GMTin + 2 GMTout
    // 4 -> DTTF + 0 GMTin + 1 GMTout
    // 5 -> DTTF + 0 GMTin + 2 GMTout
    // 6 -> DTTF + > 1 GMTin || > 2 GMTout
    int L1Path, L1PathAll;    
    // same bx as dttf
    if ((nGmtInSameBx < 2) && (nGmtOutSameBx < 3))
    {
      if ( !nGmtInSameBx && nGmtOutSameBx ) L1Path = 3 + nGmtOutSameBx;
      else L1Path = nGmtInSameBx + nGmtOutSameBx;
    }
    else L1Path = 6;
    // all bx    
    if ((nGmtIn < 2) && (nGmtOut < 3))
    {
      if ( !nGmtIn && nGmtOut ) L1PathAll = 3 + nGmtOut;
      else L1PathAll = nGmtIn + nGmtOut;
    }
    else L1PathAll = 6;
/// --- ///    
    
    
    /// loop over MBCollection
    const L1ITMu::MBLTVectorRef & mbltContainer = mbtrack.getStubs();

    L1ITMu::MBLTVectorRef::const_iterator st = mbltContainer.begin();
    L1ITMu::MBLTVectorRef::const_iterator stend = mbltContainer.end();
    
    for ( ; st != stend; ++st ) {

      const L1ITMu::MBLTCollection & mbltStation = (*st)->second;
    
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
      size_t rpcInListSize = mbltStation.getRpcInner().size();
      size_t rpcOutListSize = mbltStation.getRpcOuter().size();
      
      for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {
        
        const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);

        int dtquality = dt->getDTData().qualityCode;
        
        if (dtquality < 0) continue; // skipping theta chambers
            
//         _dtQualityNewVsL1PathAllBx->Fill(L1PathAll,dtquality);
              
        if ( dttf_bx == dt->getDTData().bx ) { 
          _dtQualityNewVsL1Path->Fill(L1Path,dtquality);
        }   
        
      }
      
    } // end loop over MBStations 
    
/// FILL HISTOS    
    /// ALL BX
    nGmtIn  = l1muregcand.size();
    nGmtOut = l1gmtextcand.size(); 
    
    /// SAME BX as DTTF   

    _nGmtIn->Fill(nGmtInSameBx);
    _nGmtOut->Fill(nGmtOutSameBx);
    _nGmtInVsGmtOutPerDTTF->Fill(nGmtInSameBx,nGmtOutSameBx);
    _nGmtInVsDttfQual->Fill(nGmtInSameBx,dttf_qual);
    _nGmtOutVsDttfQual->Fill(nGmtOutSameBx,dttf_qual);
    _dttfQualityVsL1Path->Fill(L1Path,dttf_qual);
   
    nTotGmtIn       += nGmtIn;
    nTotGmtInSameBx += nGmtInSameBx;    
    nTotGmtOut       += nGmtOut;
    nTotGmtOutSameBx += nGmtOutSameBx;    
  }
  
  _nDttfVsGmtIn     ->Fill(nMBTr,nTotGmtInSameBx);
  _nDttfVsGmtInAllBx->Fill(nMBTr,nTotGmtIn);

  _nDttfVsGmtOut     ->Fill(nMBTr,nTotGmtOutSameBx);
  _nDttfVsGmtOutAllBx->Fill(nMBTr,nTotGmtOut);
  
  _nGmtInVsGmtOut     ->Fill(nTotGmtInSameBx,nTotGmtOutSameBx);
  _nGmtInVsGmtOutAllBx->Fill(nTotGmtIn,nTotGmtOut);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelTrackPlots);
