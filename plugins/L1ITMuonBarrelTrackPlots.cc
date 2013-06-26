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
  
  TH1F * _deltaPhiDttfGmtIn;
  TH1F * _deltaPhiDttfGmtOut;
  
//   TH1F * _deltaRDttfGmtOut;
  
  TH2F * _dttfQualityVsCategory;
  
  TH1F * _nGmtIn;
  TH1F * _nGmtOut;
  TH2F * _nGmtInVsGmtOutPerDTTF;
  
  TH2F * _nGmtInVsDttfQual;
  TH2F * _nGmtOutVsDttfQual; 
  
  TH2F * _dtQualityNewVsCategory;
  
  TH2F * _nDttfVsGmtIn; // effect of the sorting

  TH2F * _nDttfVsGmtOut; // effect of the sorting
  
  TH2F * _nGmtInVsGmtOut; // effect of the sorting
  
  TH2F * _positionGmtOutNoGmtIn; // position (eta/phi) of GMT out matched witha DTTF but with NO-GMTin (ARE THEY CSC-like?)

  TH2F * _GmtQualityVsEtaNoGmtIn;
  
  TH2F * _GmtOutQualCat3;
  TH2F * _DttfQualVsGmtOutQualCat3;                
  TH2F * _DttfQualVsGmtOutHiQualCat3;                
  TH2F * _DttfQualVsGmtOutLoQualCat3;                

  TH1F * _GmtOutQualCat4;
  TH2F * _DttfQualVsGmtOutQualCat4;                

  TH2F * _GmtOutQualCat5;
  TH2F * _DttfQualVsGmtOutQualCat5;                
  TH2F * _DttfQualVsGmtOutHiQualCat5;                
  TH2F * _DttfQualVsGmtOutLoQualCat5;                
};


L1ITMuonBarrelTrackPlots::L1ITMuonBarrelTrackPlots(const PSet& p)
  : _mbltCollectionInput( p.getParameter<edm::InputTag>("MBLTCollection") ),
    _mbTracksCollectionInput( p.getParameter<edm::InputTag>("MBTracksCollection") ),
    _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") )
{

  _deltaPhiDttfGmtIn       = _fs->make<TH1F>( "deltaPhiDttfGmtIn", "deltaPhiDttfGmtIn", 200, -0.00001, 0.00001 );
  _deltaPhiDttfGmtOut      = _fs->make<TH1F>( "deltaPhiDttfGmtOut", "deltaPhiDttfGmtOut", 600, -0.3, 0.3 );

//   _deltaRDttfGmtOut      = _fs->make<TH1F>( "deltaRDttfGmtOut", "deltaRDttfGmtOut", 600, -0.3, 0.3 );

  _dttfQualityVsCategory = _fs->make<TH2F>( "dttfQual", "DTTF quality", 7, 0, 7, 7, 1, 8);
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 5, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 6, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 7, "#splitline{DTTF + > 1 GMTin}{or > 2 GMTout}");
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 1, "T34" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 4, "T234" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 5, "T134" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 7, "T1234" );
 
  _nGmtIn       = _fs->make<TH1F>( "nGmtIn",  "number of GMTin per DTTF", 6, 0, 6 );
  _nGmtOut      = _fs->make<TH1F>( "nGmtOut", "number of GMTout per DTTF", 6, 0, 6 );
  _nGmtInVsGmtOutPerDTTF      = _fs->make<TH2F>( "nGmtInVsGmtOutPerDTTF", "number of GMTin Vs number of GMTout per DTTF", 6, 0, 6 , 6, 0, 6 );
  _nGmtInVsGmtOutPerDTTF->SetXTitle("n GMT input per DTTF");
  _nGmtInVsGmtOutPerDTTF->SetYTitle("n GMT output per DTTF");
  
  _nGmtInVsDttfQual = _fs->make<TH2F>( "nGmtInVsDttfQual", "number of GMTin per DTTF vs DTTF quality", 6, 0, 6, 7, 1, 8);
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtInVsDttfQual->SetXTitle("n GMT input per DTTF");

  _nGmtOutVsDttfQual = _fs->make<TH2F>( "nGmtOutVsDttfQual", "number of GMTout per DTTF vs DTTF quality", 6, 0, 6, 7, 1, 8);
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtOutVsDttfQual->SetXTitle("n GMT output per DTTF");
  
  _dtQualityNewVsCategory = _fs->make<TH2F>( "dtQualityNewVsCategory", "dt quality", 7, 0, 7, 7, 0, 7);
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 4, "#splitline{DTTF + 1 GMTin}{+ 2 GMTout}");
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 5, "#splitline{DTTF + 0 GMTin}{+ 1 GMTout}");
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 6, "#splitline{DTTF + 0 GMTin}{+ 2 GMTout}");  
  _dtQualityNewVsCategory->GetXaxis()->SetBinLabel( 7, "#splitline{DTTF + > 1 GMTin}{or > 2 GMTout}");
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 1, "HI" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 2, "HO" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 3, "HI+RPC" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 4, "HO+RPC" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 5, "#splitline{(HI+HO)}{+RPC@bx0}" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 6, "(LL || HL)" );
  _dtQualityNewVsCategory->GetYaxis()->SetBinLabel( 7, "HH" );
  
  _nDttfVsGmtIn = _fs->make<TH2F>( "nDttfVsGmtIn", "number of DTTF vs number of GMTin", 10, 0, 10, 6, 0, 6);
  _nDttfVsGmtIn->SetXTitle("n DTTF");
  _nDttfVsGmtIn->SetYTitle("n GMTin");

  _nDttfVsGmtOut = _fs->make<TH2F>( "nDttfVsGmtOut", "number of DTTF vs number of GMTout", 10, 0, 10, 6, 0, 6);
  _nDttfVsGmtOut->SetXTitle("n DTTF");
  _nDttfVsGmtOut->SetYTitle("n GMTout");

  _nGmtInVsGmtOut = _fs->make<TH2F>( "nGmtInVsGmtOut", "number of GMTin vs number of GMTout", 10, 0, 10, 6, 0, 6);
  _nGmtInVsGmtOut->SetXTitle("n GMTin");
  _nGmtInVsGmtOut->SetYTitle("n GMTout");

  
  _positionGmtOutNoGmtIn = _fs->make<TH2F>( "positionGmtOutNoGmtIn", "position(#eta-#phi) of DTTF-matched GMTout with no GMTin", 30, -1.5, 1.5, 64, 0, 6.4);
  _positionGmtOutNoGmtIn->SetXTitle("GMTout #eta");
  _positionGmtOutNoGmtIn->SetYTitle("GMTout #phi");
  
  _GmtQualityVsEtaNoGmtIn = _fs->make<TH2F>( "GmtQualityVsEtaNoGmtIn", "GMTout quality", 30, -1.5, 1.5 ,7, 1, 8);
  _GmtQualityVsEtaNoGmtIn->SetXTitle("GMTout #eta");
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtQualityVsEtaNoGmtIn->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );
  
  _GmtOutQualCat3 = _fs->make<TH2F>( "GmtOutQualCat3", "GMTout quality for cat. 3", 7, 1, 8 ,7, 1, 8);
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtOutQualCat3->GetXaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtOutQualCat3->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );

  _GmtOutQualCat4 = _fs->make<TH1F>( "GmtOutQualCat4", "GMTout quality for cat. 4", 7, 1, 8);
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtOutQualCat4->GetXaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );
  
  _GmtOutQualCat5 = _fs->make<TH2F>( "GmtOutQualCat5", "GMTout quality for cat. 5", 7, 1, 8 ,7, 1, 8);
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtOutQualCat5->GetXaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _GmtOutQualCat5->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutQualCat3 = _fs->make<TH2F>( "DttfQualVsGmtOutQualCat3", "GMTout quality for cat. 3", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutQualCat3->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutQualCat3->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutHiQualCat3 = _fs->make<TH2F>( "DttfQualVsGmtOutHiQualCat3", "GMTout quality for cat. 3", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutHiQualCat3->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutHiQualCat3->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutLoQualCat3 = _fs->make<TH2F>( "DttfQualVsGmtOutLoQualCat3", "GMTout quality for cat. 3", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutLoQualCat3->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutLoQualCat3->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  
  
  _DttfQualVsGmtOutQualCat4 = _fs->make<TH2F>( "DttfQualVsGmtOutHiQualCat4", "GMTout quality for cat. 4", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutQualCat4->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutQualCat4->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutQualCat5 = _fs->make<TH2F>( "DttfQualVsGmtOutQualCat5", "GMTout quality for cat. 5", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutQualCat5->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutQualCat5->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutHiQualCat5 = _fs->make<TH2F>( "DttfQualVsGmtOutHiQualCat5", "GMTout quality for cat. 5", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutHiQualCat5->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutHiQualCat5->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  

  _DttfQualVsGmtOutLoQualCat5 = _fs->make<TH2F>( "DttfQualVsGmtOutLoQualCat5", "GMTout quality for cat. 5", 7, 1, 8 ,7, 1, 8);
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 1, "T34" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 2, "T23+T24" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 4, "T234" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 5, "T134" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 6, "T123+T124" );
  _DttfQualVsGmtOutLoQualCat5->GetXaxis()->SetBinLabel( 7, "T1234" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 1, "halo #mu" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 2, "VLQT1" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 3, "VLQT2" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 4, "VLQT3" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 5, "RCP unconf." );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 6, "#splitline{DT or CSC}{unconf.}" );
  _DttfQualVsGmtOutLoQualCat5->GetYaxis()->SetBinLabel( 7, "#splitline{DT/RPC or}{CSC/RPC unconf.}" );  
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

  int nTotGmtIn=0,  // nGmtIn per event
      nTotGmtOut=0, // nGmtOut per event
      nMBTr=0;      // nMBTr (Dttf) per event

  /// JP : loop over MBtracks
  edm::Handle<L1ITMu::MBTrackCollection> mbtrContainer;
  iEvent.getByLabel( _mbTracksCollectionInput, mbtrContainer);

  L1ITMu::MBTrackCollection::const_iterator tr = mbtrContainer->begin();
  L1ITMu::MBTrackCollection::const_iterator trend = mbtrContainer->end();

  nMBTr = mbtrContainer->size();
  
  for ( ; tr != trend; ++tr ) {
    
    const L1ITMu::MBTrack & mbtrack = *tr;
    
    int nGmtIn=0,
        nGmtOut=0;
    
    /// DTTF info
    const L1MuDTTrackCand& dttf = mbtrack.parent();
    
    int dttf_bx      =  dttf.bx();
    int dttf_qual    =  dttf.quality();
    int phi_local = dttf.phi_packed();
    if ( phi_local > 15 ) phi_local -= 32;
    double dttf_phi_global = (phi_local*(M_PI/72.))+((M_PI/6.)*dttf.scNum());
    if(dttf_phi_global < 0) dttf_phi_global+=2*M_PI;
    if(dttf_phi_global > 2*M_PI) dttf_phi_global-=2*M_PI;
    
    /// get GMTin (L1MuRegionalCand)
    const std::vector<L1MuRegionalCand> l1muregcand = mbtrack.getAssociatedGMTin(); 
    std::vector<L1MuRegionalCand>::const_iterator igmtin  = l1muregcand.begin();
    std::vector<L1MuRegionalCand>::const_iterator igmtinend = l1muregcand.end();
    
    for ( ; igmtin != igmtinend; ++igmtin ) {
      
      const L1MuRegionalCand & gmtin = *igmtin;
      
      if ( dttf_bx != gmtin.bx() ) continue;
      double dph = reco::deltaPhi( dttf_phi_global, gmtin.phiValue() );      
      _deltaPhiDttfGmtIn->Fill(dph);      
      
      vGmtIn.push_back(gmtin);        

    } // end loop over GMTin
    
    /// get GMTout (L1MuGMTExtendedCand)
    const std::vector<L1MuGMTExtendedCand> l1gmtextcand = mbtrack.getAssociatedGMTout(); 
    std::vector<L1MuGMTExtendedCand>::const_iterator igmtout  = l1gmtextcand.begin();
    std::vector<L1MuGMTExtendedCand>::const_iterator igmtoutend = l1gmtextcand.end();
    
    for ( ; igmtout != igmtoutend; ++igmtout ) {
      
      const L1MuGMTExtendedCand & gmtout = *igmtout;

      if ( dttf_bx != gmtout.bx() ) continue;
      double dph = reco::deltaPhi( dttf_phi_global, gmtout.phiValue() );      
      _deltaPhiDttfGmtOut->Fill(dph);
      
      vGmtOut.push_back(gmtout);        
      
    } // end loop over GMTout
    
    nGmtIn  = vGmtIn.size();
    nGmtOut = vGmtOut.size();
/// --- CATEGORY --- ///    
    // 0 -> DTTF only (GMTin=0 / GMTout=0)
    // 1 -> DTTF + 1 GMTin + 0 GMTout
    // 2 -> DTTF + 1 GMTin + 1 GMTout
    // 3 -> DTTF + 1 GMTin + 2 GMTout
    // 4 -> DTTF + 0 GMTin + 1 GMTout
    // 5 -> DTTF + 0 GMTin + 2 GMTout
    // 6 -> DTTF + > 1 GMTin || > 2 GMTout
    int Category = -9;          
    if ((nGmtIn == 0)&&(nGmtOut == 0)) Category = 0;
    if ((nGmtIn == 1)&&(nGmtOut == 0)) Category = 1;
    if ((nGmtIn == 1)&&(nGmtOut == 1)) Category = 2;
    if ((nGmtIn == 1)&&(nGmtOut == 2)) Category = 3;
    if ((nGmtIn == 0)&&(nGmtOut == 1)) Category = 4;
    if ((nGmtIn == 0)&&(nGmtOut == 2)) Category = 5;
    if ((nGmtIn  > 1)||(nGmtOut  > 2)) Category = 6;
//     if ((nGmtIn < 2) && (nGmtOut < 3))
//     {
//       if ( !nGmtIn && nGmtOut ) Category = 3 + nGmtOut;
//       else Category = nGmtIn + nGmtOut;
//     }
//     else Category = 6;
/// ---          --- ///
    
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
        if ( dttf_bx != dt->getDTData().bx ) continue;            
          
        _dtQualityNewVsCategory->Fill(Category,dtquality);
      }
      
    } // end loop over MBStations 
    
    if (nGmtIn == 0){
      
      std::vector<L1MuGMTExtendedCand>::const_iterator it  = vGmtOut.begin();
      std::vector<L1MuGMTExtendedCand>::const_iterator itend = vGmtOut.end();
      for( ; it != itend; ++it){        
        const L1MuGMTExtendedCand & gmtout = *it;
        _positionGmtOutNoGmtIn->Fill(gmtout.etaValue(),gmtout.phiValue());
        _GmtQualityVsEtaNoGmtIn->Fill(gmtout.etaValue(),gmtout.quality());
      }          
    }
    
/// FILL HISTOS    
    /// SAME BX as DTTF   

    _nGmtIn->Fill(nGmtIn);
    _nGmtOut->Fill(nGmtOut);
    _nGmtInVsGmtOutPerDTTF->Fill(nGmtIn,nGmtOut);
    _nGmtInVsDttfQual->Fill(nGmtIn,dttf_qual);
    _nGmtOutVsDttfQual->Fill(nGmtOut,dttf_qual);
    _dttfQualityVsCategory->Fill(Category,dttf_qual);
   
    nTotGmtIn       += nGmtIn;
    nTotGmtOut       += nGmtOut;
    
    if (Category==3) // 1GMTin 2GMTout
    {
      L1MuGMTExtendedCand GmtOut_1 = vGmtOut.at(0);
      L1MuGMTExtendedCand GmtOut_2 = vGmtOut.at(1);
      
      int Qual_1 = GmtOut_1.quality();
      int Qual_2 = GmtOut_2.quality();
      
      _DttfQualVsGmtOutQualCat3->Fill(dttf_qual,Qual_1); _DttfQualVsGmtOutQualCat3->Fill(dttf_qual,Qual_2);

      if (Qual_1 < Qual_2) {
        _GmtOutQualCat3->Fill(Qual_1,Qual_2);        
        _DttfQualVsGmtOutHiQualCat3->Fill(dttf_qual,Qual_2);        
        _DttfQualVsGmtOutLoQualCat3->Fill(dttf_qual,Qual_1);        
      }
      else {
        _GmtOutQualCat3->Fill(Qual_2,Qual_1);
        _DttfQualVsGmtOutHiQualCat3->Fill(dttf_qual,Qual_1);                
        _DttfQualVsGmtOutLoQualCat3->Fill(dttf_qual,Qual_2);                
      }
    }

    if (Category==4) // 0GMTin 1GMTout
    {
      L1MuGMTExtendedCand GmtOut = vGmtOut.at(0);
      
      int Qual = GmtOut.quality();
      _GmtOutQualCat4->Fill(Qual);
      _DttfQualVsGmtOutQualCat4->Fill(dttf_qual,Qual);        
    }    
    
    if (Category==5) // 1GMTin 2GMTout
    {
      L1MuGMTExtendedCand GmtOut_1 = vGmtOut.at(0);
      L1MuGMTExtendedCand GmtOut_2 = vGmtOut.at(1);
      
      int Qual_1 = GmtOut_1.quality();
      int Qual_2 = GmtOut_2.quality();

      _DttfQualVsGmtOutQualCat5->Fill(dttf_qual,Qual_1); _DttfQualVsGmtOutQualCat5->Fill(dttf_qual,Qual_2);

      if (Qual_1 < Qual_2) {
        _GmtOutQualCat5->Fill(Qual_1,Qual_2);        
        _DttfQualVsGmtOutHiQualCat5->Fill(dttf_qual,Qual_2);        
        _DttfQualVsGmtOutLoQualCat5->Fill(dttf_qual,Qual_1);        
      }
      else {
        _GmtOutQualCat5->Fill(Qual_2,Qual_1);
        _DttfQualVsGmtOutHiQualCat5->Fill(dttf_qual,Qual_1);                
        _DttfQualVsGmtOutLoQualCat5->Fill(dttf_qual,Qual_2);                        
      }
    }
    
  }
  
  _nDttfVsGmtIn     ->Fill(nMBTr,nTotGmtIn);

  _nDttfVsGmtOut     ->Fill(nMBTr,nTotGmtOut);
  
  _nGmtInVsGmtOut     ->Fill(nTotGmtIn,nTotGmtOut);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelTrackPlots);
