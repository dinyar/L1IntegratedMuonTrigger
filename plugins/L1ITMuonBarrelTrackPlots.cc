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

#include "DataFormats/MuonDetId/interface/DTChamberId.h"

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
  double GetPhiGlobal( const L1MuDTTrackCand& );

  edm::InputTag _mbltCollectionInput;
  edm::InputTag _mbTracksCollectionInput;
  edm::InputTag _l1itDtPhiChInput;
  edm::Service<TFileService> _fs;
  
  TH1F * _deltaPhiDttfGmtIn; // DeltaPhi (DTTF - GMTinput) -> effect of the "phi_packed -> phi_global" conversion
  TH1F * _deltaPhiDttfGmtOut; // DeltaPhi (DTTF - GMToutput) -> should always be 0

  TH1F * _nGmtInPerTrack; // Number of GMTinput per MBTrack
  TH1F * _nGmtOutPerTrack; // Number of GMToutput per MBTrack
  TH2F * _nGmtInVsGmtOutPerTrack; // Number of GMTinput vs Number of GMToutput per MBTrack
  
  TH1F * _nTrack; // Number of MBTrack/DTTF per event
  TH1F * _nGmtIn; // Number of GMTinput per event
  TH1F * _nGmtOut; // Number of GMToutput per event
  TH2F * _nTrackVsGmtIn; // Number of MBTrack/DTTF vs Number of GMTinput per event
  TH2F * _nTrackVsGmtOut; // Number of MBTrack/DTTF vs Number of GMToutput per event
  TH2F * _nGmtInVsGmtOut; // Number of GMTinput vs Number of GMToutput per event
  
  TH2F * _nGmtInVsDttfQual; // Number of GMTinput per MBTrack vs DTTF quality
  TH2F * _nGmtOutVsDttfQual; // Number of GMToutput per MBTrack vs DTTF quality 

  TH2F * _dttfQualityVsCategory; // DTTF quality vs MBTrack category (only DTTF / DTTF + GMTin / DTTF + GMTout)
  TH2F * _dtQualityVsCategory;  // DTPrimitive (for each DTPrim in a DTTF) quality vs MBTrack category (only DTTF / DTTF + GMTin / DTTF + GMTout)

  /// NEW
  
  TH1F * _deltaBxDtRpcIn; // DeltaBx : DTPrimitive - closest (in DeltaPhi) RPCIn hit (if any)
  TH1F * _deltaBxDtRpcOut; // DeltaBx : DTPrimitive - closest (in DeltaPhi) RPCOut hit (if any)
  TH1F * _deltaBxDtRpc; // DeltaBx : DTPrimitive - closest (in DeltaPhi) RPCIn hit (if any) ; DTPrimitive - closest (in DeltaPhi) RPCOut hit (if any)
  
  TH2F * _deltaBxDtRpcInOutVsDttfQuality; // deltaBxDtRpcIn vs DTTF quality
  TH2F * _deltaBxDtRpcInVsDttfQuality; // deltaBxDtRpcOut vs DTTF quality
  TH2F * _deltaBxDtRpcOutVsDttfQuality; // deltaBxDtRpc vs DTTF quality

  TH2F * _deltaBxDtRpcInOutVsDtQuality; // deltaBxDtRpcIn vs DT Primitive quality
  TH2F * _deltaBxDtRpcInVsDtQuality; // deltaBxDtRpcOut vs DT Primitive quality
  TH2F * _deltaBxDtRpcOutVsDtQuality; // deltaBxDtRpc vs DT Primitive quality

  TH1F * _nRpcIn; // nRpcIn per DT Trigger Primitive (for the DT primitives accessed by the DTTF address)
  TH1F * _nRpcOut; // nRpcOut per DT Trigger Primitive (for the DT primitives accessed by the DTTF address)
  TH1F * _nRpc; // nRpcIn + nRpcOut per DT Trigger Primitive (for the DT primitives accessed by the DTTF address)
  
  TH1F * _DTQuality[7]; // DT TriggerPrimitive quality per each DTTF quality
  TH1F * _deltaBxRpcIn[7]; // deltaBxDtRpcIn per each DTTF quality
  TH1F * _deltaBxRpcOut[7]; // deltaBxDtRpcOut per each DTTF quality
  TH1F * _deltaBxRpc[7]; // deltaBxDtRpc per each DTTF quality

  TH2F * _dtQualityVsConfirm; // DT TriggerPrimitive quality vs Confirmation ( DT unconfirmed / DT + RPC out / DT + RPC in / DT + RPC in+out )
  
  TH2F * _deltaBxDtRpcInVsOut; // deltaBxDtRpcIn vs deltaBxDtRpcOut per each confirmed DT
};


L1ITMuonBarrelTrackPlots::L1ITMuonBarrelTrackPlots(const PSet& p)
  : _mbltCollectionInput( p.getParameter<edm::InputTag>("MBLTCollection") ),
    _mbTracksCollectionInput( p.getParameter<edm::InputTag>("MBTracksCollection") ),
    _l1itDtPhiChInput( p.getParameter<edm::InputTag>("L1ITDTChambPhContainer") )
{

  _deltaPhiDttfGmtIn       = _fs->make<TH1F>( "deltaPhiDttfGmtIn", "deltaPhiDttfGmtIn", 200, -0.00001, 0.00001 );
  _deltaPhiDttfGmtIn->SetXTitle("#Delta#phi");

  _deltaPhiDttfGmtOut      = _fs->make<TH1F>( "deltaPhiDttfGmtOut", "deltaPhiDttfGmtOut", 600, -0.3, 0.3 );
  _deltaPhiDttfGmtOut->SetXTitle("#Delta#phi");
  
  _dttfQualityVsCategory = _fs->make<TH2F>( "dttfQual", "DTTF quality", 4, 0, 3, 7, 1, 8);
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dttfQualityVsCategory->GetXaxis()->SetBinLabel( 4, "#splitline{Other}{cases}");
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 1, "T34" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 4, "T234" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 5, "T134" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _dttfQualityVsCategory->GetYaxis()->SetBinLabel( 7, "T1234" );
 
  _nGmtInPerTrack = _fs->make<TH1F>( "nGmtInPerTrack",  "number of GMTin per MBTrack", 4, 0, 4 );
  _nGmtInPerTrack->SetXTitle("n GMT input per MBTrack");

  _nGmtOutPerTrack      = _fs->make<TH1F>( "nGmtOutPerTrack", "number of GMTout per MBTrack", 4, 0, 4 );
  _nGmtOutPerTrack->SetXTitle("n GMT output per MBTrack");

  _nGmtInVsGmtOutPerTrack      = _fs->make<TH2F>( "nGmtInVsGmtOutPerTrack", "number of GMTin Vs number of GMTout per MBTrack", 4, 0, 4 , 4, 0, 4 );
  _nGmtInVsGmtOutPerTrack->SetXTitle("n GMT input per MBTrack");
  _nGmtInVsGmtOutPerTrack->SetYTitle("n GMT output per MBTrack");

  _nTrack       = _fs->make<TH1F>( "nTrack",  "number of MBTrack per Event", 10, 0, 10 );
  _nTrack->SetXTitle("n MBTrack per Event");
  
  _nGmtIn       = _fs->make<TH1F>( "nGmtIn",  "number of GMTin per Event", 4, 0, 4 );
  _nGmtIn->SetXTitle("n GMT input per Event");

  _nGmtOut      = _fs->make<TH1F>( "nGmtOut", "number of GMTout per Event", 4, 0, 4 );  
  _nGmtOut->SetXTitle("n GMT output per Event");
  
  _nTrackVsGmtIn = _fs->make<TH2F>( "nTrackVsGmtIn", "number of MBTrack vs number of GMTin per Event", 10, 0, 10, 4, 0, 4);
  _nTrackVsGmtIn->SetXTitle("n MBTrack per Event");
  _nTrackVsGmtIn->SetYTitle("n GMT input per Event");

  _nTrackVsGmtOut = _fs->make<TH2F>( "nTrackVsGmtOut", "number of MBTrack vs number of GMTout per Event", 10, 0, 10, 4, 0, 4);
  _nTrackVsGmtOut->SetXTitle("n MBTrack per Event");
  _nTrackVsGmtOut->SetYTitle("n GMT output per Event");

  _nGmtInVsGmtOut = _fs->make<TH2F>( "nGmtInVsGmtOut", "number of GMTin vs number of GMTout per Event", 4, 0, 4, 4, 0, 4);
  _nGmtInVsGmtOut->SetXTitle("n GMT input");
  _nGmtInVsGmtOut->SetYTitle("n GMT output");

  
  _nGmtInVsDttfQual = _fs->make<TH2F>( "nGmtInVsDttfQual", "number of GMTin per MBTrack vs DTTF quality", 4, 0, 4, 7, 1, 8);
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtInVsDttfQual->SetXTitle("n GMT input per MBTrack");

  _nGmtOutVsDttfQual = _fs->make<TH2F>( "nGmtOutVsDttfQual", "number of GMTout per MBTrack vs DTTF quality", 4, 0, 4, 7, 1, 8);
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtOutVsDttfQual->SetXTitle("n GMT output per MBTrack");
  
  _dtQualityVsCategory = _fs->make<TH2F>( "dtQualityVsCategory", "dt quality", 4, 0, 4, 7, 0, 7);
  _dtQualityVsCategory->GetXaxis()->SetBinLabel( 1, "DTTF only");
  _dtQualityVsCategory->GetXaxis()->SetBinLabel( 2, "#splitline{DTTF + 1 GMTin}{+ 0 GMTout}");
  _dtQualityVsCategory->GetXaxis()->SetBinLabel( 3, "#splitline{DTTF + 1 GMTin}{+ 1 GMTout}");
  _dtQualityVsCategory->GetXaxis()->SetBinLabel( 4, "#splitline{Other}{cases}");
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 1, "LI" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 2, "LO" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 3, "HI" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 4, "HO" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 5, "LL" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 6, "HL" );
  _dtQualityVsCategory->GetYaxis()->SetBinLabel( 7, "HH" );
  
    _deltaBxDtRpcIn  = _fs->make<TH1F>( "deltaBxDtRpcIn", "deltaBxDtRpcIn", 10, -3, 7);
  _deltaBxDtRpcIn->SetXTitle("#Delta Bx");

  _deltaBxDtRpcOut = _fs->make<TH1F>( "deltaBxDtRpcOut", "deltaBxDtRpcOut", 10, -3, 7);
  _deltaBxDtRpcOut->SetXTitle("#Delta Bx");

  _deltaBxDtRpc    = _fs->make<TH1F>( "deltaBxDtRpc", "deltaBxDtRpc", 10, -3, 7);
  _deltaBxDtRpc->SetXTitle("#Delta Bx");

  _deltaBxDtRpcInOutVsDttfQuality    = _fs->make<TH2F>( "deltaBxDtRpcInOutVsDttfQuality", "deltaBxDtRpc vs DTTF quality", 10, -3, 7, 7, 1, 8);
  _deltaBxDtRpcInOutVsDttfQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 1, "T34" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 4, "T234" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 5, "T134" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _deltaBxDtRpcInOutVsDttfQuality->GetYaxis()->SetBinLabel( 7, "T1234" );

  _deltaBxDtRpcInOutVsDtQuality    = _fs->make<TH2F>( "deltaBxDtRpcInOutVsDtQuality", "deltaBxDtRpc vs DT quality", 10, -3, 7, 7, 0, 7);
  _deltaBxDtRpcInOutVsDtQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 1, "LI" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 2, "LO" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 3, "HI" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 4, "HO" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 5, "LL" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 6, "HL" );
  _deltaBxDtRpcInOutVsDtQuality->GetYaxis()->SetBinLabel( 7, "HH" );

  _deltaBxDtRpcInVsDttfQuality    = _fs->make<TH2F>( "deltaBxDtRpcInVsDttfQuality", "deltaBxDtRpcIn vs DTTF quality", 10, -3, 7, 7, 1, 8);
  _deltaBxDtRpcInVsDttfQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 1, "T34" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 4, "T234" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 5, "T134" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _deltaBxDtRpcInVsDttfQuality->GetYaxis()->SetBinLabel( 7, "T1234" );

  _deltaBxDtRpcInVsDtQuality    = _fs->make<TH2F>( "deltaBxDtRpcInVsDtQuality", "deltaBxDtRpcIn vs DT quality", 10, -3, 7, 7, 0, 7);
  _deltaBxDtRpcInVsDtQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 1, "LI" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 2, "LO" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 3, "HI" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 4, "HO" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 5, "LL" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 6, "HL" );
  _deltaBxDtRpcInVsDtQuality->GetYaxis()->SetBinLabel( 7, "HH" );
  
  _deltaBxDtRpcOutVsDttfQuality    = _fs->make<TH2F>( "deltaBxDtRpcOutVsDttfQuality", "deltaBxDtRpcOut vs DTTF quality", 10, -3, 7, 7, 1, 8);
  _deltaBxDtRpcOutVsDttfQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 1, "T34" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 4, "T234" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 5, "T134" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _deltaBxDtRpcOutVsDttfQuality->GetYaxis()->SetBinLabel( 7, "T1234" );

  _deltaBxDtRpcOutVsDtQuality    = _fs->make<TH2F>( "deltaBxDtRpcOutVsDtQuality", "deltaBxDtRpcOut vs DT quality", 10, -3, 7, 7, 0, 7);
  _deltaBxDtRpcOutVsDtQuality->SetXTitle("#Delta Bx");  
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 1, "LI" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 2, "LO" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 3, "HI" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 4, "HO" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 5, "LL" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 6, "HL" );
  _deltaBxDtRpcOutVsDtQuality->GetYaxis()->SetBinLabel( 7, "HH" );
  
  _nRpcIn  = _fs->make<TH1F>( "nRpcIn", "nRpcIn", 10, 0, 10);
  _nRpcIn->SetXTitle("n matched RPCin");

  _nRpcOut  = _fs->make<TH1F>( "nRpcOut", "nRpcOut", 10, 0, 10);
  _nRpcOut->SetXTitle("n matched RPCout");

  _nRpc  = _fs->make<TH1F>( "nRpc", "nRpc", 10, 0, 10);
  _nRpc->SetXTitle("n matched RPC");
  
  for (unsigned i = 0 ; i < 7 ; ++i){

    _deltaBxRpcIn[i]  = _fs->make<TH1F>( Form("deltaBxRpcIn_dttf%d",i+1), Form("deltaBxDtRpcIn - DTTF quality %d",i+1), 10, -3, 7);
    _deltaBxRpcIn[i]->SetXTitle("#Delta Bx");
    
    _deltaBxRpcOut[i]  = _fs->make<TH1F>( Form("deltaBxRpcOut_dttf%d",i+1), Form("deltaBxDtRpcOut - DTTF quality %d",i+1), 10, -3, 7);
    _deltaBxRpcOut[i]->SetXTitle("#Delta Bx");
    
    _deltaBxRpc[i]  = _fs->make<TH1F>( Form("deltaBxRpc_dttf%d",i+1), Form("deltaBxDtRpc - DTTF quality %d",i+1), 10, -3, 7);
    _deltaBxRpc[i]->SetXTitle("#Delta Bx");
        
    _DTQuality[i] = _fs->make<TH1F>( Form("DTQuality_dttf%d",i+1), Form("DT quality - DTTF quality %d",i+1), 7, 0, 7);
    _DTQuality[i]->GetXaxis()->SetBinLabel( 1, "LI" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 2, "LO" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 3, "HI" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 4, "HO" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 5, "LL" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 6, "HL" );
    _DTQuality[i]->GetXaxis()->SetBinLabel( 7, "HH" );
    
  }

  _dtQualityVsConfirm = _fs->make<TH2F>( "dtQualityVsConfirm", "dt quality", 4, 0, 4, 7, 0, 7);
  _dtQualityVsConfirm->GetXaxis()->SetBinLabel( 1, "DT only");
  _dtQualityVsConfirm->GetXaxis()->SetBinLabel( 2, "DT + RPCin");
  _dtQualityVsConfirm->GetXaxis()->SetBinLabel( 3, "DT + RPCout");
  _dtQualityVsConfirm->GetXaxis()->SetBinLabel( 4, "DT + 2 RPC");
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 1, "LI" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 2, "LO" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 3, "HI" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 4, "HO" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 5, "LL" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 6, "HL" );
  _dtQualityVsConfirm->GetYaxis()->SetBinLabel( 7, "HH" );

  _deltaBxDtRpcInVsOut = _fs->make<TH2F>( "deltaBxDtRpcInVsOut", "deltaBxDtRpcInVsOut", 10, -3, 7, 10, -3, 7);
  _deltaBxDtRpcInVsOut->SetYTitle("#Delta Bx");
  _deltaBxDtRpcInVsOut->SetYTitle("#Delta Bx");
}
  
L1ITMuonBarrelTrackPlots::~L1ITMuonBarrelTrackPlots()
{
  
}

double L1ITMuonBarrelTrackPlots::GetPhiGlobal( const L1MuDTTrackCand& dttf_){
     
  int lphi = dttf_.phi_packed();
  double gphi = -999;
  
  if ( lphi > 15 ) lphi -= 32;
  gphi = (lphi*(M_PI/72.))+((M_PI/6.)*dttf_.scNum());
  if(gphi < 0) gphi+=2*M_PI;
  if(gphi > 2*M_PI) gphi-=2*M_PI;
  
  return gphi;
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
  
  int nTotGmtIn =0,  // nGmtIn per event
      nTotGmtOut=0,  // nGmtOut per event
      nTotMBTr  =0;  // nMBTr per event
  
  /// JP : loop over MBtracks
  edm::Handle<L1ITMu::MBTrackCollection> mbtrContainer;
  iEvent.getByLabel( _mbTracksCollectionInput, mbtrContainer);

  L1ITMu::MBTrackCollection::const_iterator tr = mbtrContainer->begin();
  L1ITMu::MBTrackCollection::const_iterator trend = mbtrContainer->end();

  nTotMBTr = mbtrContainer->size();
  
  for ( ; tr != trend; ++tr ) {
    
    const L1ITMu::MBTrack & mbtrack = *tr;
    
    int nGmtIn=0,
        nGmtOut=0;
    std::vector<L1MuRegionalCand> vGmtIn;
    std::vector<L1MuGMTExtendedCand> vGmtOut;
  
    /// DTTF info
    const L1MuDTTrackCand& dttf = mbtrack.parent();
    int dttf_bx      =  dttf.bx();
    int dttf_qual    =  dttf.quality();
    double dttf_phi_global = GetPhiGlobal(dttf);
    
    
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
    // 3 -> Any other case
    int Category = 3;    
    if ((nGmtIn == 0)&&(nGmtOut == 0)){
      Category = 0;      
    } else if ((nGmtIn == 1)&&(nGmtOut == 0)) {
        Category = 1;
    } else if ((nGmtIn == 1)&&(nGmtOut == 1)) {
        Category = 2;
    }
/// ---          --- ///
    
    /// loop over MBCollection
    const L1ITMu::MBLTVectorRef & mbltPairs = mbtrack.getStubs();
    L1ITMu::MBLTVectorRef::const_iterator mbtp = mbltPairs.begin();
    L1ITMu::MBLTVectorRef::const_iterator mbtpend = mbltPairs.end();
        
    for ( ; mbtp != mbtpend; ++mbtp ) {

      const L1ITMu::MBLTContainerRef    & mbcont = (*mbtp).first;  // the MBContainer  (map<DTChamberId,MBCollection>)
      const L1ITMu::TriggerPrimitiveRef & dttp   = (*mbtp).second; // the TriggerPrimitive associated to the MBCollection
  
// // //       const DTChamberId            & dtid        = mbcont->first;  // the DTChamberId
      const L1ITMu::MBLTCollection & mbltStation = mbcont->second; // the MBCollection
    
      /// consistency check
      int station = mbltStation.station();
      if ( station < 1 || station > 4 )
        throw cms::Exception("Invalid Station")
          << "wrong station number " << station << std::endl;

      /// the actual DT primitive      
      _DTQuality[dttf_qual-1]->Fill(dttp->getDTData().qualityCode);
//       std::cout <<
//             "\n\tThe DTstub matched with the DTTF = " << 
//             "\n\t\tDTstub->wheel          = " << dttp->getDTData().wheel <<
//             "\n\t\tDTstub->station        = " << dttp->getDTData().station <<
//             "\n\t\tDTstub->sector         = " << dttp->getDTData().sector <<
//             "\n\t\tDTstub->qualityCode    = " << dttp->getDTData().qualityCode <<
//             "\n\t\tDTstub->bx             = " << dttp->getDTData().bx <<
//             "\n\t\tDTstub->Ts2TagCode     = " << dttp->getDTData().Ts2TagCode;            

      /// access all the DT primitives in the MBCollection
      size_t dtListSize = mbltStation.getDtSegments().size();

      for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {
        
        const L1ITMu::TriggerPrimitiveRef & dt = mbltStation.getDtSegments().at(iDt);

        int dtquality   = dt->getDTData().qualityCode;
        int dtcorrectbx = dt->getDTData().bx - dt->getDTData().Ts2TagCode;

        if (dtquality < 0) continue; // skipping theta chambers
//         std::cout <<
//             "\n\t\tDTstubs["<< iDt <<"]->wheel          = " << dt->getDTData().wheel <<
//             "\n\t\tDTstubs["<< iDt <<"]->station        = " << dt->getDTData().station <<
//             "\n\t\tDTstubs["<< iDt <<"]->sector         = " << dt->getDTData().sector <<
//             "\n\t\tDTstubs["<< iDt <<"]->qualityCode    = " << dt->getDTData().qualityCode <<
//             "\n\t\tDTstubs["<< iDt <<"]->bx             = " << dt->getDTData().bx <<
//             "\n\t\tDTstubs["<< iDt <<"]->Ts2TagCode     = " << dt->getDTData().Ts2TagCode;

        /// access the RPC associated to the DT primitive
        if (dt == dttp) { 

          const L1ITMu::TriggerPrimitiveList rpcins   =  mbltStation.getRpcInAssociatedStubs(iDt); // the List of RPCin TriggerPrimitive associated to the DT
          const L1ITMu::TriggerPrimitiveList rpcouts  =  mbltStation.getRpcOutAssociatedStubs(iDt); // the List of RPCout TriggerPrimitive associated to the DT

          _nRpcIn->Fill(rpcins.size()); 
          _nRpcOut->Fill(rpcouts.size());
          _nRpc->Fill(rpcins.size()+rpcouts.size());

          /// --- CONFIRM --- ///    
          // 0 -> DT only (RpcIn = 0 / RpcOut = 0)
          // 1 -> DT + RpcIn (RpcIn > 0 / RpcOut = 0)
          // 2 -> DT + RpcOut (RpcIn = 0 / RpcOut > 0)
          // 3 -> DT + RpcIn + RpcOut (RpcIn > 0 / RpcOut > 0)                  
          int Confirm = 0;
          if (rpcins.size() && rpcouts.size()){
            Confirm = 3;
            
            int rpcbxIn = rpcins.front()->getBX();
            int rpcbxOut = rpcouts.front()->getBX();
          
            _deltaBxDtRpcInVsOut->Fill(rpcbxIn - dtcorrectbx, rpcbxOut - dtcorrectbx);
            _deltaBxDtRpc->Fill(rpcbxIn - dtcorrectbx);   _deltaBxRpc[dttf_qual-1]->Fill(rpcbxIn - dtcorrectbx); 
            _deltaBxDtRpc->Fill(rpcbxOut - dtcorrectbx);    _deltaBxRpc[dttf_qual-1]->Fill(rpcbxOut - dtcorrectbx);                  
            _deltaBxDtRpcInOutVsDttfQuality->Fill(rpcbxIn - dtcorrectbx, dttf_qual);
            _deltaBxDtRpcInOutVsDttfQuality->Fill(rpcbxOut - dtcorrectbx, dttf_qual);
            _deltaBxDtRpcInOutVsDtQuality->Fill(rpcbxIn - dtcorrectbx, dtquality);
            _deltaBxDtRpcInOutVsDtQuality->Fill(rpcbxOut - dtcorrectbx, dtquality);                                 
            
          } else if ( rpcins.size() ) {
              Confirm = 1;
              
              int rpcbxIn = rpcins.front()->getBX();
              
              _deltaBxDtRpcIn->Fill(rpcbxIn - dtcorrectbx); _deltaBxRpcIn[dttf_qual-1]->Fill(rpcbxIn - dtcorrectbx);
              _deltaBxDtRpcInVsDttfQuality->Fill(rpcbxIn - dtcorrectbx, dttf_qual);
              _deltaBxDtRpcInVsDtQuality->Fill(rpcbxIn - dtcorrectbx, dtquality);
              
          } else if ( rpcouts.size() ) {
              Confirm = 2;
              
              int rpcbxOut = rpcouts.front()->getBX();          
              
              _deltaBxDtRpcOut->Fill(rpcbxOut - dtcorrectbx); _deltaBxRpcOut[dttf_qual-1]->Fill(rpcbxOut - dtcorrectbx);                      
              _deltaBxDtRpcOutVsDttfQuality->Fill(rpcbxOut - dtcorrectbx, dttf_qual);
              _deltaBxDtRpcOutVsDtQuality->Fill(rpcbxOut - dtcorrectbx, dtquality);                                 
              
          }
          _dtQualityVsConfirm->Fill(Confirm,dtquality);
          _dtQualityVsCategory->Fill(Category,dtquality);
                                           
        } // end if dt == dttp
                                 
      } // end loop over DT TriggerPrimitives
                     
    } // end loop over MBStations 
       
    /// fill histos per MBTrack   
    _nGmtInPerTrack->Fill(nGmtIn);
    _nGmtOutPerTrack->Fill(nGmtOut);    
    _nGmtInVsGmtOutPerTrack->Fill(nGmtIn,nGmtOut);
    
    _nGmtInVsDttfQual->Fill(nGmtIn,dttf_qual);
    _nGmtOutVsDttfQual->Fill(nGmtOut,dttf_qual);
    _dttfQualityVsCategory->Fill(Category,dttf_qual);
       
    nTotGmtIn  += nGmtIn;
    nTotGmtOut += nGmtOut;
    
    vGmtIn.clear();
    vGmtOut.clear();
  }

  /// fill histos per Event   
  _nTrack ->Fill(nTotMBTr);    
  _nGmtIn ->Fill(nTotGmtIn);  
  _nGmtOut->Fill(nTotGmtOut);
  
  _nTrackVsGmtIn ->Fill(nTotMBTr,nTotGmtIn);
  _nTrackVsGmtOut->Fill(nTotMBTr,nTotGmtOut);  
  _nGmtInVsGmtOut->Fill(nTotGmtIn,nTotGmtOut);
   
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMuonBarrelTrackPlots);
