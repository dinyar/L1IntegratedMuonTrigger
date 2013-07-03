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

  
  
  TH2F * _nGmtInVsDttfQual;
  TH2F * _nGmtOutVsDttfQual; 

  TH2F * _dttfQualityVsCategory;  
  TH2F * _dtQualityNewVsCategory;
  
  
  
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
  _deltaPhiDttfGmtIn->SetXTitle("#Delta#phi");

  _deltaPhiDttfGmtOut      = _fs->make<TH1F>( "deltaPhiDttfGmtOut", "deltaPhiDttfGmtOut", 600, -0.3, 0.3 );
  _deltaPhiDttfGmtOut->SetXTitle("#Delta#phi");
  
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
 
  _nGmtInPerTrack = _fs->make<TH1F>( "nGmtInPerTrack",  "number of GMTin per MBTrack", 6, 0, 6 );
  _nGmtInPerTrack->SetXTitle("n GMT input per MBTrack");

  _nGmtOutPerTrack      = _fs->make<TH1F>( "nGmtOutPerTrack", "number of GMTout per MBTrack", 6, 0, 6 );
  _nGmtOutPerTrack->SetXTitle("n GMT output per MBTrack");

  _nGmtInVsGmtOutPerTrack      = _fs->make<TH2F>( "nGmtInVsGmtOutPerTrack", "number of GMTin Vs number of GMTout per MBTrack", 6, 0, 6 , 6, 0, 6 );
  _nGmtInVsGmtOutPerTrack->SetXTitle("n GMT input per MBTrack");
  _nGmtInVsGmtOutPerTrack->SetYTitle("n GMT output per MBTrack");

  _nTrack       = _fs->make<TH1F>( "nTrack",  "number of MBTrack per Event", 10, 0, 10 );
  _nTrack->SetXTitle("n MBTrack per Event");
  
  _nGmtIn       = _fs->make<TH1F>( "nGmtIn",  "number of GMTin per Event", 6, 0, 6 );
  _nGmtIn->SetXTitle("n GMT input per Event");

  _nGmtOut      = _fs->make<TH1F>( "nGmtOut", "number of GMTout per Event", 6, 0, 6 );  
  _nGmtOut->SetXTitle("n GMT output per Event");
  
  _nTrackVsGmtIn = _fs->make<TH2F>( "nTrackVsGmtIn", "number of MBTrack vs number of GMTin per Event", 10, 0, 10, 6, 0, 6);
  _nTrackVsGmtIn->SetXTitle("n MBTrack per Event");
  _nTrackVsGmtIn->SetYTitle("n GMT input per Event");

  _nTrackVsGmtOut = _fs->make<TH2F>( "nTrackVsGmtOut", "number of MBTrack vs number of GMTout per Event", 10, 0, 10, 6, 0, 6);
  _nTrackVsGmtOut->SetXTitle("n MBTrack per Event");
  _nTrackVsGmtOut->SetYTitle("n GMT output per Event");

  _nGmtInVsGmtOut = _fs->make<TH2F>( "nGmtInVsGmtOut", "number of GMTin vs number of GMTout per Event", 10, 0, 10, 6, 0, 6);
  _nGmtInVsGmtOut->SetXTitle("n GMT input");
  _nGmtInVsGmtOut->SetYTitle("n GMT output");

  
  _nGmtInVsDttfQual = _fs->make<TH2F>( "nGmtInVsDttfQual", "number of GMTin per MBTrack vs DTTF quality", 6, 0, 6, 7, 1, 8);
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtInVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtInVsDttfQual->SetXTitle("n GMT input per MBTrack");

  _nGmtOutVsDttfQual = _fs->make<TH2F>( "nGmtOutVsDttfQual", "number of GMTout per MBTrack vs DTTF quality", 6, 0, 6, 7, 1, 8);
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 1, "T34" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 2, "T23+T24" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 3, "T12+T13+T14" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 4, "T234" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 5, "T134" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 6, "T123+T124" );
  _nGmtOutVsDttfQual->GetYaxis()->SetBinLabel( 7, "T1234" );
  _nGmtOutVsDttfQual->SetXTitle("n GMT output per MBTrack");
  
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
/// ---          --- ///
    
    /// loop over MBCollection
    const L1ITMu::MBLTVectorRef & mbltPairs = mbtrack.getStubs();
    L1ITMu::MBLTVectorRef::const_iterator mbtp = mbltPairs.begin();
    L1ITMu::MBLTVectorRef::const_iterator mbtpend = mbltPairs.end();
        
    for ( ; mbtp != mbtpend; ++mbtp ) {

      const L1ITMu::MBLTContainerRef    & mbcont = (*mbtp).first;
      const L1ITMu::TriggerPrimitiveRef & dttp   = (*mbtp).second;
  
      const DTChamberId            & dtid        = mbcont->first;
      const L1ITMu::MBLTCollection & mbltStation = mbcont->second;
    
      /// useful index
      int station = mbltStation.station();
      
      /// consistency check
      if ( station < 1 || station > 4 )
        throw cms::Exception("Invalid Station")
          << "wrong station number " << station << std::endl;

      /// size for primitives vectors
      size_t dtListSize = mbltStation.getDtSegments().size();
//       size_t rpcInListSize = mbltStation.getRpcInner().size();
//       size_t rpcOutListSize = mbltStation.getRpcOuter().size();      

      /// the actual DT primitive
//       std::cout <<
//             "\n\tThe DTstub matched with the DTTF = " << 
//             "\n\t\tDTstub->wheel          = " << dttp->getDTData().wheel <<
//             "\n\t\tDTstub->station        = " << dttp->getDTData().station <<
//             "\n\t\tDTstub->sector         = " << dttp->getDTData().sector <<
//             "\n\t\tDTstub->qualityCode    = " << dttp->getDTData().qualityCode <<
//             "\n\t\tDTstub->bx             = " << dttp->getDTData().bx <<
//             "\n\t\tDTstub->Ts2TagCode     = " << dttp->getDTData().Ts2TagCode;            

      /// the list of all the DT primitives in the MBCollection
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
                 
        _dtQualityNewVsCategory->Fill(Category,dtquality);
      }
           
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
