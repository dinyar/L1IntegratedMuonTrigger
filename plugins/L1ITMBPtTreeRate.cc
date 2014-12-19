// 
// Class: L1ITMBPtTreeRate 
//
// Info:  Info: Create a tree with the pt values of phi bend in, phi bend out and delta phi of a sample for rate studies
//
// Author: Gian Luca Pinna Angioni
//

#include <stdlib.h>

#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrackFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/ChambPairId.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/PtProvider.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace L1ITMu;
 

typedef struct {Float_t PtBIn,PtBOut,PtDPhi,GMTPt,Wheel_,Qual1_,Qual2_;} PTs_;

// --------------------------------------------------
// Analyzer class, 
// performs tree filling,
// --------------------------------------------------

class L1ITMBPtTreeRate : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtTreeRate( const edm::ParameterSet & );
  ~L1ITMBPtTreeRate();

  //void endJob();
  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);

private:

  TH1F *_hDist; 
  TH1F *_hRateDPhi; 
  TH1F *_hRateDttf; 

  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;
 
  std::map<int, int> NGen;

  TTree *theTree;
  PTs_  Pts;  

};


L1ITMBPtTreeRate::L1ITMBPtTreeRate(const edm::ParameterSet& p) : 
  _mbTracksTag( p.getParameter<edm::InputTag>("MBTracksCollection") ),
  _genParticlesTag( p.getParameter<edm::InputTag>("GenParticlesCollection"))

{

  theTree = _fs->make<TTree>("ElderTree","ElderTree");
  theTree->Branch("Pts" , &Pts,"PtBIn:PtBOut:PtDPhi:GMTPt:Wheel_:Qual1_:Qual2_");

}
 


float L1ITMBPtTreeRate::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}


float L1ITMBPtTreeRate::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);  
}
 
L1ITMBPtTreeRate::~L1ITMBPtTreeRate()
{

}

// void L1ITMBPtTreeRate::endJob()
// {
//   histo->draw();
// }

void L1ITMBPtTreeRate::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtTreeRate]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
    return;
  }

  L1ITMu::MBTrackCollection::const_iterator trIt  = mbTracks->begin();
  L1ITMu::MBTrackCollection::const_iterator trEnd = mbTracks->end();

  for ( ; trIt != trEnd; ++trIt ) {
    
    const L1ITMu::MBTrack & mbTrack = *trIt;

    const std::vector<L1MuRegionalCand> gmtRegCand = mbTrack.getAssociatedGMTin(); 
    
    const L1MuDTTrackCand& dttf = mbTrack.parent();    


    if(dttf.bx() != 0) continue; // CB using only in time DTTF tracks 


    // CB using GMT pt (after sorting) as DTTF has no valid ptValue()
    float gmtInPt = mbTrack.getAssociatedGMTin().size() == 1  ?
                    mbTrack.getAssociatedGMTin()[0].ptValue() : -10.;    


    int sector = dttf.scNum();
    int wheel  = dttf.whNum();

  
    // if (/*wheel!=-1 || */ sector!=0) continue; // CB hack for test just using a given sector
    //if (abs(wheel)==3) continue; // CB hack for test just using a given sector
    
    const L1ITMu::MBLTVectorRef & muonBarrelPrimitives = mbTrack.getStubs();

    L1ITMu::MBLTVectorRef::const_iterator mbPrimIt  = muonBarrelPrimitives.begin();
    L1ITMu::MBLTVectorRef::const_iterator mbPrimEnd = muonBarrelPrimitives.end();

    const L1ITMu::TriggerPrimitive * dtBestPrims[2] = { 0, 0 }; //CB using only MB1 and MB2
     

    // bool first = 1; // GP Add all chambers combinations 
    //  Int_t chamber[2]{0,0}; 

    for ( ; mbPrimIt != mbPrimEnd; ++mbPrimIt ) {
    
      const L1ITMu::TriggerPrimitiveRef & dtMatch = (*mbPrimIt).second;

      int iSt = dtMatch->getDTData().station - 1;
      if ( iSt > 1) continue; //CB using only MB1 and MB2
      if (dtMatch->getDTData().qualityCode == -1) {
	std::cout << "[L1ITMBPtTreeRate]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
	continue;
      }

      dtBestPrims[iSt] = dtMatch.get();
      
      //  std:: cout<<"dtbestPrimis "<<dtBestPrims[iSt]<<std::endl;  //GP Add all chambers combinations. Save chambers informationnumbers in a vector. 
      // std::cout<<"ist "<<iSt<<std::endl;  
      //   if(first==1){
      // 	chamber[0]=iSt;
      // 	first=0;
      // 	  }
      //       else chamber[1]=iSt;
    }
    
    
    std::stringstream PlaceName;
    PlaceName << "Wh" << wheel << "Sc" << sector;
    
    
    if (dtBestPrims[0] && dtBestPrims[1]) { // has MB1 and MB2      
      
      int mb1Obj = objFromPrim(dtBestPrims[0]);
      int mb2Obj = objFromPrim(dtBestPrims[1]);
      
      NGen[1]=1;
      
      if (gmtInPt >= 0.){ 	
	
	Float_t DeltaPhi =  GetDeltaPhi(dtBestPrims[0],dtBestPrims[1]);  
	Float_t PhiBendIn =  GetPhiBend(dtBestPrims[0]); 
	Float_t PhiBendOut =  GetPhiBend(dtBestPrims[1]);   
	
	Float_t ptBendIn = -1;
	Float_t ptBendOut = -1;
	Float_t ptDeltaPhi = -1;
	
	std::string BendInName = "PhiBendIn"+PlaceName.str();
	ptBendIn =  PtProvider::ptLut(BendInName.c_str(), PhiBendIn, mb1Obj,mb2Obj);    
	
	std::string BendOutName = "PhiBendOut"+PlaceName.str();
	ptBendOut =  PtProvider::ptLut(BendOutName.c_str(), PhiBendOut, mb1Obj,mb2Obj);    
	
	std::string DeltaPhiName = "Phi"+PlaceName.str();
	ptDeltaPhi =  PtProvider::ptLut(DeltaPhiName.c_str(), DeltaPhi, mb1Obj,mb2Obj);    	    

	  Pts.PtBIn  = ptBendIn;
	  Pts.PtBOut = ptBendOut;
	  Pts.PtDPhi = ptDeltaPhi;
	  Pts.GMTPt = gmtInPt;
	  Pts.Wheel_ = wheel;
	  Pts.Qual1_ = mb1Obj;
	  Pts.Qual2_ = mb2Obj;

	  theTree->Fill();	
      }	  
    } 
  }
}

					  
int L1ITMBPtTreeRate::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
{

  int qual = prim->getDTData().qualityCode;

  if (qual >=4)
    return ChambPairId::DTCORR;
  else if (qual == 0 || qual == 2)
    return ChambPairId::DTIN;
  else if (qual == 1 || qual == 3)
    return ChambPairId::DTOUT;

  return ChambPairId::NONE;

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtTreeRate);
