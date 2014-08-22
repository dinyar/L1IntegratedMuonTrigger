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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrack.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBTrackFwd.h"


#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

using namespace L1ITMu;

typedef edm::ParameterSet PSet;
int nptbin=50;
class L1ITMBDTTFStudies : public edm::EDAnalyzer {
public:
  L1ITMBDTTFStudies( const PSet& );
  ~L1ITMBDTTFStudies();
  void analyze( const edm::Event&, const edm::EventSetup& );
  void findBestGmtDt(
				      reco::GenParticleCollection::const_iterator& genIt,
				      edm::Handle<L1ITMu::MBTrackCollection>& mbTracks,
				      float& phi,
				      float& eta,
				      float& pt,
				      int& quality,
				      float& dr
				      );

private:
  void endJob();
  void beginJob();
  edm::InputTag _mbTracksTag;
  edm::InputTag _mbTracksTag2;
  edm::InputTag _genParticlesTag;
  edm::Service<TFileService> _fs;
  TH1F * _bxGmtDt;
  TH1F * _phiGen;  TH1F * _etaGen;  TH1F * _ptGen; TH2F * _phiVsEtaGen;
  TH1F * _phiGenNum;  TH1F * _etaGenNum;  TH1F * _ptGenNum; TH2F * _phiVsEtaGenNum;
  TH1F * _phiGenNumBx0;  TH1F * _etaGenNumBx0;  TH1F * _ptGenNumBx0; TH2F * _phiVsEtaGenNumBx0;
  TH1F * _phiGenNumBx0T2;  TH1F * _etaGenNumBx0T2;  TH1F * _ptGenNumBx0T2; TH2F * _phiVsEtaGenNumBx0T2;
  

  TH1F * _phiTrack;  TH1F * _etaTrack;  TH1F * _ptTrack;  TH1F * _qualTrack;
  TH1F * _dPhiTrackGen;  TH1F * _dEtaTrackGen;   TH1F * _dRTrackGen;  TH2F * _dPhiDEtaTrackGen;
  TH2F * _phiTrackVsPhiGen; TH2F * _etaTrackVsEtaGen; TH2F * _dPhiTrackGenVsPtGen;
  TH1F * _dPtTrackGen;
  
  TH2F * _ptGenVsEtaGen; TH2F * _ptGenVsEtaGenNumBx0T1; TH2F * _ptGenVsEtaGenNumBx0T2;
  TH2F * _ptGenVsEtaGenNumBx0T1NotT2; TH2F * _ptGenVsEtaGenNumBx0T2NotT1;
  TH2F * _ptGenVsEtaGenNumBx0T1AndT2;

  TH2F * _qualT2vsQualT1; TH2F * _ptT2vsPtT1;  TH2F * _dptT2vsDptT1;
};
double pig = acos(-1);
double phistep=pig/144.0;
L1ITMBDTTFStudies::L1ITMBDTTFStudies(const PSet& p):
  _mbTracksTag(p.getParameter<edm::InputTag>("mbTracksCollection")),
  _mbTracksTag2(p.getParameter<edm::InputTag>("mbTracksCollection2")),
  _genParticlesTag(p.getParameter<edm::InputTag>("genParticlesCollection")){

  _bxGmtDt = _fs->make<TH1F>("bxGmtDt","bxGmtDt",7,-3.5,3.5);
  
  _phiGen = _fs->make<TH1F>("phiGen","phiGen",144,-pig,pig);
  _etaGen = _fs->make<TH1F>("etaGen","etaGen",64,-1.2,1.2);
  _ptGen = _fs->make<TH1F>("ptGen","ptGen",nptbin,0,100);
  _phiVsEtaGen = _fs->make<TH2F>("phiVsEtaGen","phiVsEtaGen",64,-1.2,1.2,144,-pig,pig);

  _phiGenNum = _fs->make<TH1F>("phiGenNum","phiGenNum",144,-pig,pig);
  _etaGenNum = _fs->make<TH1F>("etaGenNum","etaGenNum",64,-1.2,1.2);
  _ptGenNum = _fs->make<TH1F>("ptGenNum","ptGenNum",nptbin,0,100);
  _phiVsEtaGenNum = _fs->make<TH2F>("phiVsEtaGenNum","phiVsEtaGenNum",64,-1.2,1.2,144,-pig,pig);

  _phiGenNumBx0 = _fs->make<TH1F>("phiGenNumBx0","phiGenNumBx0",144,-pig,pig);
  _etaGenNumBx0 = _fs->make<TH1F>("etaGenNumBx0","etaGenNumBx0",64,-1.2,1.2);
  _ptGenNumBx0 = _fs->make<TH1F>("ptGenNumBx0","ptGenNumBx0",nptbin,0,100);
  _phiVsEtaGenNumBx0 = _fs->make<TH2F>("phiVsEtaGenNumBx0","phiVsEtaGenNumBx0",64,-1.2,1.2,144,-pig,pig);

  _phiGenNumBx0T2 = _fs->make<TH1F>("phiGenNumBx0T2","phiGenNumBx0T2",144,-pig,pig);
  _etaGenNumBx0T2 = _fs->make<TH1F>("etaGenNumBx0T2","etaGenNumBx0T2",64,-1.2,1.2);
  _ptGenNumBx0T2 = _fs->make<TH1F>("ptGenNumBx0T2","ptGenNumBx0T2",nptbin,0,100);
  _phiVsEtaGenNumBx0T2 = _fs->make<TH2F>("phiVsEtaGenNumBx0T2","phiVsEtaGenNumBx0T2",64,-1.2,1.2,144,-pig,pig);
  
  _phiTrack = _fs->make<TH1F>("phiTrack","phiTrack",144,-pig,pig);
  _etaTrack = _fs->make<TH1F>("etaTrack","etaTrack",64,-1.2,1.2);
  _ptTrack = _fs->make<TH1F>("ptTrack","ptTrack",nptbin,0,100);
  _qualTrack = _fs->make<TH1F>("qualTrack","qualTrack",8,0,8);
  _dPhiTrackGen = _fs->make<TH1F>("dPhiTrackGen","dPhiTrackGen",nptbin,-pig/4.0,pig/4.0);
  _dEtaTrackGen = _fs->make<TH1F>("dEtaTrackGen","dEtaTrackGen",nptbin,-0.5,0.5);
  _dRTrackGen = _fs->make<TH1F>("dRTrackGen","dRTrackGen",nptbin,0,2);
  _dPhiDEtaTrackGen = _fs->make<TH2F>("dPhiDEtaTrackGen","dPhiDEtaTrackGen",nptbin,-0.5,0.5,nptbin,-pig/4.0,pig/4.0);

  _phiTrackVsPhiGen = _fs->make<TH2F>("phiTrackVsPhiGen","phiTrackVsPhiGen",144,-pig,pig,144,-pig,pig);
  _etaTrackVsEtaGen = _fs->make<TH2F>("etaTrackVsEtaGen","etaTrackVsEtaGen",64,-1.2,1.2,64,-1.2,1.2);
  _dPhiTrackGenVsPtGen = _fs->make<TH2F>("dPhiTrackGenVsPtGen","dPhiTrackGenVsPtGen",nptbin,0,100,nptbin,-pig/4.0,pig/4.0);

  _ptGenVsEtaGen = _fs->make<TH2F>("ptGenVsEtaGen","ptGenVsEtaGen",64,-1.2,1.2, nptbin,0,100);
  _ptGenVsEtaGenNumBx0T1 = _fs->make<TH2F>("ptGenVsEtaGenNumBx0T1","ptGenVsEtaGenNumBx0T1",64,-1.2,1.2, nptbin,0,100);
  _ptGenVsEtaGenNumBx0T2 = _fs->make<TH2F>("ptGenVsEtaGenNumBx0T2","ptGenVsEtaGenNumBx0T2",64,-1.2,1.2, nptbin,0,100);
  _ptGenVsEtaGenNumBx0T2NotT1 = _fs->make<TH2F>("ptGenVsEtaGenNumBx0T2NotT1","ptGenVsEtaGenNumBx0T2NotT1",64,-1.2,1.2, nptbin,0,100);
  _ptGenVsEtaGenNumBx0T1NotT2 = _fs->make<TH2F>("ptGenVsEtaGenNumBx0T1NotT2","ptGenVsEtaGenNumBx0T1NotT2",64,-1.2,1.2, nptbin,0,100);
  _ptGenVsEtaGenNumBx0T1AndT2 = _fs->make<TH2F>("ptGenVsEtaGenNumBx0T1AndT2","ptGenVsEtaGenNumBx0T1AndT2",64,-1.2,1.2, nptbin,0,100);
  _qualT2vsQualT1 = _fs->make<TH2F>("qualT2vsQualT1","qualT2vsQualT1",8,0,8,8,0,8);
  _ptT2vsPtT1 =  _fs->make<TH2F>("ptT2vsPtT1","ptT2vsPtT1",141,0,141,141,0,141);
  _dptT2vsDptT1 = _fs->make<TH2F>("dptT2vsDptT1","dptT2vsDptT1",nptbin,-3,3,nptbin,-3,3);
}

L1ITMBDTTFStudies::~L1ITMBDTTFStudies(){}
void L1ITMBDTTFStudies::beginJob(){}
void L1ITMBDTTFStudies::endJob(){}

void L1ITMBDTTFStudies::findBestGmtDt(
				      reco::GenParticleCollection::const_iterator& genIt,
				      edm::Handle<L1ITMu::MBTrackCollection>& mbTracks,
				      float& phi,
				      float& eta,
				      float& pt,
				      int& quality,
				      float& dr
				      ){
    if (abs(genIt->pdgId()) == 13 && genIt->eta()>=-10000) {
      L1ITMu::MBTrackCollection::const_iterator trIt = mbTracks->begin();
      L1ITMu::MBTrackCollection::const_iterator trEnd = mbTracks->end();
      float bestDRBx0 = 999.;
      float bestPhi=999.;
      float bestEta=999.;
      float bestPt=999.;
      int bestQual=999;
      float genPhi = genIt->phi();
      float genEta = genIt->eta();
      
      bool hadnegbx=false;
    
      for ( ; trIt != trEnd; ++trIt ) {
	const L1ITMu::MBTrack & mbTrack = *trIt;
	const std::vector<L1MuRegionalCand> gmtRegCand = mbTrack.getAssociatedGMTin();
	if ( gmtRegCand.size() == 1) {
      	  float gmtPhi = gmtRegCand.at(0).phiValue() + phistep;
	  float newGmtPhi = (gmtPhi>pig)?gmtPhi-2.0*pig:gmtPhi;
	  float gmtEta = gmtRegCand.at(0).etaValue();
	  float gmtPt = gmtRegCand.at(0).ptValue();
	  float gmtQual = gmtRegCand.at(0).quality();
	  int gbx = gmtRegCand.at(0).bx();
	  if (gbx==-1 || gbx==-2) hadnegbx=true;
	  float dR = reco::deltaR(gmtEta,newGmtPhi,genEta,genPhi);
	  if (dR < 1. && dR < bestDRBx0 && gbx==0 && !hadnegbx) { // CB get it from CFG
	    bestDRBx0 = dR;
	    bestPhi=newGmtPhi;
	    bestEta=gmtEta;
	    bestPt=gmtPt;
	    bestQual=gmtQual;
	  }
	}
      }
      phi=bestPhi;
      eta=bestEta;
      pt=bestPt;
      quality=bestQual;
      dr=bestDRBx0;
    }
    return;
}

void L1ITMBDTTFStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks2;
  iEvent.getByLabel( _mbTracksTag2, mbTracks2);
  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);
  if (!gens.isValid() || !mbTracks.isValid() || !mbTracks2.isValid()) {
    std::cout << "[L1ITMBPtStudies]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
    return;
  }
  reco::GenParticleCollection::const_iterator genIt = gens->begin();
  reco::GenParticleCollection::const_iterator genEnd = gens->end();
  for(; genIt!=genEnd; ++genIt) {
    if (abs(genIt->pdgId()) == 13 && genIt->eta()>=0) {
      float genPhi = genIt->phi();
      float genEta = genIt->eta();
      float genPt = genIt->pt();
      if (genPt>10)_phiGen->Fill(genPhi);
      if (genPt>10)_etaGen->Fill(genEta);
      _ptGen->Fill(genPt);
      _phiVsEtaGen->Fill(genEta,genPhi);
     

      //special things once
      L1ITMu::MBTrackCollection::const_iterator trIt = mbTracks->begin();
      L1ITMu::MBTrackCollection::const_iterator trEnd = mbTracks->end();
      float bestDR = 999.;
      float bestDRBx0 = 999.;
      //bool hadnegbx=false;
      for ( ; trIt != trEnd; ++trIt ) {
	const L1ITMu::MBTrack & mbTrack = *trIt;
	const std::vector<L1MuRegionalCand> gmtRegCand = mbTrack.getAssociatedGMTin();
	if ( gmtRegCand.size() == 1) {
      	  float gmtPhi = gmtRegCand.at(0).phiValue() + phistep;
	  float newGmtPhi = (gmtPhi>pig)?gmtPhi-2.0*pig:gmtPhi;
	  float gmtEta = gmtRegCand.at(0).etaValue();
	  float gmtPt = gmtRegCand.at(0).ptValue();
	  float gmtQual = gmtRegCand.at(0).quality();
	  int gbx = gmtRegCand.at(0).bx();
	  _bxGmtDt->Fill(gbx);
	  //if (gbx==-1 || gbx==-2) hadnegbx=true;
	  _phiTrack->Fill(newGmtPhi); _etaTrack->Fill(gmtEta);_ptTrack->Fill(gmtPt);_qualTrack->Fill(gmtQual);
       	  float dR = reco::deltaR(gmtEta,newGmtPhi,genEta,genPhi);
	  float dPhi = reco::deltaPhi(newGmtPhi,genPhi);
	  float dEta = gmtEta - genEta;
	  _dPhiTrackGen->Fill(dPhi);_dEtaTrackGen->Fill(dEta);_dPhiDEtaTrackGen->Fill(dEta,dPhi);
	  _phiTrackVsPhiGen->Fill(genPhi,newGmtPhi);_etaTrackVsEtaGen->Fill(genEta,gmtEta);_dPhiTrackGenVsPtGen->Fill(genPt,dPhi);
	  if (dR < 1. && dR < bestDR) { // CB get it from CFG
	    bestDR = dR;
	  }
	  if (dR < 1. && dR < bestDRBx0 && gbx==0) { // CB get it from CFG
	    bestDRBx0 = dR;
	  }
	}
      }
      if (bestDR<900) _dRTrackGen->Fill(bestDR);
      if (bestDR<1.){//matched
	if (genPt>10)_phiGenNum->Fill(genPhi);	if (genPt>10)_etaGenNum->Fill(genEta);
	_ptGenNum->Fill(genPt);_phiVsEtaGenNum->Fill(genEta,genPhi);
      }



      //normal things here
      float phi1,eta1,pt1,dr1,phi2,eta2,pt2,dr2;
      int qual1, qual2;
      findBestGmtDt(genIt, mbTracks, phi1, eta1, pt1, qual1, dr1);
      findBestGmtDt(genIt, mbTracks2, phi2, eta2, pt2, qual2, dr2);
      _ptGenVsEtaGen->Fill(genEta,genPt);
      if (dr1<1){// one mbtrack was matched
	_ptGenVsEtaGenNumBx0T1->Fill(genEta,genPt);
	_ptGenNumBx0->Fill(genPt);	
	if (genPt>10){
	  _phiGenNumBx0->Fill(genPhi);	
	  _etaGenNumBx0->Fill(genEta);
	  _phiVsEtaGenNumBx0->Fill(genEta,genPhi);
	}
      }
      if (dr2<1){// one mbtrack2 was matched
	_ptGenVsEtaGenNumBx0T2->Fill(genEta,genPt);
	_ptGenNumBx0T2->Fill(genPt);	
	if (genPt>10){
	  _phiGenNumBx0T2->Fill(genPhi);	
	  _etaGenNumBx0T2->Fill(genEta);
	  _phiVsEtaGenNumBx0T2->Fill(genEta,genPhi);
	}
      }
      if (dr2<1 && dr1>1){// phase1 found track not found by legacy
	_ptGenVsEtaGenNumBx0T2NotT1->Fill(genEta,genPt);
      }
      if (dr1<1 && dr2>1){// legacy found track not found by phase1
	_ptGenVsEtaGenNumBx0T1NotT2->Fill(genEta,genPt);
      }
      if (dr1<1 && dr2<1){// legacy and phase1 both found track
	_ptGenVsEtaGenNumBx0T1AndT2->Fill(genEta,genPt);
	_qualT2vsQualT1->Fill(qual1,qual2);
	_ptT2vsPtT1->Fill(pt1,pt2); 
	_dptT2vsDptT1->Fill((pt1-genPt)/genPt,(pt2-genPt)/genPt);

      }
      
      // todo
      // - study events where mbtrack2 is efficient and mbtrack is not
      // - study the Dpt wrt gen for mbtrack and mbtrack2
      // - study the correlation between mbtrack and mbtrack2 pt assignement
    }
 }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBDTTFStudies);


