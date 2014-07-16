// 
// Class: L1ITMBPtRes
//
// Info: Performs Resolution plots and resolution scatter plots of phi bend and Delta phi variables.
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
#include "TProfile.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace L1ITMu;

// --------------------------------------------------
// Class holding plot (book, fill, draw and save)
// methods for a given chamb object pair
// --------------------------------------------------


class ChPairPlotterRes {
  
public:

  ChPairPlotterRes(int wh, int sec, int inCh, int outCh, 
		   int inChObj, int outChObj, TFileService * fs);

  ChPairPlotterRes(ChambPairId id, TFileService * fs);

  ChPairPlotterRes(const ChPairPlotterRes & plotter);

  int rawId() const { return _id.rawId(); };
  std::string name() const { return _id.name(); };

  const TH1 * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void draw() const ;

  void book(TFileService * fs);

  void fillPtLut(const L1ITMu::TriggerPrimitive * in, 
		 const L1ITMu::TriggerPrimitive * out, float pt); 

  void fillResolution(float dttfPt, float phiBInPt, float phiBOutpt, float DeltaphiPt, float pt,  Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);
 

private:

  ChambPairId _id;

  std::vector<int> _PtResCuts;
  std::map<std::string, TH1 *>  _hPlots;

};


ChPairPlotterRes::ChPairPlotterRes(int wh, int sec, int inCh, int outCh, 
				   int inChObj, int outChObj, TFileService * fs) :
  _id(wh,sec,inCh,outCh,inChObj,outChObj) 
{ 
  
    _PtResCuts.push_back(30); 
  _PtResCuts.push_back(35);
 _PtResCuts.push_back(40);
 _PtResCuts.push_back(45);


  book(fs); 
  
}


ChPairPlotterRes::ChPairPlotterRes(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  // GP gen pt cuts for resolution plots 
  _PtResCuts.push_back(30);
  _PtResCuts.push_back(35);
  _PtResCuts.push_back(40);
  _PtResCuts.push_back(45);

  book(fs); 
  
}

  
  ChPairPlotterRes::ChPairPlotterRes(const ChPairPlotterRes & plotter) :
    _id(plotter._id),_PtResCuts(plotter._PtResCuts),
    _hPlots(plotter._hPlots)
  {
    
  }
  
 
//GP Method for fill scatter plot for correlation studies 
 
  void ChPairPlotterRes::fillPtLut(const L1ITMu::TriggerPrimitive * in, 
				   const L1ITMu::TriggerPrimitive * out, float pt) 
  { 
    
    float inPhiValue  = in->getCMSGlobalPhi();  
    float outPhiValue = out->getCMSGlobalPhi();
    float inPhiBendValue  = fabs(in->getDTData().bendingAngle);
    float outPhiBendValue = fabs(out->getDTData().bendingAngle);
    
    float deltaPhi     = fabs(inPhiValue-outPhiValue); // CB for now it is fabs, change to mu+ mu- 
    
    _hPlots["phiBendInvsphiBendOut"]->Fill(outPhiBendValue,inPhiBendValue);
    _hPlots["phiBendInvsdPhi"]->Fill(deltaPhi,inPhiBendValue);
    _hPlots["phiBendOutvsdPhi"]->Fill(deltaPhi,outPhiBendValue);
    
  }
  
  

void ChPairPlotterRes::fillResolution(float dttfPt, float phiBInPt, float phiBOutPt, float DeltaphiPt, float pt, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName) 
{ 
  
  
  _hPlots["hDeltaPhiPtvsphiBOutPt"] ->Fill(DeltaphiPt,phiBOutPt);
  
  _hPlots["hDeltaPhiPtvsphiBInPt"] ->Fill(DeltaphiPt,phiBInPt);
  
  _hPlots["hGMTPtResolvsPt"] ->Fill(pt,(pt - dttfPt)/ pt);
  
  _hPlots["phiBInPtResolvsPt"] ->Fill(pt,(pt - phiBInPt)/ pt);
  
  _hPlots["phiBOutPtResolvsPt"] ->Fill(pt,(pt - phiBOutPt)/ pt);
  
  _hPlots["DeltaphiPtResolvsPt"] ->Fill(pt,(pt - DeltaphiPt)/ pt);
  
  
  std::vector<int>::const_iterator PtResCutIt = _PtResCuts.begin();
  std::vector<int>::const_iterator PtResCutEnd = _PtResCuts.end();
  
  
  
  for (; PtResCutIt != PtResCutEnd; ++PtResCutIt)
    {
      int PtResCut = (*PtResCutIt);
      
      if(pt<=PtResCut){
	std::stringstream ptCut; ptCut << PtResCut;
      std::string PtTag = ptCut.str();
   
      PtTag="Cut_"+PtTag;

      _hPlots["GMTptResol" + PtTag]->Fill((pt - dttfPt)/ pt);
      
      _hPlots["PhiBInptResol" + PtTag]->Fill((pt - phiBInPt)/ pt);
      
      _hPlots["PhiBOutptResol" + PtTag]->Fill((pt - phiBOutPt)/ pt);
      
      _hPlots["DeltaPhiptResol" + PtTag]->Fill((pt - DeltaphiPt)/ pt);
      

      _hPlots["DeltaBInPtResol"+PtTag]->Fill(((pt - DeltaphiPt)/pt),((pt - phiBInPt)/pt));
      _hPlots["DeltaBOutPtResol"+PtTag]->Fill(((pt - DeltaphiPt)/pt),((pt - phiBOutPt)/pt));
      _hPlots["BInBOutPtResol"+PtTag]->Fill(((pt - phiBInPt)/pt),((pt - phiBOutPt)/pt));
      
      }
  }
}


void ChPairPlotterRes::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
  
 
   
  TFileDirectory baseFolder = fs->mkdir(hDir.c_str());
      

  //GP Scatter plots for correlation studies

   
  _hPlots["phiBendInvsphiBendOut" ] = baseFolder.make<TH2F>(("hdPhiBendInvsphiBendOut"  + hName).c_str(), 
							    ("obj #Delta#phi_{b}In vs #Delta#phi_{b}Out for " + hName +
							     ";#Delta#phi_{b} Out;#Delta#phi_{b} In").c_str(),
							    100, 0, 100, 100, 0, 100);
  
  _hPlots["phiBendInvsdPhi" ] = baseFolder.make<TH2F>(("hdPhiBendInvsdPhi"  + hName).c_str(), 
						      ("obj #Delta#phi_{b}In vs dPhi for " + hName +
						       ";#Delta#phi;#Delta#phi_{b} In").c_str(),
						      100, -0, .05, 100, 0, 100);
  
  _hPlots["phiBendOutvsdPhi" ] = baseFolder.make<TH2F>(("hdPhiBendOutvsdPhi"  + hName).c_str(), 
						       ("obj #Delta#phi_{b}Out vs dPhi for " + hName +
  						";#Delta#phi;#Delta#phi_{b} Out").c_str(),		
						100, -0, .05, 100, 0, 100);
  

  

  TFileDirectory folderRes = baseFolder.mkdir("Resolutions");
  

  _hPlots["hGMTPtResolvsPt" ] = folderRes.make<TH2F>(("hGMTPtResolvsPt" + hName ).c_str(), 
  					     (" obj (dttf_[pt] - pt) / pt vs pt  " + hName + 
  					      "; p_{T}; DTTF  mu p_{T}").c_str(), 
  					    60, 0., 120., 120, -2.5,2.5);
  


  _hPlots["hDeltaPhivsGMTPt"]  = folderRes.make<TH2F>(("hDeltaPhivsGMTPt" + hName ).c_str(), 
  					     (" Delta phi vs dttf_[pt] " + hName + 
  					      " ;#Delta #phi; DTTF  mu p_{T}").c_str(), 
  					      280, 0.,140., 100, 0.,.05);



  _hPlots["hGMTPtResolvsDeltaPhiPt" ] = folderRes.make<TH2F>(("hGMTPtResolvsDeltaPhiPt" + hName ).c_str(), 
  					     ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
  					      ";GEN  mu p_{T};(GEN  mu p_{T} - DTTF  mu p_{T}) / GEN  mu p_{T}").c_str(), 
  					     60, 0.,120., 120, -2.5, 2.5);
  
  



  _hPlots["hDeltaPhiPtvsphiBInPt"]  = folderRes.make<TH2F>(("hDeltaPhiPtvsphiBInPt" + hName ).c_str(), 
							(" #Delta #phi mu p_{T} vs #phi B In p_{t} " + hName + 
							 "; #Delta #phi mu  p_{T}; #phi B In mu p_{T}").c_str(), 
							60, 0.,120., 60, 0.,120.);



  _hPlots["hDeltaPhiPtvsphiBOutPt"]  = folderRes.make<TH2F>(("hDeltaPhiPtvsphiBOutPt" + hName ).c_str(), 
							(" #Delta #phi mu p_{T} vs #phi B Out p_{t} " + hName + 
							 "; #Delta #phi mu  p_{T}; #phi B Out mu p_{T}").c_str(), 
							60, 0.,120., 60, 0.,120.);

  
  _hPlots["phiBInPtResolvsPt"] = folderRes.make<TH2F>(("hPhiBendInPtResolvsPt" + hName ).c_str(), 
						      ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
						       ";GEN  mu p_{T};(GEN  mu p_{T} - #phi_{b}In  mu p_{T}) / GEN  mu p_{T}").c_str(), 
						      60, 0.,120., 120, -2.5, 2.5);
  
  _hPlots["phiBOutPtResolvsPt"] = folderRes.make<TH2F>(("hPhiBendOutPtResolvsPt" + hName ).c_str(), 
						       ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							";GEN  mu p_{T};(GEN  mu p_{T} - #phi_{b}Out  mu p_{T}) / GEN  mu p_{T}").c_str(), 
						       60, 0.,120., 120, -2.5, 2.5);
  
  _hPlots["DeltaphiPtResolvsPt"] = folderRes.make<TH2F>(("hDeltaPhiPtResolvsPt" + hName ).c_str(), 
							("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							 ";GEN  mu p_{T};(GEN  mu p_{T} - #Delta #phi  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							60, 0.,120., 120, -2.5, 2.5);


 std::vector<int>::const_iterator PtResCutIt = _PtResCuts.begin();
  std::vector<int>::const_iterator PtResCutEnd = _PtResCuts.end();


 for (; PtResCutIt != PtResCutEnd; ++PtResCutIt)
    {
      int PtResCut = (*PtResCutIt);

      std::stringstream ptCut; ptCut << PtResCut;
      std::string PtTag = ptCut.str();
      PtTag="Cut_"+PtTag;

      
    TFileDirectory folderCut = folderRes.mkdir(PtTag.c_str());
    
    _hPlots["DeltaBInPtResol"+PtTag] = folderCut.make<TH2F>(("DeltaBInPtResol" + hName+PtTag ).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt " + hName + PtTag+
							     ";(GEN  mu p_{T} - #Delta #phi pt mu p_{T}) / GEN  mu p_{T};(GEN  mu p_{T} - #phi B In pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -2.5, 2.5, 120, -2.5, 2.5);
    
    _hPlots["DeltaBOutPtResol"+PtTag] = folderCut.make<TH2F>(("DeltaBOutPtResol" + hName+PtTag ).c_str(), 
							     ("obj (dttf_[pt] - pt) / pt " + hName + PtTag+
							      ";(GEN  mu p_{T} - #Delta #phi pt mu p_{T}) / GEN  mu p_{T};(GEN  mu p_{T} - #phi B Out pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
							     120, -2.5, 2.5, 120, -2.5, 2.5);
    
    
    _hPlots["BInBOutPtResol"+PtTag] = folderCut.make<TH2F>(("BInBOutPtResol" + hName+PtTag ).c_str(), 
							   ("obj (dttf_[pt] - pt) / pt " + hName + PtTag+
							    ";(GEN  mu p_{T} - #phi B In pt mu p_{T}) / GEN  mu p_{T};(GEN  mu p_{T} - #phi B Out pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
							   120, -2.5, 2.5, 120, -2.5, 2.5);
    
    
    
      _hPlots["GMTptResol" + PtTag] = folderCut.make<TH1F>(("hGMTPtResol" + hName + PtTag).c_str(), 
    							 ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
    							  ";(GEN  mu p_{T} - DTTF  mu p_{T}) / GEN  mu p_{T}").c_str(), 
    							 120, -5., 1.);
    

         
    _hPlots["DeltaPhiptResol" + PtTag] = folderCut.make<TH1F>(("hDeltaPhiPtResol" + hName + PtTag).c_str(), 
							      ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							       ";(GEN  mu p_{T} - #Delta #phi  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							      120, -5., 1.);
    
    
    _hPlots["PhiBInptResol" + PtTag] = folderCut.make<TH1F>(("hPhiBInPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - #phi_{b}In  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);
    
    _hPlots["PhiBOutptResol" + PtTag] = folderCut.make<TH1F>(("hPhiBOutPtResol" + hName + PtTag).c_str(), 
							     ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							      ";(GEN  mu p_{T} - #phi_{b}Out  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							     120, -5., 1.);
    }
}


void ChPairPlotterRes::draw() const
{

  gStyle->SetPalette(52);
  
  
  std::map<std::string, TH1 *>::const_iterator hPlotsIt = _hPlots.begin();
  std::map<std::string, TH1 *>::const_iterator hPlotsEnd = _hPlots.end();

  for (; hPlotsIt != hPlotsEnd ; ++hPlotsIt)
    {

      TH1 * hHisto = hPlotsIt->second;

      if ( hPlotsIt->first.find("ptResol") == std::string::npos ) continue;
 
      std::string histoName = hPlotsIt->first;
      
      std::string cName = "c" + histoName + name();
      
      system(std::string("mkdir -p plots/Resolution/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();
      
      hHisto->Draw("");

      c->SaveAs(("plots/Resolution/" + dttfId().name() + "/" + cName+".pdf").c_str());
  
    }

}




// --------------------------------------------------
// Analyzer class, holds all ChPairPlotterRess 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtRes : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtRes( const edm::ParameterSet & );
  ~L1ITMBPtRes();

  void endJob();
  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  ChPairPlotterRes* getPlotter(int dttfRawId, int chPairRawId); 


private:


  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChPairPlotterRes*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

  std::map<int, int> NGen;
};


L1ITMBPtRes::L1ITMBPtRes(const edm::ParameterSet& p) : 
  _mbTracksTag( p.getParameter<edm::InputTag>("MBTracksCollection") ),
  _genParticlesTag( p.getParameter<edm::InputTag>("GenParticlesCollection") )

{

  // CB: combinations to be plotted go here 
  // for now just booking DT IN, OUT CORR
  // and ch 1 and 2 of DTTF wh 1 sec 1
  
  _mbPtChambObjects.push_back(ChambPairId::DTIN);
  _mbPtChambObjects.push_back(ChambPairId::DTOUT);
  _mbPtChambObjects.push_back(ChambPairId::DTCORR);

 

 for (int wheel = -3; wheel <=3; ++wheel) {
    if (wheel == 0) continue;
    for (int sector = 0; sector <=0; ++sector) { // Change for all Sector
      for (int inCh = 1; inCh <=1; ++inCh) { //Change for all Chambers
	for (int outCh = inCh + 1; outCh <=2; ++outCh) { //Change for all Chambers
	  
	  std::vector<int>::const_iterator inObjIt = _mbPtChambObjects.begin();
	  std::vector<int>::const_iterator objEnd  = _mbPtChambObjects.end();
	  
	  for (; inObjIt != objEnd; ++inObjIt) {
	    std::vector<int>::const_iterator outObjIt = _mbPtChambObjects.begin();
	    for (; outObjIt != objEnd; ++outObjIt) {
	      
	      ChambPairId id(wheel, sector, inCh, outCh, (*inObjIt), (*outObjIt));
	      ChPairPlotterRes *plotter = new ChPairPlotterRes(id,_fs.operator->());  

	      histos[DTTFId(wheel,sector).rawId()][id.rawId()] = plotter;
	      
	    }
	  }
	}
      }
    }
  }  
}
 


float L1ITMBPtRes::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}




float L1ITMBPtRes::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);
  
}


L1ITMBPtRes::~L1ITMBPtRes()
{

}

void L1ITMBPtRes::endJob()
{
    
  std::map<int,std::map<int,ChPairPlotterRes*> >::iterator dttfIt = histos.begin();
  std::map<int,std::map<int,ChPairPlotterRes*> >::iterator dttfEnd = histos.end();
  
  for(; dttfIt!=dttfEnd; ++dttfIt) {
    
    std::map<int,ChPairPlotterRes*>::iterator chambPairIt = dttfIt->second.begin();
    std::map<int,ChPairPlotterRes*>::iterator chambPairEnd = dttfIt->second.end();
    
    for(; chambPairIt!=chambPairEnd; ++chambPairIt) {
      
      ChPairPlotterRes *plotter = chambPairIt->second;

      plotter->draw();
      delete plotter;
      
    }
  }
  
}

void L1ITMBPtRes::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{


  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtRes]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
    return;
  }

  L1ITMu::MBTrackCollection::const_iterator trIt  = mbTracks->begin();
  L1ITMu::MBTrackCollection::const_iterator trEnd = mbTracks->end();

  for ( ; trIt != trEnd; ++trIt ) {
    
    const L1ITMu::MBTrack & mbTrack = *trIt;

    const std::vector<L1MuRegionalCand> gmtRegCand = mbTrack.getAssociatedGMTin(); 
    
    const reco::GenParticle *bestGen=0;

    if ( gmtRegCand.size() == 1) {
      
      float gmtPhi = gmtRegCand.at(0).phiValue();
      float gmtEta = gmtRegCand.at(0).etaValue();

      float bestDR = 999.;

      reco::GenParticleCollection::const_iterator genIt  = gens->begin();
      reco::GenParticleCollection::const_iterator genEnd = gens->end();
      
      for(; genIt!=genEnd; ++genIt) {
	if (abs(genIt->pdgId()) == 13 ) {
	  
	  float genPhi = genIt->phi();
	  float genEta = genIt->eta();
	  
	  float dR = deltaR(gmtEta,gmtPhi,genEta,genPhi);
		  
	  if (dR < 1. && dR < bestDR) { // CB get it from CFG 
	    bestDR = dR;
	    bestGen = &(*genIt);
	  }
	}	
      }
    }
    
    if (!bestGen) continue; //CB skip if GEN matching is missing or != 1 GMT in matched

    const L1MuDTTrackCand& dttf = mbTrack.parent();    
    if(dttf.bx() != 0) continue; // CB using only in time DTTF tracks 

    // CB using GMT pt (after sorting) as DTTF has no valid ptValue()
    float gmtInPt = mbTrack.getAssociatedGMTin().size() == 1  ?
                    mbTrack.getAssociatedGMTin()[0].ptValue() : -10.;    

    int sector = dttf.scNum();
    int wheel  = dttf.whNum();

    if (/*wheel!=-1 || */ sector!=0) continue; // CB hack for test just using a given sector
    
    int dttfRawId  = DTTFId(wheel,sector).rawId();
    
    const L1ITMu::MBLTVectorRef & muonBarrelPrimitives = mbTrack.getStubs();

    L1ITMu::MBLTVectorRef::const_iterator mbPrimIt  = muonBarrelPrimitives.begin();
    L1ITMu::MBLTVectorRef::const_iterator mbPrimEnd = muonBarrelPrimitives.end();

    const L1ITMu::TriggerPrimitive * dtBestPrims[2] = { 0, 0 }; //CB using only MB1 and MB2
    
    for ( ; mbPrimIt != mbPrimEnd; ++mbPrimIt ) {
      
      
      const L1ITMu::TriggerPrimitiveRef & dtMatch = (*mbPrimIt).second;
      
      int iSt = dtMatch->getDTData().station - 1;
      if ( iSt > 1) continue; //CB using only MB1 and MB2
      if (dtMatch->getDTData().qualityCode == -1) {
	std::cout << "[L1ITMBPtRes]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
	continue;
      }
      
      dtBestPrims[iSt] = dtMatch.get();
      
    }
    
    std::stringstream PlaceName;
    PlaceName << "Wh" << wheel << "Sc" << sector;
    
    if (dtBestPrims[0] && dtBestPrims[1]) { // has MB1 and MB2
      
      int mb1Obj = objFromPrim(dtBestPrims[0]);
      int mb2Obj = objFromPrim(dtBestPrims[1]);
      
      int chPairRawId = ChambPairId(wheel,sector,1,2,mb1Obj,mb2Obj).rawId();
      
      ChPairPlotterRes *plotter = getPlotter(dttfRawId,chPairRawId);
      
      NGen[1]=1;


      if (plotter && gmtInPt >= 0.){
	
	Float_t GenPt  = bestGen->pt();
	
	Float_t DeltaPhi =  GetDeltaPhi(dtBestPrims[0],dtBestPrims[1]);  
	Float_t PhiBendIn =  GetPhiBend(dtBestPrims[0]); 
	Float_t PhiBendOut =  GetPhiBend(dtBestPrims[1]);   
	
	Float_t ptBendIn = -1;
	Float_t ptBendOut = -1;
	Float_t ptDeltaPhi = -1;
	
	std::string BendInName = "PhiBendIn"+PlaceName.str();
	ptBendIn =  PtProvider::ptLut(BendInName.c_str(), PhiBendIn, mb1Obj,mb2Obj);    
	
	std::string BendOutName = "PhiBendOut"+PlaceName.str();
	ptBendOut = PtProvider::ptLut(BendOutName.c_str(), PhiBendOut, mb1Obj,mb2Obj);    
	
	std::string DeltaPhiName = "Phi"+PlaceName.str();
	ptDeltaPhi = PtProvider::ptLut(DeltaPhiName.c_str(), DeltaPhi, mb1Obj,mb2Obj); 
	
	
	
	std::string QualInString, QualOutString;
	
	switch(mb1Obj){
	  
	case 0: QualInString = "DTIN";
	  break;
	case 1: QualInString = "DTCORR";
	  break;
	case 2: QualInString = "DTDIRR";
	  break;
	case 3: QualInString = "DTOUT";
	  break;
	case 4: QualInString = "NONE";
	  break;
	}
	
	switch(mb2Obj){
	  
	case 0: QualOutString = "DTIN";
	  break;
	case 1: QualOutString = "DTCORR";
	  break;
	case 2: QualOutString = "DTDIRR";
	  break;
	case 3: QualOutString = "DTOUT";
	  break;
	case 4: QualOutString = "NONE";
	  break;
	}
	
	
	//GP print information for events with bad resolution 
	
	//  if(((GenPt - ptDeltaPhi)/GenPt)<-1) std::cout<<"Location  "<<PlaceName.c_str()<<" "<<(QualInString+QualOutString).c_str()<<std::endl<<std::endl<<"Pt Gen "<<GenPt<<std::endl<<"GMT pt "<<gmtInPt<<std::endl<<"Delta phi Pt "<<ptDeltaPhi<<"  Delta Phi "<<DeltaPhi<<std::endl<<"Phi B In Pt "<<ptBendIn<<"  Phi Bend In "<<PhiBendIn<<std::endl<<"Phi B Out Pt "<<ptBendOut<<"   Phi Ben In "<<PhiBendOut<<std::endl<<std::endl;;
	

       	
	//GP Fill scatter plot for correlation studies
	
	  plotter->fillPtLut(dtBestPrims[0],
				     dtBestPrims[1],
				     GenPt); 
	
	
	plotter->fillResolution(gmtInPt, ptBendIn, ptBendOut, ptDeltaPhi, GenPt,mb1Obj,mb2Obj,PlaceName.str()); 	          
	
      }
    }
  }
}


int L1ITMBPtRes::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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


ChPairPlotterRes* L1ITMBPtRes::getPlotter(int dttfRawId, int chPairRawId) 
{
  // CB optimise 
  if (histos.find(dttfRawId) == histos.end()) { 
    std::cout << "[L1ITMBPtRes]::getPlotter : DTTF not in histos. Return 0.\n";
    return 0;
  }
  if (histos[dttfRawId].find(chPairRawId) == histos[dttfRawId].end()) {
    std::cout << "[L1ITMBPtRes]::getPlotter : ChambPair not in histos. Return 0.\n";
    return 0;
  }
  return histos[dttfRawId][chPairRawId];
}




#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtRes);
