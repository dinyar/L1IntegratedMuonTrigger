// 
// Class: L1ITMBPtLutPlots
//
// Info: Performs GEN-DTTF matching and computes LUT plots  
//
// Author: 
//

#include <stdlib.h>

#include <memory>
#include <string>
#include <sstream>
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

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace L1ITMu;


// --------------------------------------------------
// Id class to handle DTTF identification
// --------------------------------------------------

class DTTFId {

public:

  DTTFId(int wh, int sec);
  DTTFId(const DTTFId & id);

  ~DTTFId() { };

  int rawId() const { return (_wh+3) + 10*_sec; } ;
  std::string name() const;

private:
  
  int _wh, _sec;

};


DTTFId::DTTFId(int wh, int sec) : _wh(wh), _sec(sec) 
{
 
}


DTTFId::DTTFId(const DTTFId& id) 
{

  _wh = id._wh;
  _sec = id._sec;
  
}

std::string DTTFId::name() const
{ 
  
  std::stringstream name;
  
  name << "Wh" << _wh << "Sc" << _sec;
  
  return name.str();
  
}


// --------------------------------------------------
// Id class to handle chamber pair and the origin
// objects used to build them (e.g. inner RPC layer)
// --------------------------------------------------

class ChambPairId {

public:

  enum chamb_objects {DTIN=0, DTCORR, DTDIR, DTOUT, NONE};

  ChambPairId(int wh, int sec, int inCh, int outCh, int inChObj, int outChObj);
  ChambPairId(const ChambPairId & id);

  ~ChambPairId() { };
  
  int rawId() const;
  std::string name() const;

  const DTTFId & dttfId() const { return _dttfId; };

  std::string inObjName()  const { return MBPtChambObjectName[_inChObj]; }; 
  std::string outObjName() const { return MBPtChambObjectName[_outChObj]; };
  
private:
  
  DTTFId _dttfId;
  int _inCh, _outCh;
  int _inChObj, _outChObj;

  std::string MBPtChambObjectName[5];
  
};


ChambPairId::ChambPairId(int wh, int sec, int inCh, int outCh, int inChObj, int outChObj) :
  _dttfId(wh,sec), _inCh(inCh), _outCh(outCh), _inChObj(inChObj), _outChObj(outChObj) 
{

  MBPtChambObjectName[DTIN]   = "DTIN";
  MBPtChambObjectName[DTCORR] = "DTCORR";
  MBPtChambObjectName[DTDIR]  = "DTDIR";
  MBPtChambObjectName[DTOUT]  = "DTOUT";
  MBPtChambObjectName[NONE]   = "NONE";
    
}


ChambPairId::ChambPairId(const ChambPairId & id) :
  _dttfId(id._dttfId)
{
  
  _inCh = id._inCh;
  _outCh = id._outCh;
  _inChObj = id._inChObj;
  _outChObj = id._outChObj;

  MBPtChambObjectName[DTIN]   =   id.MBPtChambObjectName[DTIN];
  MBPtChambObjectName[DTCORR] =   id.MBPtChambObjectName[DTCORR];
  MBPtChambObjectName[DTDIR]  =   id.MBPtChambObjectName[DTDIR];
  MBPtChambObjectName[DTOUT]  =   id.MBPtChambObjectName[DTOUT];
  MBPtChambObjectName[NONE]   =   id.MBPtChambObjectName[NONE];
  
}


int ChambPairId::rawId() const
{ 

  int id = _dttfId.rawId() + 1000*_inCh + 10000*_outCh 
           + 100000*_inChObj + 10000000*_outChObj;
  
  return id;
  
}


std::string ChambPairId::name() const
{ 
  
  std::stringstream name;
  
  name << _dttfId.name() 
       << "inCh" << _inCh << "outCh" <<_outCh 
       << inObjName() << outObjName();
  
  return name.str();
  
}


// --------------------------------------------------
// Class holding plot (book, fill, draw and save)
// methods for a given chamb object pair
// --------------------------------------------------

class ChambPairPlotter {
  
public:

  ChambPairPlotter(int wh, int sec, int inCh, int outCh, 
		   int inChObj, int outChObj, TFileService * fs);

  ChambPairPlotter(ChambPairId id, TFileService * fs);

  ChambPairPlotter(const ChambPairPlotter & plotter);

  int rawId() const { return _id.rawId(); };
  std::string name() const { return _id.name(); };

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void book(TFileService * fs);
  void fill(const L1ITMu::TriggerPrimitive * in, const L1ITMu::TriggerPrimitive * out, float pt);
  void draw() const;

private:

  ChambPairId _id;

  TProfile * _pdPhivsPt;
  TH2F     * _hdPhivsPt;

  TProfile * _pdPhiBendvsPt;
  TH2F     * _hdPhiBendvsPt;

  TProfile * _pdPhivsPtCenter;
  TH2F     * _hdPhivsPtCenter;

  TProfile * _pdPhiBendvsPtCenter;
  TH2F     * _hdPhiBendvsPtCenter;

  TProfile * _pdPhivsPtBorder;
  TH2F     * _hdPhivsPtBorder;

  TProfile * _pdPhiBendvsPtBorder;
  TH2F     * _hdPhiBendvsPtBorder;
};


ChambPairPlotter::ChambPairPlotter(int wh, int sec, int inCh, int outCh, 
				   int inChObj, int outChObj, TFileService * fs) :
  _id(wh,sec,inCh,outCh,inChObj,outChObj) 
{ 
  
  book(fs); 
  
}


ChambPairPlotter::ChambPairPlotter(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  book(fs); 
  
}


ChambPairPlotter::ChambPairPlotter(const ChambPairPlotter & plotter) :
  _id(plotter._id), _pdPhivsPt(plotter._pdPhivsPt), _hdPhivsPt(plotter._hdPhivsPt), _pdPhiBendvsPt(plotter._pdPhiBendvsPt), _hdPhiBendvsPt(plotter._hdPhiBendvsPt), _pdPhivsPtCenter(plotter._pdPhivsPtCenter), _hdPhivsPtCenter(plotter._hdPhivsPtCenter), _pdPhiBendvsPtCenter(plotter._pdPhiBendvsPtCenter), _hdPhiBendvsPtCenter(plotter._hdPhiBendvsPtCenter), _pdPhivsPtBorder(plotter._pdPhivsPtBorder), _hdPhivsPtBorder(plotter._hdPhivsPtBorder), _pdPhiBendvsPtBorder(plotter._pdPhiBendvsPtBorder), _hdPhiBendvsPtBorder(plotter._hdPhiBendvsPtBorder)
{
  
}

void ChambPairPlotter::fill(const L1ITMu::TriggerPrimitive * in, const L1ITMu::TriggerPrimitive * out, float pt) 
{ 
  
  float inPhiValue  = in->getCMSGlobalPhi();  
  float outPhiValue = out->getCMSGlobalPhi();
  float inPhiBendValue  = in->getDTData().bendingAngle;
  float outPhiBendValue = out->getDTData().bendingAngle;
  
  float inRelPhiValue = fmod(inPhiValue,(M_PI/6)); // phi relative to the center of the inner chamber (-pi/12;pi/12)
  if (inRelPhiValue > (M_PI/12) ) inRelPhiValue -= (M_PI/6);

  float outRelPhiValue = fmod(outPhiValue,(M_PI/6));  // phi relative to the center of the outner chamber (-pi/12;pi/12)
  if (outRelPhiValue> (M_PI/12) ) outRelPhiValue -= (M_PI/6);  
  
  float deltaPhi     = fabs(inPhiValue-outPhiValue); // CB for now it is fabs, change to mu+ mu- 
  float deltaPhiBend = fabs(inPhiBendValue-outPhiBendValue); // CB for now it is fabs, change to mu+ mu- 

  _pdPhivsPt->Fill(pt,deltaPhi);
  _hdPhivsPt->Fill(pt,deltaPhi);

  _pdPhiBendvsPt->Fill(pt,deltaPhiBend);
  _hdPhiBendvsPt->Fill(pt,deltaPhiBend);  
  
  if ( ( fabs(inRelPhiValue) < M_PI/24 ) && ( fabs(outRelPhiValue) < M_PI/24 ) ){ // fill ch-center histos
    _pdPhivsPtCenter->Fill(pt,deltaPhi);
    _hdPhivsPtCenter->Fill(pt,deltaPhi);
    
    _pdPhiBendvsPtCenter->Fill(pt,deltaPhiBend);
    _hdPhiBendvsPtCenter->Fill(pt,deltaPhiBend);  
  } 
  else { // fill ch-border histos
    _pdPhivsPtBorder->Fill(pt,deltaPhi);
    _hdPhivsPtBorder->Fill(pt,deltaPhi);
    
    _pdPhiBendvsPtBorder->Fill(pt,deltaPhiBend);
    _hdPhiBendvsPtBorder->Fill(pt,deltaPhiBend);  
  }
   
}

void ChambPairPlotter::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
  
  TFileDirectory folder  = fs->mkdir(hDir.c_str());
  
  _pdPhivsPt = folder.make<TProfile>(("pdPhivsPt"+hName).c_str(), 
                                     ("obj #Delta#phi vs pt for "+hName).c_str(), 
                                     60, -0.5, 119.5, 0, .05);

  _hdPhivsPt = folder.make<TH2F>(("hdPhivsPt"+hName).c_str(), 
                                 ("obj #Delta#phi vs pt for "+hName).c_str(),
                                 60, -0.5, 119.5, 100, -0, .05);

  _pdPhiBendvsPt = folder.make<TProfile>(("pdPhiBendvsPt"+hName).c_str(), 
                                     ("obj #Delta#phi_{b} vs pt for "+hName).c_str(), 
                                     60, -0.5, 119.5, 0, 100);

  _hdPhiBendvsPt = folder.make<TH2F>(("hdPhiBendvsPt"+hName).c_str(), 
                                 ("obj #Delta#phi_{b} vs pt for "+hName).c_str(),
                                 60, -0.5, 119.5, 100, 0, 100);

  _pdPhivsPtCenter = folder.make<TProfile>(("pdPhivsPtCenter"+hName).c_str(), 
                                     ("obj #Delta#phi vs pt for "+hName+" - center").c_str(), 
                                     60, -0.5, 119.5, 0, .05);

  _hdPhivsPtCenter = folder.make<TH2F>(("hdPhivsPtCenter"+hName).c_str(), 
                                 ("obj #Delta#phi vs pt for "+hName+" - center").c_str(),
                                 60, -0.5, 119.5, 100, -0, .05);

  _pdPhiBendvsPtCenter = folder.make<TProfile>(("pdPhiBendvsPtCenter"+hName).c_str(), 
                                     ("obj #Delta#phi_{b} vs pt for "+hName+" - center").c_str(), 
                                     60, -0.5, 119.5, 0, 100);

  _hdPhiBendvsPtCenter = folder.make<TH2F>(("hdPhiBendvsPtCenter"+hName).c_str(), 
                                 ("obj #Delta#phi_{b} vs pt for "+hName+" - center").c_str(),
                                 60, -0.5, 119.5, 100, 0, 100);
  
  _pdPhivsPtBorder = folder.make<TProfile>(("pdPhivsPtBorder"+hName).c_str(), 
                                     ("obj #Delta#phi vs pt for "+hName+" - border").c_str(), 
                                     60, -0.5, 119.5, 0, .05);

  _hdPhivsPtBorder = folder.make<TH2F>(("hdPhivsPtBorder"+hName).c_str(), 
                                 ("obj #Delta#phi vs pt for "+hName+" - border").c_str(),
                                 60, -0.5, 119.5, 100, -0, .05);

  _pdPhiBendvsPtBorder = folder.make<TProfile>(("pdPhiBendvsPtBorder"+hName).c_str(), 
                                     ("obj #Delta#phi_{b} vs pt for "+hName+" - border").c_str(), 
                                     60, -0.5, 119.5, 0, 100);

  _hdPhiBendvsPtBorder = folder.make<TH2F>(("hdPhiBendvsPtBorder"+hName).c_str(), 
                                 ("obj #Delta#phi_{b} vs pt for "+hName+" - border").c_str(),
                                 60, -0.5, 119.5, 100, 0, 100);
  
  std::string title = _id.inObjName() + " - " + _id.outObjName();
  
  _pdPhivsPt->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhivsPt->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhivsPt->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhivsPt->GetYaxis()->SetTitle(title.c_str());

  _pdPhiBendvsPt->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhiBendvsPt->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhiBendvsPt->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhiBendvsPt->GetYaxis()->SetTitle(title.c_str());

  _pdPhivsPtCenter->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhivsPtCenter->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhivsPtCenter->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhivsPtCenter->GetYaxis()->SetTitle(title.c_str());

  _pdPhiBendvsPtCenter->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhiBendvsPtCenter->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhiBendvsPtCenter->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhiBendvsPtCenter->GetYaxis()->SetTitle(title.c_str());

  _pdPhivsPtBorder->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhivsPtBorder->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhivsPtBorder->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhivsPtBorder->GetYaxis()->SetTitle(title.c_str());

  _pdPhiBendvsPtBorder->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _pdPhiBendvsPtBorder->GetYaxis()->SetTitle(title.c_str());
  
  _hdPhiBendvsPtBorder->GetXaxis()->SetTitle("GEN  mu p_{T}");
  _hdPhiBendvsPtBorder->GetYaxis()->SetTitle(title.c_str());
}


void ChambPairPlotter::draw() const 
{ 

  std::string cName;
  
  cName = "cPhi" + name();

  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphi = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphi->cd();
  cphi->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhivsPt->Draw("colz");
  _pdPhivsPt->Draw("sameP");
  
  cphi->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  ///
  
  cName = "cPhiBend" + name();
  
  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphibend = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphibend->cd();
  cphibend->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhiBendvsPt->Draw("colz");
  _pdPhiBendvsPt->Draw("sameP");
  
  cphibend->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  ///

  cName = "cPhiCenter" + name();

  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphicent = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphicent->cd();
  cphicent->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhivsPtCenter->Draw("colz");
  _pdPhivsPtCenter->Draw("sameP");
  
  cphicent->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  ///
  
  cName = "cPhiBendCenter" + name();
  
  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphibendcent = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphibendcent->cd();
  cphibendcent->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhiBendvsPtCenter->Draw("colz");
  _pdPhiBendvsPtCenter->Draw("sameP");
  
  cphibendcent->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  ///

  cName = "cPhiBorder" + name();

  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphibord = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphibord->cd();
  cphibord->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhivsPtBorder->Draw("colz");
  _pdPhivsPtBorder->Draw("sameP");
  
  cphibord->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  ///
  
  cName = "cPhiBendBorder" + name();
  
  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());

  TCanvas * cphibendbord = new TCanvas(cName.c_str(),cName.c_str(),500,500);
  cphibendbord->cd();
  cphibendbord->SetGrid();  

  gStyle->SetPalette(52);
  _hdPhiBendvsPtBorder->Draw("colz");
  _pdPhiBendvsPtBorder->Draw("sameP");
  
  cphibendbord->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

}


// --------------------------------------------------
// Analyzer class, holds all ChambPairPlotters 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtLutPlots : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtLutPlots( const edm::ParameterSet & );
  ~L1ITMBPtLutPlots();

  void analyze( const edm::Event &, const edm::EventSetup & );  

  void endJob();

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  ChambPairPlotter* getPlotter(int dttfRawId, int chPairRawId); 

private:



  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::map< int, std::map<int, ChambPairPlotter*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

};


L1ITMBPtLutPlots::L1ITMBPtLutPlots(const edm::ParameterSet& p) : 
  _mbTracksTag( p.getParameter<edm::InputTag>("MBTracksCollection") ),
  _genParticlesTag( p.getParameter<edm::InputTag>("GenParticlesCollection") )

{

  // CB: combinations to be plotted go here 
  // for now just booking DT IN, OUT CORR
  // and ch 1 and 2 of DTTF wh 1 sec 1
  
  std::vector<int> mbPtChambObjects;   

  mbPtChambObjects.push_back(ChambPairId::DTIN);
  mbPtChambObjects.push_back(ChambPairId::DTOUT);
  mbPtChambObjects.push_back(ChambPairId::DTCORR);

  for (int wheel = -1; wheel <=-1; ++wheel) {
    for (int sector = 0; sector <=0; ++sector) {
      for (int inCh = 1; inCh <=1; ++inCh) {
	for (int outCh = inCh + 1; outCh <=2; ++outCh) {
	  
	  std::vector<int>::const_iterator inObjIt = mbPtChambObjects.begin();
	  std::vector<int>::const_iterator objEnd  = mbPtChambObjects.end();
	  
	  for (; inObjIt != objEnd; ++inObjIt) {
	    std::vector<int>::const_iterator outObjIt = mbPtChambObjects.begin();
	    for (; outObjIt != objEnd; ++outObjIt) {
	      
	      ChambPairId id(wheel, sector, inCh, outCh, (*inObjIt), (*outObjIt));
	      ChambPairPlotter *plotter = new ChambPairPlotter(id,_fs.operator->());  

	      histos[DTTFId(wheel,sector).rawId()][id.rawId()] = plotter;
	      
	    }
	  }
	}
      }
    }
  }
  
}
  
L1ITMBPtLutPlots::~L1ITMBPtLutPlots()
{

}




void L1ITMBPtLutPlots::endJob()
{

  std::map<int,std::map<int,ChambPairPlotter*> >::iterator dttfIt  = histos.begin();
  std::map<int,std::map<int,ChambPairPlotter*> >::iterator dttfEnd = histos.end();
  
  for(; dttfIt!=dttfEnd; ++dttfIt) {
    
    std::map<int,ChambPairPlotter*>::iterator chambPairIt  = dttfIt->second.begin();
    std::map<int,ChambPairPlotter*>::iterator chambPairEnd = dttfIt->second.end();
    
    for(; chambPairIt!=chambPairEnd; ++chambPairIt) {
      
      ChambPairPlotter *plotter = chambPairIt->second;

      plotter->draw();
      delete plotter;
      
    }
  }
  
}

void L1ITMBPtLutPlots::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtLutPlots]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
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

    int sector = dttf.scNum();
    int wheel  = dttf.whNum();

    if (wheel!=-1 || sector!=0) continue; // CB hack for test just using a given sector
    
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
	std::cout << "[L1ITMBPtLutPlots]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
	continue;
      }

      dtBestPrims[iSt] = dtMatch.get();

    }

    if (dtBestPrims[0] && dtBestPrims[1]) { // has MB1 and MB2
      
      int mb1Obj = objFromPrim(dtBestPrims[0]);
      int mb2Obj = objFromPrim(dtBestPrims[1]);
      
      int chPairRawId = ChambPairId(wheel,sector,1,2,mb1Obj,mb2Obj).rawId();
      
      ChambPairPlotter *plotter = getPlotter(dttfRawId,chPairRawId);

      if (plotter) plotter->fill(dtBestPrims[0], 
                                 dtBestPrims[1],
                                 bestGen->pt());
    }
    
  }
  
}


int L1ITMBPtLutPlots::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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


ChambPairPlotter* L1ITMBPtLutPlots::getPlotter(int dttfRawId, int chPairRawId) 
{
  // CB optimise 
  if (histos.find(dttfRawId) == histos.end()) { 
    std::cout << "[L1ITMBPtLutPlots]::getPlotter : DTTF not in histos. Return 0.\n";
    return 0;
  }
  if (histos[dttfRawId].find(chPairRawId) == histos[dttfRawId].end()) {
    std::cout << "[L1ITMBPtLutPlots]::getPlotter : ChambPair not in histos. Return 0.\n";
    return 0;
  }
  
  return histos[dttfRawId][chPairRawId];

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtLutPlots);
