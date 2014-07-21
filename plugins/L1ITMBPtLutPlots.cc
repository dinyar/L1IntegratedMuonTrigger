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

#include "L1Trigger/L1IntegratedMuonTrigger/interface/ChambPairId.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace L1ITMu;

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

  const TProfile * getProfile(std::string name) const { return _pPlots.find(name)->second; }
  const TH1      * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void book(TFileService * fs);

  void fillPtLut(const L1ITMu::TriggerPrimitive * in, 
		 const L1ITMu::TriggerPrimitive * out, float pt);

  void fillEfficicency(float dttfPt, float pt);
  void draw() const;

private:

  ChambPairId _id;

  std::vector<int> _dttfPtCuts;
  
  std::map<std::string, TEfficiency *> _ePlots;
  std::map<std::string, TProfile *>    _pPlots;
  std::map<std::string, TH1      *>    _hPlots;


};


ChambPairPlotter::ChambPairPlotter(int wh, int sec, int inCh, int outCh, 
				   int inChObj, int outChObj, TFileService * fs) :
  _id(wh,sec,inCh,outCh,inChObj,outChObj) 
{ 
  
  // CB DTTF pt cuts for plots (now one per object)  
  _dttfPtCuts.push_back(12);
  _dttfPtCuts.push_back(16);
  _dttfPtCuts.push_back(20);
  _dttfPtCuts.push_back(24);

  book(fs); 
  
}


ChambPairPlotter::ChambPairPlotter(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  // CB DTTF pt cuts for plots (now one per object)  
  _dttfPtCuts.push_back(12);
  _dttfPtCuts.push_back(16);
  _dttfPtCuts.push_back(20);
  _dttfPtCuts.push_back(24);

  book(fs); 
  
}


ChambPairPlotter::ChambPairPlotter(const ChambPairPlotter & plotter) :
  _id(plotter._id), _dttfPtCuts(plotter._dttfPtCuts), _ePlots(plotter._ePlots), 
  _pPlots(plotter._pPlots), _hPlots(plotter._hPlots)
{

}

void ChambPairPlotter::fillPtLut(const L1ITMu::TriggerPrimitive * in, 
				 const L1ITMu::TriggerPrimitive * out, float pt) 
{ 
  
  float inPhiValue  = in->getCMSGlobalPhi();  
  float outPhiValue = out->getCMSGlobalPhi();
  float inPhiBendValue  = fabs(in->getDTData().bendingAngle);
  float outPhiBendValue = fabs(out->getDTData().bendingAngle);
  
  int inPhiBin  = acos(cos(inPhiValue))  / (M_PI/24.); // phi bin (0,48)
  int outPhiBin = acos(cos(outPhiValue)) / (M_PI/24.); // phi bin (0;48)

  int inRelPhiBin  = (inPhiBin+1)%4;  // phi bin(0,3);
  int outRelPhiBin = (outPhiBin+1)%4; // phi bin(0,3);

  float deltaPhi     = fabs(inPhiValue-outPhiValue); // CB for now it is fabs, change to mu+ mu- 

  _pPlots["dPhivsPt"]->Fill(pt,deltaPhi);
  _hPlots["dPhivsPt"]->Fill(pt,deltaPhi);

  _pPlots["phiBendInvsPt"]->Fill(pt,inPhiBendValue);
  _hPlots["phiBendInvsPt"]->Fill(pt,inPhiBendValue);  

  _pPlots["phiBendOutvsPt"]->Fill(pt,outPhiBendValue);
  _hPlots["phiBendOutvsPt"]->Fill(pt,outPhiBendValue);  
  
  if ( inRelPhiBin<2 && outRelPhiBin<2 )
    { // fill ch-center histos
      _pPlots["dPhivsPtCenter"]->Fill(pt,deltaPhi);
      _hPlots["dPhivsPtCenter"]->Fill(pt,deltaPhi);
    
      _pPlots["phiBendInvsPtCenter"]->Fill(pt,inPhiBendValue);
      _hPlots["phiBendInvsPtCenter"]->Fill(pt,inPhiBendValue);  

      _pPlots["phiBendOutvsPtCenter"]->Fill(pt,outPhiBendValue);
      _hPlots["phiBendOutvsPtCenter"]->Fill(pt,outPhiBendValue);    
    } 
  else if ( inRelPhiBin>=2 && outRelPhiBin>=2 ) 
    { // fill ch-border histos
      _pPlots["dPhivsPtBorder"]->Fill(pt,deltaPhi);
      _hPlots["dPhivsPtBorder"]->Fill(pt,deltaPhi);
    
      _pPlots["phiBendInvsPtBorder"]->Fill(pt,inPhiBendValue);
      _hPlots["phiBendInvsPtBorder"]->Fill(pt,inPhiBendValue);  

      _pPlots["phiBendOutvsPtBorder"]->Fill(pt,outPhiBendValue);
      _hPlots["phiBendOutvsPtBorder"]->Fill(pt,outPhiBendValue);  
    }
   
}


void ChambPairPlotter::fillEfficicency(float dttfPt, float pt) 
{ 

  std::vector<int>::const_iterator dttfPtCutIt  =  _dttfPtCuts.begin();
  std::vector<int>::const_iterator dttfPtCutEnd =  _dttfPtCuts.end();
  
  for (; dttfPtCutIt != dttfPtCutEnd; ++dttfPtCutIt)
    {
      int dttfPtCut = (*dttfPtCutIt);

      std::stringstream ptCut; ptCut << dttfPtCut;
      std::string dttfPtTag = ptCut.str();

      _ePlots["effvsPt" + dttfPtTag]->Fill(dttfPt >= dttfPtCut,pt);
      if (pt > dttfPtCut) _hPlots["ptResol" + dttfPtTag]->Fill((pt - dttfPt)/ pt);
    }

}


void ChambPairPlotter::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
  
  TFileDirectory baseFolder = fs->mkdir(hDir.c_str());
      
  std::string phiReg[3] = { "", "Center", "Border" };

  for (int iPhiReg = 0; iPhiReg < 3; ++ iPhiReg )
    {

      std::string reg = phiReg[iPhiReg];

      TFileDirectory folder = baseFolder.mkdir(reg.c_str());
      
      _pPlots["dPhivsPt" + reg] = folder.make<TProfile>(("pdPhivsPt" + reg + hName).c_str(), 
							("obj #Delta#phi vs pt for " + hName +
							 ";GEN  mu p_{T};#Delta#phi").c_str(), 
							60, 0.,120, 0, .05);

      _hPlots["dPhivsPt" + reg] = folder.make<TH2F>(("hdPhivsPt" + reg + hName).c_str(), 
						    ("obj #Delta#phi vs pt for "+ hName +
						     ";GEN  mu p_{T};#Delta#phi").c_str(),
						    60, 0.,120, 100, -0, .05);

      _pPlots["phiBendInvsPt" + reg] = folder.make<TProfile>(("pdPhiBendInvsPt" + reg + hName).c_str(), 
							     ("obj #Delta#phi_{b}In vs pt for " + hName +
							      ";GEN  mu p_{T};#Delta#phi_{b}In").c_str(), 
							     60, 0.,120, 0, 100);
      
      _hPlots["phiBendInvsPt" + reg] = folder.make<TH2F>(("hdPhiBendInvsPt" + reg + hName).c_str(), 
							 ("obj #Delta#phi_{b}In vs pt for " + hName +
							  ";GEN  mu p_{T};#Delta#phi_{b}In").c_str(),
							 60, 0.,120, 100, 0, 100);


      _pPlots["phiBendOutvsPt" + reg] = folder.make<TProfile>(("pdPhiBendOutvsPt" + reg + hName).c_str(), 
							      ("obj #Delta#phi_{b}Out vs pt for " + hName +
							       ";GEN  mu p_{T};#Delta#phi_{b}Out").c_str(), 
							      60, 0.,120, 0, 100);
      
      _hPlots["phiBendOutvsPt" + reg] = folder.make<TH2F>(("hdPhiBendOutvsPt" + reg + hName).c_str(), 
							  ("obj #Delta#phi_{b}Out vs pt for " + hName +
							   ";GEN  mu p_{T};#Delta#phi_{b}Out").c_str(),
							  60, 0.,120, 100, 0, 100);

    }

  TFileDirectory folderEff = baseFolder.mkdir("Efficiency");

  std::vector<int>::const_iterator dttfPtCutIt  =  _dttfPtCuts.begin();
  std::vector<int>::const_iterator dttfPtCutEnd =  _dttfPtCuts.end();
  
  for (; dttfPtCutIt != dttfPtCutEnd; ++dttfPtCutIt)
    {
      int dttfPtCut = (*dttfPtCutIt);

      std::stringstream ptCut; ptCut << dttfPtCut;
      std::string dttfPtTag = ptCut.str();
  
      _ePlots["effvsPt" + dttfPtTag] = folderEff.make<TEfficiency>(("eEffvsPt" + hName + dttfPtTag).c_str(), 
								   ("obj efficiency vs pt for " + hName + dttfPtTag + 
								    ";GEN  mu p_{T};p_{T} assign. eff.").c_str(), 
								   60, 0.,120);

      _hPlots["ptResol" + dttfPtTag] = folderEff.make<TH1F>(("hPtResol" + hName + dttfPtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + dttfPtTag + 
							     ";(DTTF  mu p_{T} - GEN  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);
    }


}


void ChambPairPlotter::draw() const 
{ 

  gStyle->SetPalette(52);
  
  std::map<std::string, TProfile *>::const_iterator pPlotsIt  =  _pPlots.begin();
  std::map<std::string, TProfile *>::const_iterator pPlotsEnd =  _pPlots.end();

  for (; pPlotsIt != pPlotsEnd; ++pPlotsIt)
    {
 
      TProfile * pHisto = pPlotsIt->second;

      std::string histoName = pPlotsIt->first;
      TH1 * hHisto = _hPlots.find(histoName)->second;
      
      std::string cName = "c" + histoName + name();
      
      system(std::string("mkdir -p plots/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();  
      
      hHisto->Draw("colz");
      pHisto->Draw("sameP");
  
      c->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());
    }

  std::map<std::string, TH1 *>::const_iterator hPlotsIt  =  _hPlots.begin();
  std::map<std::string, TH1 *>::const_iterator hPlotsEnd =  _hPlots.end();

  for (; hPlotsIt != hPlotsEnd ; ++hPlotsIt)
    {

      TH1 * hHisto = hPlotsIt->second;

      if ( hPlotsIt->first.find("ptResol") == std::string::npos )  continue;
 
      std::string histoName = hPlotsIt->first;
      
      std::string cName = "c" + histoName + name();
      
      system(std::string("mkdir -p plots/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();  
      
      hHisto->Draw("");

      c->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());
    }

  std::map<std::string, TEfficiency *>::const_iterator ePlotsIt  =  _ePlots.begin();
  std::map<std::string, TEfficiency *>::const_iterator ePlotsEnd =  _ePlots.end();

  for (; ePlotsIt != ePlotsEnd ; ++ePlotsIt)
    {

      TEfficiency * histo = ePlotsIt->second;
      std::string histoName = ePlotsIt->first;
      
      std::string cName = "c" + histoName + name();
      
      system(std::string("mkdir -p plots/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();  
      
      histo->Draw("E0PA");
      
      c->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());
    }

  TProfile * pCenter = static_cast<TProfile *>(getProfile("dPhivsPtCenter")->Clone("center"));
  TProfile * pBorder = static_cast<TProfile *>(getProfile("dPhivsPtBorder")->Clone("border"));
      
  std::string cName = "cdPhivsPtComparison" + name();
      
  system(std::string("mkdir -p plots/" + dttfId().name()).c_str());
      
  TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
  c->cd();
  c->SetGrid();  
      
  pCenter->SetLineColor(1);
  pBorder->SetLineColor(2);

  pCenter->Draw("P");
  pBorder->Draw("sameP");

  c->SaveAs(("plots/" + dttfId().name() + "/" + cName+".pdf").c_str());

  delete pCenter;
  delete pBorder;

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

  void drawChambPair(DTTFId dttfId, int inCh, int outCh);

private:



  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChambPairPlotter*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

};


L1ITMBPtLutPlots::L1ITMBPtLutPlots(const edm::ParameterSet& p) : 
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
    for (int sector = 0; sector <=0; ++sector) {
      for (int inCh = 1; inCh <=1; ++inCh) {
	for (int outCh = inCh + 1; outCh <=2; ++outCh) {
	  
	  std::vector<int>::const_iterator inObjIt = _mbPtChambObjects.begin();
	  std::vector<int>::const_iterator objEnd  = _mbPtChambObjects.end();
	  
	  for (; inObjIt != objEnd; ++inObjIt) {
	    std::vector<int>::const_iterator outObjIt = _mbPtChambObjects.begin();
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

  for (int wheel = -3; wheel <=3; ++wheel) {
    if (wheel == 0) continue;
    for (int sector = 0; sector <=0; ++sector) {
      for (int inCh = 1; inCh <=1; ++inCh) {
	for (int outCh = inCh + 1; outCh <=2; ++outCh) {
	  drawChambPair(DTTFId(wheel,sector), inCh, outCh); 
	}
      }
    }
  }
  
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

      if (plotter) 
	{
	  plotter->fillPtLut(dtBestPrims[0], 
			     dtBestPrims[1],
			     bestGen->pt());

	  if (gmtInPt >= 0.)
	    {	
	      plotter->fillEfficicency(gmtInPt,bestGen->pt());
	    }

	}
	  
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


void  L1ITMBPtLutPlots::drawChambPair(DTTFId dttfId, int inCh, int outCh) 
{

  std::vector<TProfile *> histos;

  std::vector<int>::const_iterator inObjIt = _mbPtChambObjects.begin();
  std::vector<int>::const_iterator objEnd  = _mbPtChambObjects.end();
	  
  for (; inObjIt != objEnd; ++inObjIt) {
    std::vector<int>::const_iterator outObjIt = _mbPtChambObjects.begin();
    for (; outObjIt != objEnd; ++outObjIt) {
      
      ChambPairId id(dttfId, inCh, outCh, (*inObjIt), (*outObjIt));
      ChambPairPlotter * plotter = getPlotter(dttfId.rawId(),id.rawId());
      TProfile *histo = (TProfile*)(plotter->getProfile("dPhivsPt")->Clone(plotter->name().c_str()));
      histos.push_back(histo);

    }
  }

  // CB fixme should be chambers not DTTF id 
  std::string cName = "c" + dttfId.name() + "dPhivsPtComparison"; 

  system(std::string("mkdir -p plots/" + dttfId.name()).c_str());
      
  TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
  c->cd();
  c->SetGrid();  
      
  std::vector<TProfile *>::iterator histoIt  = histos.begin();
  std::vector<TProfile *>::iterator histoEnd = histos.end();

  for (int iColor = 1; histoIt != histoEnd; ++histoIt, ++iColor)
    {
      (*histoIt)->SetLineColor(iColor);
      (*histoIt)->Draw(iColor == 1 ? "P" :"sameP");
    }
      
  c->SaveAs(("plots/" + dttfId.name() + "/" + cName+".pdf").c_str());

  histoIt  = histos.begin();
  histoEnd = histos.end();

  for (; histoIt != histoEnd; ++histoIt)
    {
      delete (*histoIt);
    }  

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtLutPlots);
