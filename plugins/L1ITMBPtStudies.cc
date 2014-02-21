// 
// Class: L1ITMBPtStudies 
//
// Info: Performs GEN-DTTF matching and computes LUT plots  
//
// Author: 
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


class ChPairPlotter {
  
public:

  ChPairPlotter(int wh, int sec, int inCh, int outCh, 
		   int inChObj, int outChObj, TFileService * fs);

  ChPairPlotter(ChambPairId id, TFileService * fs);

  ChPairPlotter(const ChPairPlotter & plotter);

  int rawId() const { return _id.rawId(); };
  std::string name() const { return _id.name(); };

  const TH1      * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void book(TFileService * fs);

  void fillPtLut(const L1ITMu::TriggerPrimitive * in, 
		 const L1ITMu::TriggerPrimitive * out, float pt);

  void fillResolution(float dttfPt, float phiBInPt, float phiBOutpt, float DeltaphiPt, float pt,  Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);

  void SetWeight( Float_t & whGMT, Float_t & whBIn, Float_t &whBOut, Float_t & whDPhi, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);
 
  void draw() const;

private:

  ChambPairId _id;

  std::vector<int> _dttfPtCuts;
  std::map<std::string, TH1 *>    _hPlots;


};


ChPairPlotter::ChPairPlotter(int wh, int sec, int inCh, int outCh, 
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


ChPairPlotter::ChPairPlotter(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  // CB DTTF pt cuts for plots (now one per object)  
  _dttfPtCuts.push_back(12);
  _dttfPtCuts.push_back(16);
  _dttfPtCuts.push_back(20);
  _dttfPtCuts.push_back(24);

  book(fs); 
  
}


ChPairPlotter::ChPairPlotter(const ChPairPlotter & plotter) :
  _id(plotter._id), _dttfPtCuts(plotter._dttfPtCuts),
 _hPlots(plotter._hPlots)
{

}

void ChPairPlotter::fillPtLut(const L1ITMu::TriggerPrimitive * in, 
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

void ChPairPlotter::SetWeight( Float_t & whGMT, Float_t & whBIn, Float_t &whBOut, Float_t & whDPhi, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName){
    
  std::string  FileName="../test/Macros/Weight/"+PlaceName+".txt";
  
  //std::cout<<"FileName "<<FileName.c_str()<<std::endl;
  std::ifstream datafileIn(FileName.c_str(),std::ifstream::in);
  
  std::string line;
  
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
   
  while (getline(datafileIn,line)){
    std::stringstream readline; 
    std::string D1, D2;
    readline<<line;
    
    readline>>D1>>D2>>whGMT>>whDPhi>>whBIn>>whBOut;
    //std::cout<<"line "<<line.c_str()<<std::endl<<"D1 "<<D1.c_str()<<" D2 "<<D2.c_str()<<std::endl;
    if(D1 != QualInString && D2 != QualOutString ) continue;
   
    // std::cout<<"D1 "<<D1.c_str()<<" D2 " <<D2.c_str()<<" GMT "<<whGMT<<std::endl;
    break; 
  }
  whGMT=1./(whGMT*whGMT);
  whDPhi=1./(whDPhi*whDPhi);
  whBIn=1./(whBIn*whBIn);
  whBOut=1./(whBOut*whBOut);
}



void ChPairPlotter::fillResolution(float dttfPt, float phiBInPt, float phiBOutPt, float DeltaphiPt, float pt, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName) 
{ 

  _hPlots["hGMTPtResolvsPt"] ->Fill(pt,(pt - dttfPt)/ pt);

  if(phiBInPt != -1)  _hPlots["phiBInPtResolvsPt"] ->Fill(pt,(pt - phiBInPt)/ pt);
  if(phiBOutPt != -1)  _hPlots["phiBOutPtResolvsPt"] ->Fill(pt,(pt - phiBOutPt)/ pt);
  if(DeltaphiPt != -1) _hPlots["DeltaphiPtResolvsPt"] ->Fill(pt,(pt - DeltaphiPt)/ pt);
  
  Float_t whGmt = 1;
  Float_t whBIn = 1;
  Float_t whBOut = 1;
  Float_t whDPhi = 1;
  
  ChPairPlotter::SetWeight( whGmt,  whBIn, whBOut, whDPhi,  mb1Obj, mb2Obj, PlaceName);
  //std::cout<<"Dt1 "<<mb1Obj<<" Dt2 "<<mb2Obj<<" "<< whGmt<<" "<< whBIn<<" "<<whBOut<<" "<<whDPhi<<std::endl;
    
    Float_t SumWh = whGmt+whBIn+whBOut+whDPhi;
    Float_t SumWhOut = whGmt+whBOut+whDPhi;
    Float_t SumWhIn = whGmt+whBIn+whDPhi;
    
  for(Int_t Cut = 35; Cut<=80; Cut+=15){
    
    if(pt<=Cut && DeltaphiPt != -1 && phiBInPt != -1 && phiBOutPt != -1){
    std::stringstream ptCut; ptCut << Cut;
    std::string PtTag = ptCut.str();
    PtTag="Cut_"+PtTag;
   
    _hPlots["GMTptResol" + PtTag]->Fill((pt - dttfPt)/ pt);
    
    _hPlots["PhiBInptResol" + PtTag]->Fill((pt - phiBInPt)/ pt);
    
    _hPlots["PhiBOutptResol" + PtTag]->Fill((pt - phiBOutPt)/ pt);
    
    _hPlots["DeltaPhiptResol" + PtTag]->Fill((pt - DeltaphiPt)/ pt);
    
    _hPlots["MeanptResol" + PtTag]->Fill((pt - (phiBOutPt+phiBInPt+DeltaphiPt+dttfPt)/4.)/pt);

    _hPlots["WeightMeanptResol" + PtTag]->Fill((pt - (phiBOutPt*whBOut+phiBInPt*whBIn+DeltaphiPt*whDPhi+dttfPt*whGmt)/SumWh)/pt);

    _hPlots["WeightMeanInptResol" + PtTag]->Fill((pt - (phiBInPt*whBIn+DeltaphiPt*whDPhi+dttfPt*whGmt)/SumWhIn)/pt);

    _hPlots["WeightMeanOutptResol" + PtTag]->Fill((pt - (phiBOutPt*whBOut+DeltaphiPt*whDPhi+dttfPt*whGmt)/SumWhOut)/pt);

    }
  }
}



void ChPairPlotter::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
  
  TFileDirectory baseFolder = fs->mkdir(hDir.c_str());
      

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
						     ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
						      ";GEN  mu p_{T};(GEN  mu p_{T} - DTTF  mu p_{T}) / GEN  mu p_{T}").c_str(), 
						     60, -0.5, 119.5, 120, -2.5, 2.5);
  
  
  _hPlots["phiBInPtResolvsPt"] = folderRes.make<TH2F>(("hPhiBendInPtResolvsPt" + hName ).c_str(), 
						      ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
						       ";GEN  mu p_{T};(GEN  mu p_{T} - #phi_{b}In  mu p_{T}) / GEN  mu p_{T}").c_str(), 
						      60, -0.5, 119.5, 120, -2.5, 2.5);
  
  _hPlots["phiBOutPtResolvsPt"] = folderRes.make<TH2F>(("hPhiBendOutPtResolvsPt" + hName ).c_str(), 
						       ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							";GEN  mu p_{T};(GEN  mu p_{T} - #phi_{b}Out  mu p_{T}) / GEN  mu p_{T}").c_str(), 
						       60, -0.5, 119.5, 120, -2.5, 2.5);
  
  _hPlots["DeltaphiPtResolvsPt"] = folderRes.make<TH2F>(("hDeltaPhiPtResolvsPt" + hName ).c_str(), 
							("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							 ";GEN  mu p_{T};(GEN  mu p_{T} - #Delta #phi  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							60, -0.5, 119.5, 120, -2.5, 2.5);



  for(Int_t Cut = 35; Cut<=80; Cut+=15){


      std::stringstream ptCut; ptCut << Cut;
      std::string PtTag = ptCut.str();
      PtTag="Cut_"+PtTag;


      _hPlots["GMTptResol" + PtTag] = folderRes.make<TH1F>(("hGMTPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - DTTF  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);


       _hPlots["DeltaPhiptResol" + PtTag] = folderRes.make<TH1F>(("hDeltaPhiPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - #Delta #phi  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);


       _hPlots["PhiBInptResol" + PtTag] = folderRes.make<TH1F>(("hPhiBInPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - #phi_{b}In  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);

       _hPlots["PhiBOutptResol" + PtTag] = folderRes.make<TH1F>(("hPhiBOutPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - #phi_{b}Out  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);


       _hPlots["MeanptResol" + PtTag] = folderRes.make<TH1F>(("hMeanPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - MeanLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);


       _hPlots["WeightMeanptResol" + PtTag] = folderRes.make<TH1F>(("hWMeanPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - WeightMeanLut  mu p_{T}) / GEN  mu p_{T}").c_str(),  //FIXME
							    120, -5., 1.);

       _hPlots["WeightMeanInptResol" + PtTag] = folderRes.make<TH1F>(("hWMeanInPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - WeightMeanInLut  mu p_{T}) / GEN  mu p_{T}").c_str(),  //FIXME
							    120, -5., 1.);

       
       _hPlots["WeightMeanOutptResol" + PtTag] = folderRes.make<TH1F>(("hWMeanOutPtResol" + hName + PtTag).c_str(), 
								      ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
								       ";(GEN  mu p_{T} - WeightMeanOutInLut  mu p_{T}) / GEN  mu p_{T}").c_str(),  //FIXME
								      120, -5., 1.);
  }
}


void ChPairPlotter::draw() const 
{ 

  gStyle->SetPalette(52);

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
}


// --------------------------------------------------
// Analyzer class, holds all ChPairPlotters 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtStudies : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtStudies( const edm::ParameterSet & );
  ~L1ITMBPtStudies();

  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  float ptLut(std::string LutName,float PhiValue, Int_t QualIn, Int_t QualOut);


  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  ChPairPlotter* getPlotter(int dttfRawId, int chPairRawId); 


private:



  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChPairPlotter*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

};


L1ITMBPtStudies::L1ITMBPtStudies(const edm::ParameterSet& p) : 
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
	      ChPairPlotter *plotter = new ChPairPlotter(id,_fs.operator->());  

	      histos[DTTFId(wheel,sector).rawId()][id.rawId()] = plotter;
	      
	    }
	  }
	}
      }
    }
  }
  
}
 


float L1ITMBPtStudies::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}




float L1ITMBPtStudies::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);
  
}




float L1ITMBPtStudies::ptLut(std::string LutName,float PhiValue, Int_t QualIn, Int_t QualOut)
{

  bool PhiFound = 0;
  LutName="../test/my_luts/"+LutName;

  std::ifstream datafileIn(LutName.c_str(),std::ifstream::in);

  Float_t PtLut = -1;
  Float_t Phi1 =1;
  Float_t Phi2 =0;
  std::string line;

  std::string QualInString, QualOutString;

  switch(QualIn){

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

 switch(QualOut){

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

while (getline(datafileIn,line)){
		  
		
		  std::stringstream readline; 
		  std::string A, D1, D2;
		  readline<<line;
		  readline>>A>>A>>D1>>D2>>PtLut>>Phi2>>A;
    
		  if(D1 != QualInString && D2 != QualOutString ) continue;
    
		  if(PhiValue>=Phi2&&PhiValue<Phi1){
		    PhiFound=1;
		    // std::cout<<" Found   pt = "<<ptLut<<std::endl<<std::endl;
		    break;
		  }
		  
		  Phi1=Phi2; 
 }

 if(PhiFound)return PtLut;
 else return -1.;
}




 
L1ITMBPtStudies::~L1ITMBPtStudies()
{

}


void L1ITMBPtStudies::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtStudies]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
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
	std::cout << "[L1ITMBPtStudies]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
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
      
      ChPairPlotter *plotter = getPlotter(dttfRawId,chPairRawId);

      if (plotter) 
	{
	  plotter->fillPtLut(dtBestPrims[0], 
			     dtBestPrims[1],
			     bestGen->pt());


 if (gmtInPt >= 10.) 
   {	
     
     Float_t DeltaPhi =  GetDeltaPhi(dtBestPrims[0],dtBestPrims[1]);  
     Float_t PhiBendIn =  GetPhiBend(dtBestPrims[0]); 
     Float_t PhiBendOut =  GetPhiBend(dtBestPrims[1]);   
     
     Float_t ptBendIn = -1;
     Float_t ptBendOut = -1;
     Float_t ptDeltaPhi = -1;
     
     std::string BendInName = "PhiBendIn"+PlaceName.str();
     ptBendIn =  ptLut(BendInName.c_str(), PhiBendIn, mb1Obj,mb2Obj);    
     
     std::string BendOutName = "PhiBendOut"+PlaceName.str();
     ptBendOut =  ptLut(BendOutName.c_str(), PhiBendOut, mb1Obj,mb2Obj);    
     
     std::string DeltaPhiName = "Phi"+PlaceName.str();
     ptDeltaPhi =  ptLut(DeltaPhiName.c_str(), DeltaPhi, mb1Obj,mb2Obj);    		
     
     plotter->fillResolution(gmtInPt, ptBendIn, ptBendOut, ptDeltaPhi, bestGen->pt(),mb1Obj,mb2Obj,PlaceName.str()); 	          
         
   }
 
	}	  
    }
    
  }
  
}


int L1ITMBPtStudies::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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


ChPairPlotter* L1ITMBPtStudies::getPlotter(int dttfRawId, int chPairRawId) 
{
  // CB optimise 
  if (histos.find(dttfRawId) == histos.end()) { 
    std::cout << "[L1ITMBPtStudies]::getPlotter : DTTF not in histos. Return 0.\n";
    return 0;
  }
  if (histos[dttfRawId].find(chPairRawId) == histos[dttfRawId].end()) {
    std::cout << "[L1ITMBPtStudies]::getPlotter : ChambPair not in histos. Return 0.\n";
    return 0;
  }
  return histos[dttfRawId][chPairRawId];
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtStudies);
