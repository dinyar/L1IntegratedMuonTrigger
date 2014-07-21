// 
// Class: L1ITMBPtComb
//
// Info: Performs Combination of pt measuraments of phi bend in, phi bend out and delta phi.
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


class ChPairPlotterComb{
  
public:

  ChPairPlotterComb(int wh, int sec, int inCh, int outCh, 
		   int inChObj, int outChObj, TFileService * fs);

  ChPairPlotterComb(ChambPairId id, TFileService * fs);

  ChPairPlotterComb(const ChPairPlotterComb & plotter);

  int rawId() const { return _id.rawId(); };
  std::string name() const { return _id.name(); };

  const TH1 * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void draw() const ;

  void book(TFileService * fs);

  void fillResolution(float dttfPt, float phiBInPt, float phiBOutpt, float DeltaphiPt, float pt,  Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);


private:

  ChambPairId _id;

  std::vector<int> _PtResCuts;
  std::map<std::string, TH1 *>  _hPlots;

};


ChPairPlotterComb::ChPairPlotterComb(int wh, int sec, int inCh, int outCh, 
				   int inChObj, int outChObj, TFileService * fs) :
  _id(wh,sec,inCh,outCh,inChObj,outChObj) 
{ 
  
  // GP Gen pt cuts for resolution plots 
  _PtResCuts.push_back(30);
  _PtResCuts.push_back(35);
  _PtResCuts.push_back(40);
  _PtResCuts.push_back(45);
  
  book(fs); 
  
}


ChPairPlotterComb::ChPairPlotterComb(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  // GP Gen pt cuts for resolution plots  
  _PtResCuts.push_back(30);
  _PtResCuts.push_back(35);
  _PtResCuts.push_back(40);
  _PtResCuts.push_back(45);
  
  book(fs); 
  
}


ChPairPlotterComb::ChPairPlotterComb(const ChPairPlotterComb & plotter) :
 _id(plotter._id), _PtResCuts(plotter._PtResCuts),
 _hPlots(plotter._hPlots)
{

}



void ChPairPlotterComb::fillResolution(float dttfPt, float phiBInPt, float phiBOutPt, float DeltaphiPt, float pt, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName) 
{ 
  
  Float_t whWMRMSBIn = 1;
  Float_t whWMRMSBOut = 1;
  Float_t whWMRMSDPhi = 1;
  
  Float_t whWMSigmaBIn = 1;
  Float_t whWMSigmaBOut = 1;
  Float_t whWMSigmaDPhi = 1;

  Float_t whWMInRMSBIn = 1;
  Float_t whWMInRMSBOut = 1;
  Float_t whWMInRMSDPhi = 1;
  
  Float_t whWMInSigmaBIn = 1;
  Float_t whWMInSigmaBOut = 1;
  Float_t whWMInSigmaDPhi = 1;
 
  Float_t whWMOutRMSBIn = 1;
  Float_t whWMOutRMSBOut = 1;
  Float_t whWMOutRMSDPhi = 1;  
  
  Float_t whWMOutSigmaBIn = 1;
  Float_t whWMOutSigmaBOut = 1;
  Float_t whWMOutSigmaDPhi = 1;
  
  
  PtProvider::SetWeight( whWMRMSBIn, whWMRMSBOut, whWMRMSDPhi,  mb1Obj, mb2Obj, PlaceName,"WMean","RMS");
  PtProvider::SetWeight( whWMSigmaBIn, whWMSigmaBOut, whWMSigmaDPhi,  mb1Obj, mb2Obj, PlaceName,"WMean","Sigma");
  
  PtProvider::SetWeight( whWMInRMSBIn, whWMInRMSBOut, whWMInRMSDPhi,  mb1Obj, mb2Obj, PlaceName,"WMeanIn","RMS");
  PtProvider::SetWeight( whWMInSigmaBIn, whWMInSigmaBOut, whWMInSigmaDPhi,  mb1Obj, mb2Obj, PlaceName,"WMeanIn","Sigma");
  
  PtProvider::SetWeight( whWMOutRMSBIn, whWMOutRMSBOut, whWMOutRMSDPhi,  mb1Obj, mb2Obj, PlaceName,"WMeanOut","RMS");
  PtProvider::SetWeight( whWMOutSigmaBIn, whWMOutSigmaBOut, whWMOutSigmaDPhi,  mb1Obj, mb2Obj, PlaceName,"WMeanOut","Sigma"); 
  
  //chose the minimum pt between phi bend pt and delta phi pt. Use only phi bend of correlated chambers.
    
 Float_t MinPt=DeltaphiPt;

  if(mb1Obj==1&&mb2Obj==1){
  
    if(phiBInPt<MinPt) MinPt=phiBInPt;
    if(phiBOutPt<MinPt) MinPt=phiBOutPt;
  }
  
  else if(mb1Obj==1&&mb2Obj!=1){
    if(phiBInPt<MinPt) MinPt=phiBInPt;
  }
  
  else if(mb1Obj!=1&&mb2Obj==1){
    if(phiBOutPt<MinPt) MinPt=phiBOutPt;
  }

  

  Float_t MinPt2;
  
  if(mb1Obj==1 &&  mb2Obj==1){
    if(phiBOutPt<DeltaphiPt && phiBOutPt<phiBInPt && DeltaphiPt<phiBInPt ) MinPt2 = (phiBOutPt+DeltaphiPt)/2.;
    else if (phiBInPt<DeltaphiPt && phiBInPt<phiBOutPt && DeltaphiPt<phiBInPt) MinPt2 = (phiBOutPt+DeltaphiPt)/2.;
    else  MinPt2 = (phiBOutPt+phiBInPt)/2.;
  }
  
  else if(mb1Obj==1 &&  mb2Obj!=1) MinPt2 = (phiBInPt+DeltaphiPt)/2.;
  else if(mb1Obj!=1 &&  mb2Obj==1) MinPt2 = (phiBOutPt+DeltaphiPt)/2.;
  else MinPt2 = DeltaphiPt;
  
  Float_t PtM = (phiBOutPt+phiBInPt+DeltaphiPt)/3.;
  Float_t PtRMSWM = (phiBOutPt*whWMRMSBOut+phiBInPt*whWMRMSBIn+DeltaphiPt*whWMRMSDPhi); 
  Float_t PtRMSWMIn = (phiBInPt*whWMInRMSBIn+DeltaphiPt*whWMInRMSDPhi);
  Float_t PtRMSWMOut = (phiBOutPt*whWMOutRMSBOut+DeltaphiPt*whWMOutRMSDPhi); 
  Float_t PtSigmaWM  = (phiBOutPt*whWMSigmaBOut+phiBInPt*whWMSigmaBIn+DeltaphiPt*whWMSigmaDPhi);
  Float_t PtSigmaWMIn = (phiBInPt*whWMInSigmaBIn+DeltaphiPt*whWMInSigmaDPhi);
  Float_t PtSigmaWMOut = (phiBOutPt*whWMOutSigmaBOut+DeltaphiPt*whWMOutSigmaDPhi);

  
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
	
	_hPlots["MeanptResol" + PtTag]->Fill((pt - PtM)/pt);

       	_hPlots["MinPtResol" + PtTag]->Fill((pt- MinPt)/pt);

	_hPlots["MinPt2Resol" + PtTag]->Fill((pt- MinPt2)/pt);
	
	//RMS weight
	_hPlots["RMSWeightMeanptResol" + PtTag]->Fill((pt - PtRMSWM )/pt);
	
	_hPlots["RMSWeightMeanInptResol" + PtTag]->Fill((pt - PtRMSWMIn)/pt);
	
	_hPlots["RMSWeightMeanOutptResol" + PtTag]->Fill((pt - PtRMSWMOut)/pt);  
	
	//Sigma Weight
	_hPlots["SigmaWeightMeanptResol" + PtTag]->Fill((pt - PtSigmaWM)/pt); 
	
	_hPlots["SigmaWeightMeanInptResol" + PtTag]->Fill((pt - PtSigmaWMIn)/pt); 
	
	_hPlots["SigmaWeightMeanOutptResol" + PtTag]->Fill((pt - PtSigmaWMOut)/pt); 
	
      }
    }
}


void ChPairPlotterComb::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
  
 
   
  TFileDirectory baseFolder = fs->mkdir(hDir.c_str());
           
     TFileDirectory folderRes = baseFolder.mkdir("Resolutions");
     
     _hPlots["hGMTPtResolvsPt" ] = folderRes.make<TH2F>(("hGMTPtResolvsPt" + hName ).c_str(), 
							("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							 ";GEN  mu p_{T};(GEN  mu p_{T} - DTTF  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							60, 0.,120., 120, -2.5, 2.5);
     
     
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
     
     _hPlots["BestPtResolvsPt"] = folderRes.make<TH2F>(("hBestPhiPtResolvsPt" + hName ).c_str(), 
						       ("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							";GEN  mu p_{T};(GEN  mu p_{T} - Best pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
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


       _hPlots["MinPt2Resol" + PtTag] = folderCut.make<TH1F>(("hMinPt2Resol" + hName + PtTag).c_str(), 
							     ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							      ";(GEN  mu p_{T} - MinPt2Lut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							     120, -5., 1.);
       
       
       _hPlots["MinPtResol" + PtTag] = folderCut.make<TH1F>(("hMinPtResol" + hName + PtTag).c_str(), 
							     ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							      ";(GEN  mu p_{T} - MinPtLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							     120, -5., 1.);
       

       _hPlots["MeanptResol" + PtTag] = folderCut.make<TH1F>(("hMeanPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - MeanLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);


       _hPlots["RMSWeightMeanptResol" + PtTag] = folderCut.make<TH1F>(("hRMSWMeanPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - RMSWeightMeanLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);

       _hPlots["RMSWeightMeanInptResol" + PtTag] = folderCut.make<TH1F>(("hRMSWMeanInPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - RMSWeightMeanInLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);

       
       _hPlots["RMSWeightMeanOutptResol" + PtTag] = folderCut.make<TH1F>(("hRMSWMeanOutPtResol" + hName + PtTag).c_str(), 
								      ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
								       ";(GEN  mu p_{T} - RMSWeightMeanOutInLut  mu p_{T}) / GEN  mu p_{T}").c_str(),
								      120, -5., 1.);


       _hPlots["SigmaWeightMeanptResol" + PtTag] = folderCut.make<TH1F>(("hSigmaWMeanPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - SigmaWeightMeanLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);

       _hPlots["SigmaWeightMeanInptResol" + PtTag] = folderCut.make<TH1F>(("hSigmaWMeanInPtResol" + hName + PtTag).c_str(), 
							    ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
							     ";(GEN  mu p_{T} - SigmaWeightMeanInLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							    120, -5., 1.);

       
       _hPlots["SigmaWeightMeanOutptResol" + PtTag] = folderCut.make<TH1F>(("hSigmaWMeanOutPtResol" + hName + PtTag).c_str(), 
								      ("obj (dttf_[pt] - pt) / pt for " + hName + PtTag + 
								       ";(GEN  mu p_{T} - SigmaWeightMeanOutInLut  mu p_{T}) / GEN  mu p_{T}").c_str(), 
							      120, -5., 1.);

  }
}


void ChPairPlotterComb::draw() const
{

  gStyle->SetPalette(52);
  
  
  std::map<std::string, TH1 *>::const_iterator hPlotsIt = _hPlots.begin();
  std::map<std::string, TH1 *>::const_iterator hPlotsEnd = _hPlots.end();

  for (; hPlotsIt != hPlotsEnd ; ++hPlotsIt)
    {

      TH1 * hHisto = hPlotsIt->second;
 
      std::string histoName = hPlotsIt->first;
      
      std::string cName = "c" + histoName + "_2step_" + name();
      
      system(std::string("mkdir -p plots/Resolution/Combinations" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();
      
      hHisto->Draw("");

      c->SaveAs(("plots/Resolution/Combinations" + dttfId().name() + "/" + cName+".pdf").c_str());
    }

}




// --------------------------------------------------
// Analyzer class, holds all ChPairPlotterCombs 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtComb : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtComb( const edm::ParameterSet & );
  ~L1ITMBPtComb();

  void endJob();
  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  ChPairPlotterComb* getPlotter(int dttfRawId, int chPairRawId); 


private:



  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChPairPlotterComb*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

  std::map<int, int> NGen;
};


L1ITMBPtComb::L1ITMBPtComb(const edm::ParameterSet& p) : 
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
	      ChPairPlotterComb *plotter = new ChPairPlotterComb(id,_fs.operator->());  

	      histos[DTTFId(wheel,sector).rawId()][id.rawId()] = plotter;
	      
	    }
	  }
	}
      }
    }
  }  
}
 


float L1ITMBPtComb::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}




float L1ITMBPtComb::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);
  
}
 
L1ITMBPtComb::~L1ITMBPtComb()
{

}

void L1ITMBPtComb::endJob()
{
  
  std::map<int,std::map<int,ChPairPlotterComb*> >::iterator dttfIt = histos.begin();
  std::map<int,std::map<int,ChPairPlotterComb*> >::iterator dttfEnd = histos.end();
  
  for(; dttfIt!=dttfEnd; ++dttfIt) {
    
    std::map<int,ChPairPlotterComb*>::iterator chambPairIt = dttfIt->second.begin();
    std::map<int,ChPairPlotterComb*>::iterator chambPairEnd = dttfIt->second.end();
    
    for(; chambPairIt!=chambPairEnd; ++chambPairIt) {
      
      ChPairPlotterComb *plotter = chambPairIt->second;

      plotter->draw();
      delete plotter;
      
    }
  }
  
}

void L1ITMBPtComb::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtComb]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
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
	std::cout << "[L1ITMBPtComb]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
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
      
      ChPairPlotterComb *plotter = getPlotter(dttfRawId,chPairRawId);
      
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
	  ptBendOut =  PtProvider::ptLut(BendOutName.c_str(), PhiBendOut, mb1Obj,mb2Obj);    
	  
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
	  
	  
	  // print information of events with bad resolution
	  
	  //   if(((GenPt - PtM)/GenPt)>0.5){
	  
	  // 	    std::cout<<"Place  "<<PlaceName.str()<<" "<<(QualInString+QualOutString).c_str()<<std::endl<<"Mean Pt "<<PtM<<std::endl<<"Pt Gen "<<GenPt<<std::endl<<"GMT pt "<<gmtInPt<<std::endl<<"Delta phi Pt "<<ptDeltaPhi<<"  Delta Phi "<<DeltaPhi<<std::endl<<"Phi B In Pt "<<ptBendIn<<"  Phi Bend In "<<PhiBendIn<<std::endl<<"Phi B Out Pt "<<ptBendOut<<"   Phi Ben In "<<PhiBendOut<<std::endl<<std::endl;
	  
	  // 	  }
	  
	 
	  plotter->fillResolution(gmtInPt, ptBendIn, ptBendOut, ptDeltaPhi, GenPt,mb1Obj,mb2Obj,PlaceName.str()); 	          
	  
       
      }	  
    }
  }
}


int L1ITMBPtComb::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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


ChPairPlotterComb* L1ITMBPtComb::getPlotter(int dttfRawId, int chPairRawId) 
{
  // CB optimise 
  if (histos.find(dttfRawId) == histos.end()) { 
    std::cout << "[L1ITMBPtComb]::getPlotter : DTTF not in histos. Return 0.\n";
    return 0;
  }
  if (histos[dttfRawId].find(chPairRawId) == histos[dttfRawId].end()) {
    std::cout << "[L1ITMBPtComb]::getPlotter : ChambPair not in histos. Return 0.\n";
    return 0;
  }
  return histos[dttfRawId][chPairRawId];
}




#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtComb);
