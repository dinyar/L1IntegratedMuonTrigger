// 
// Class: L1ITMBPtEff 
//
// Info:   Performs Efficency plot vs pt  for different algorithms 
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


class ChPairPlotterEff {
  
public:

  ChPairPlotterEff(int wh, int sec, int inCh, int outCh, 
		   int inChObj, int outChObj, TFileService * fs);

  ChPairPlotterEff(ChambPairId id, TFileService * fs);

  ChPairPlotterEff(const ChPairPlotterEff & plotter);

  int rawId() const { return _id.rawId(); };
  std::string name() const { return _id.name(); };

  const TH1 * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  const DTTFId & dttfId() const { return _id.dttfId(); };

  void draw() const ;

  void book(TFileService * fs);


  void fillEfficiency(float dttfPt, float phiBInPt, float phiBOutpt, float DeltaphiPt, float pt,  Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);

private:

  ChambPairId _id;

  std::vector<int> _dttfPtCuts;
  std::map<std::string, TH1 *>  _hPlots;
  std::map<std::string, TEfficiency *> _ePlots;


};


ChPairPlotterEff::ChPairPlotterEff(int wh, int sec, int inCh, int outCh, 
				   int inChObj, int outChObj, TFileService * fs) :
  _id(wh,sec,inCh,outCh,inChObj,outChObj) 
{ 
  
  // CB DTTF pt cuts for plots (now one per object)  
  
  _dttfPtCuts.push_back(16);
  _dttfPtCuts.push_back(18);
  _dttfPtCuts.push_back(20);
  _dttfPtCuts.push_back(22);
  _dttfPtCuts.push_back(25); 
  
  book(fs); 
  
}


ChPairPlotterEff::ChPairPlotterEff(ChambPairId id, TFileService * fs) : _id(id) 
{ 
  
  // CB DTTF pt cuts for plots (now one per object)  
  
  _dttfPtCuts.push_back(16);
  _dttfPtCuts.push_back(18);
  _dttfPtCuts.push_back(20);
  _dttfPtCuts.push_back(22);
  _dttfPtCuts.push_back(25); 
  
  book(fs); 
  
}


ChPairPlotterEff::ChPairPlotterEff(const ChPairPlotterEff & plotter) :
 _id(plotter._id), _dttfPtCuts(plotter._dttfPtCuts),
 _hPlots(plotter._hPlots)
{

}


void ChPairPlotterEff::fillEfficiency(float dttfPt, float phiBInPt, float phiBOutPt, float DeltaphiPt, float pt, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName) 
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

  Float_t PtM = (phiBOutPt+phiBInPt+DeltaphiPt)/3.;
  Float_t PtRMSWM = (phiBOutPt*whWMRMSBOut+phiBInPt*whWMRMSBIn+DeltaphiPt*whWMRMSDPhi); 
  Float_t PtRMSWMIn = (phiBInPt*whWMInRMSBIn+DeltaphiPt*whWMInRMSDPhi);
  Float_t PtRMSWMOut = (phiBOutPt*whWMOutRMSBOut+DeltaphiPt*whWMOutRMSDPhi); 
  Float_t PtSigmaWM  = (phiBOutPt*whWMSigmaBOut+phiBInPt*whWMSigmaBIn+DeltaphiPt*whWMSigmaDPhi);
  Float_t PtSigmaWMIn = (phiBInPt*whWMInSigmaBIn+DeltaphiPt*whWMInSigmaDPhi);
  Float_t PtSigmaWMOut = (phiBOutPt*whWMOutSigmaBOut+DeltaphiPt*whWMOutSigmaDPhi);

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
  default: std::cout<<"[L1ITMBPtEff]::SetWeight: WRONG value in weight file.\n";
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
  default: std::cout<<"[L1ITMBPtEff]::SetWeight: WRONG value in weight file.\n";
    break;
  }
  

  std::string FileName = "BestResult/"+PlaceName+".txt";
  std::ifstream datafileIn(FileName.c_str(),std::ifstream::in);
  std::string line;
  std::string BestTypeSigma,BestTypeRMS,BestType;

  
  while (getline(datafileIn,line)){
    std::stringstream readline; 
    std::string DT1,DT2,A;
    readline<<line;
    readline>>DT1>>DT2>>A>>BestTypeRMS>>A>>BestTypeSigma;
 
    if(DT1 != QualInString && DT2 != QualOutString ) continue;

    break; 
  }

  //GP choose combination with the higher gain in RMS or sigma of a gaussian fit
   BestType=BestTypeSigma;
  // BestType=BestTypeRMS;
  
  std::vector<int>::const_iterator dttfPtCutIt = _dttfPtCuts.begin();
  std::vector<int>::const_iterator dttfPtCutEnd = _dttfPtCuts.end();


  
   for(; dttfPtCutIt != dttfPtCutEnd; ++dttfPtCutIt){
      int dttfPtCut = (*dttfPtCutIt);
      
      std::stringstream ptCut; ptCut << dttfPtCut;
      std::string dttfPtTag = ptCut.str();

  
      _ePlots["MinPtEff"+ dttfPtTag]->Fill(MinPt>=dttfPtCut,pt);

      _ePlots["MinPt2Eff"+ dttfPtTag]->Fill(MinPt2>=dttfPtCut,pt);

      _ePlots["DeltaPhiPtEff"+ dttfPtTag]->Fill(DeltaphiPt>=dttfPtCut,pt);

      _ePlots["GMTPtEff"+ dttfPtTag]->Fill(dttfPt>=dttfPtCut,pt);

      _ePlots["PhiBendInPtEff"+ dttfPtTag]->Fill(phiBInPt>=dttfPtCut,pt);

      _ePlots["PhiBendOutPtEff"+ dttfPtTag]->Fill(phiBOutPt>=dttfPtCut,pt);


     
      if(BestType=="DeltaPhiPtResol")  _ePlots["BestPtEff"+dttfPtTag]->Fill(DeltaphiPt>=dttfPtCut,pt);
      
      else if(BestType=="PhiBInPtResol")     _ePlots["BestPtEff"+dttfPtTag]->Fill(phiBInPt>=dttfPtCut,pt);
      
      else if(BestType=="PhiBOutPtResol")    _ePlots["BestPtEff"+dttfPtTag]->Fill(phiBOutPt>=dttfPtCut,pt);
      
      else if(BestType=="MeanPtResol")       _ePlots["BestPtEff"+dttfPtTag]->Fill(PtM>=dttfPtCut,pt);
      
      else if(BestType=="RMSWMeanPtResol")   _ePlots["BestPtEff"+dttfPtTag]->Fill(PtRMSWM>=dttfPtCut,pt);

      else if(BestType=="RMSWMeanInPtResol")   _ePlots["BestPtEff"+dttfPtTag]->Fill(PtRMSWMIn>=dttfPtCut,pt);
       									    
      else if(BestType=="RMSWMeanOutPtResol") _ePlots["BestPtEff"+dttfPtTag]->Fill(PtRMSWMOut>=dttfPtCut,pt);
									    
      else if(BestType=="SigmaWMeanPtResol")  _ePlots["BestPtEff"+dttfPtTag]->Fill(PtSigmaWM>=dttfPtCut,pt);

      else if(BestType=="SigmaWMeanInPtResol") 	_ePlots["BestPtEff"+dttfPtTag]->Fill(PtSigmaWMIn>=dttfPtCut,pt);
      
      else if(BestType=="SigmaWMeanOutPtResol")  _ePlots["BestPtEff"+dttfPtTag]->Fill(PtSigmaWMOut>=dttfPtCut,pt);
      
      else {
	std::cout<<"Wrong best type, check result file or the code"<<std::endl;
	std::cout<<"qual "<<  (QualInString+QualOutString).c_str()<<" place "<<PlaceName.c_str()<<" best type = "<<BestType.c_str()<<std::endl;
	abort();
      }         
   }
}


void ChPairPlotterEff::book(TFileService * fs) 
{ 

  std::string hDir  = dttfId().name();
  std::string hName = name();
   
  TFileDirectory baseFolder = fs->mkdir(hDir.c_str());
  
  TFileDirectory folderEff = baseFolder.mkdir("Efficiency");
  
  std::vector<int>::const_iterator dttfPtCutIt = _dttfPtCuts.begin();
  std::vector<int>::const_iterator dttfPtCutEnd = _dttfPtCuts.end();
  

  for (; dttfPtCutIt != dttfPtCutEnd; ++dttfPtCutIt)
    {
      int dttfPtCut = (*dttfPtCutIt);
      
      std::stringstream ptCut; ptCut << dttfPtCut;
      std::string dttfPtTag = ptCut.str();
      
      _ePlots["GMTPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hGMTPtEff"+hName+ dttfPtTag).c_str(), 
								   ("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
								  60,-1.,119.);
      
      _ePlots["BestPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hBestPtEff"+hName+ dttfPtTag).c_str(), 
								    ("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
								   60,-1.,119.);
      
      _ePlots["MinPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hMinPtEff"+hName+ dttfPtTag).c_str(), 
								   ("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
								  60,-1.,119.);
      
      _ePlots["MinPt2Eff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hMinPt2Eff"+hName+ dttfPtTag).c_str(), 
								   ("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
								   60,-1.,119.);
      
      _ePlots["DeltaPhiPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hDeltaPhiPtEff"+hName+ dttfPtTag).c_str(), 
									("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
									60,-1.,119.);

      _ePlots["PhiBendInPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hPhiBInPtEff"+hName+ dttfPtTag).c_str(), 
									("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
									60,-1.,119.);
      
      _ePlots["PhiBendOutPtEff"+ dttfPtTag] = folderEff.make<TEfficiency>(("hPhiBOutPtEff"+hName+ dttfPtTag).c_str(), 
									("events Pt reco > pt Gen for"+hName+ dttfPtTag+";pt gen ").c_str(),  
									60,-1.,119.);

 } 

  
  TFileDirectory folderRes = baseFolder.mkdir("Resolutions");
  
 

  _hPlots["BestPtResolvsPt"] = folderRes.make<TH2F>(("hBestPhiPtResolvsPt" + hName ).c_str(), 
							("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							 ";GEN  mu p_{T};(GEN  mu p_{T} - Best pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
							60, 0.,120., 120, -2.5, 2.5);

  _hPlots["MinPtResolvsPt"] = folderRes.make<TH2F>(("hMinPhiPtResolvsPt" + hName ).c_str(), 
							("obj (dttf_[pt] - pt) / pt vs pt for " + hName + 
							 ";GEN  mu p_{T};(GEN  mu p_{T} - Best pt mu p_{T}) / GEN  mu p_{T}").c_str(), 
							60, 0.,120., 120, -2.5, 2.5);
}


void ChPairPlotterEff::draw() const
{

  gStyle->SetPalette(52);
  
  
  std::map<std::string, TH1 *>::const_iterator hPlotsIt = _hPlots.begin();
  std::map<std::string, TH1 *>::const_iterator hPlotsEnd = _hPlots.end();

  for (; hPlotsIt != hPlotsEnd ; ++hPlotsIt)
    {

      TH1 * hHisto = hPlotsIt->second;
 
      std::string histoName = hPlotsIt->first;
      
      std::string cName = "c" + histoName + "_2step_" + name();
      
      system(std::string("mkdir -p plots/Efficiency/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();
      
      hHisto->Draw("");

      c->SaveAs(("plots/Efficiency/" + dttfId().name() + "/" + cName+".pdf").c_str());
    }


 std::map<std::string, TEfficiency *>::const_iterator ePlotsIt  =  _ePlots.begin();
  std::map<std::string, TEfficiency *>::const_iterator ePlotsEnd =  _ePlots.end();

  for (; ePlotsIt != ePlotsEnd ; ++ePlotsIt)
    {

      TEfficiency * histo = ePlotsIt->second;
      std::string histoName = ePlotsIt->first;
      
      std::string cName = "c" + histoName + name();
      
      system(std::string("mkdir -p plots/Efficiency/" + dttfId().name()).c_str());
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();  
      
      histo->Draw("E0PA");
      
      c->SaveAs(("plots/Efficiency/" + dttfId().name() + "/" + cName+".pdf").c_str());
    }
}


// --------------------------------------------------
// Analyzer class, holds all ChPairPlotterEffs 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtEff : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtEff( const edm::ParameterSet & );
  ~L1ITMBPtEff();

  void endJob();
  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  ChPairPlotterEff* getPlotter(int dttfRawId, int chPairRawId); 


private:

  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChPairPlotterEff*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

  std::map<int, int> NGen;
};


L1ITMBPtEff::L1ITMBPtEff(const edm::ParameterSet& p) : 
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
	      ChPairPlotterEff *plotter = new ChPairPlotterEff(id,_fs.operator->());  

	      histos[DTTFId(wheel,sector).rawId()][id.rawId()] = plotter;      
	    }
	  }
	}
      }
    }
  }  
}
 
float L1ITMBPtEff::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}



float L1ITMBPtEff::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);  
}

L1ITMBPtEff::~L1ITMBPtEff()
{

}

void L1ITMBPtEff::endJob()
{
  
  std::map<int,std::map<int,ChPairPlotterEff*> >::iterator dttfIt = histos.begin();
  std::map<int,std::map<int,ChPairPlotterEff*> >::iterator dttfEnd = histos.end();
  
  for(; dttfIt!=dttfEnd; ++dttfIt) {
    
    std::map<int,ChPairPlotterEff*>::iterator chambPairIt = dttfIt->second.begin();
    std::map<int,ChPairPlotterEff*>::iterator chambPairEnd = dttfIt->second.end();
    
    for(; chambPairIt!=chambPairEnd; ++chambPairIt) {
      
      ChPairPlotterEff *plotter = chambPairIt->second;

      plotter->draw();
      delete plotter;
      
    }
  }
  
}

void L1ITMBPtEff::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtEff]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
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
	std::cout << "[L1ITMBPtEff]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
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
      
      ChPairPlotterEff *plotter = getPlotter(dttfRawId,chPairRawId);
      
      NGen[1]=1;

      if (plotter){

	Float_t GenPt  = bestGen->pt();

	
	if (gmtInPt >= 0.){ 	
	  
	  Float_t DeltaPhi =  GetDeltaPhi(dtBestPrims[0],dtBestPrims[1]);  
	  Float_t PhiBendIn =  GetPhiBend(dtBestPrims[0]); 
	  Float_t PhiBendOut =  GetPhiBend(dtBestPrims[1]);   
	  
	  Float_t ptBendIn = -1;
	  Float_t ptBendOut = -1;
	  Float_t ptDeltaPhi = -1;
	  
	  std::string BendInName = "PhiBendIn"+PlaceName.str();
	  ptBendIn = PtProvider::ptLut(BendInName.c_str(), PhiBendIn, mb1Obj,mb2Obj);    
	  
	  std::string BendOutName = "PhiBendOut"+PlaceName.str();
	  ptBendOut =  PtProvider::ptLut(BendOutName.c_str(), PhiBendOut, mb1Obj,mb2Obj);    
	  
	  std::string DeltaPhiName = "Phi"+PlaceName.str();
	  ptDeltaPhi =  PtProvider::ptLut(DeltaPhiName.c_str(), DeltaPhi, mb1Obj,mb2Obj);    		
	  
	  plotter->fillEfficiency(gmtInPt, ptBendIn, ptBendOut, ptDeltaPhi, GenPt,mb1Obj,mb2Obj,PlaceName.str()); 	          
	  
	}
      }	  
    }
  }
}


int L1ITMBPtEff::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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


ChPairPlotterEff* L1ITMBPtEff::getPlotter(int dttfRawId, int chPairRawId) 
{
  // CB optimise 
  if (histos.find(dttfRawId) == histos.end()) { 
    std::cout << "[L1ITMBPtEff]::getPlotter : DTTF not in histos. Return 0.\n";
    return 0;
  }
  if (histos[dttfRawId].find(chPairRawId) == histos[dttfRawId].end()) {
    std::cout << "[L1ITMBPtEff]::getPlotter : ChambPair not in histos. Return 0.\n";
    return 0;
  }
  return histos[dttfRawId][chPairRawId];
}




#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtEff);
