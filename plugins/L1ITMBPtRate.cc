// 
// Class: L1ITMBPtRate 
//
// Info:  Performs rate plots vs pt cut for several algorithms
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


class ChPairPlotterRate {
  
public:

  ChPairPlotterRate( TFileService * fs);


  const TH1 * getHisto(std::string name) const   { return _hPlots.find(name)->second; }

  void draw() const ;

  void book(TFileService * fs);

  void fillRate(float dttfPt, float phiBInPt, float phiBOutpt, float DeltaphiPt,  Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName);


private:

  std::map<std::string, TH1 *>  _hPlots;
 
};


ChPairPlotterRate::ChPairPlotterRate(TFileService * fs)
{ 
  book(fs);   
}


void ChPairPlotterRate::fillRate(float dttfPt, float phiBInPt, float phiBOutPt, float DeltaphiPt, Int_t mb1Obj, Int_t mb2Obj, std::string PlaceName) 
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
    else if (phiBInPt<DeltaphiPt && phiBInPt<phiBOutPt && DeltaphiPt<phiBOutPt) MinPt2 = (phiBInPt+DeltaphiPt)/2.;
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
  default: std::cout<<"[L1ITMBPtRate]::SetWeight: WRONG value in weight file.\n";
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
  default: std::cout<<"[L1ITMBPtRate]::SetWeight: WRONG value in weight file.\n";
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

    if(DT1 == QualInString && DT2 == QualOutString ) break; 
  }

  for(int i = 0.; i<=120; i+=2){
    
    if(DeltaphiPt>=i)  _hPlots["DeltaPhiPtRate"]->Fill(i);
    
    
    if(dttfPt>=i) _hPlots["GMTPtRate"]->Fill(i);
    
    if(mb1Obj==1) {
      if(phiBInPt>=i) _hPlots["PhiBendPtRate"]->Fill(i);
    }
    else if (mb1Obj!=1 && mb2Obj==1){
      if(phiBOutPt>=i) _hPlots["PhiBendPtRate"]->Fill(i);
    }
    else{ 
      if(DeltaphiPt>=i)  _hPlots["PhiBendPtRate"]->Fill(i);
    }
    
    if(MinPt>=i) _hPlots["MinPtRate"]->Fill(i);
    
    if(MinPt2>=i) _hPlots["MinPt2Rate"]->Fill(i);
    
    

 //GP Best weighted average algorithm

  // BestType=BestTypeSigma;    //GP choose combination with the higher gain in RMS or sigma of a gaussian fit
  BestType=BestTypeSigma;   

    if(BestType=="DeltaPhiPtResol" ){ 
      if( DeltaphiPt>=i){
	_hPlots["BestPtRate"]->Fill(i);}     
    }
    
    else if(BestType=="PhiBInPtResol"){ 
      if( phiBInPt>=i)      _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="PhiBOutPtResol" ){ 
      if(    phiBOutPt>=i)    _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="MeanPtResol" ){ 
      if(   PtM>=i)      _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="RMSWMeanPtResol" ){
      if( PtRMSWM>=i)   _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="RMSWMeanInPtResol" ){
      if( PtRMSWMIn>=i) {
       	_hPlots["BestPtRate"]->Fill(i);}
          }
    
    else if(BestType=="RMSWMeanOutPtResol" ){
      if( PtRMSWMOut>=i) _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="SigmaWMeanPtResol" ){
      if( PtSigmaWM>=i) _hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="SigmaWMeanInPtResol" ){
      if( PtSigmaWMIn>=i) 	_hPlots["BestPtRate"]->Fill(i);}
    
    else if(BestType=="SigmaWMeanOutPtResol" ){
      if( PtSigmaWMOut>=i) _hPlots["BestPtRate"]->Fill(i);}
    
    else {
      std::cout<<"Wrong best type, check result file or the code"<<std::endl;
      std::cout<<"qual "<<  (QualInString+QualOutString).c_str()<<" place "<<PlaceName.c_str()<<" best type "<<BestType<<std::endl;
      abort();
    }    
  }
}

  
void ChPairPlotterRate::book(TFileService * fs) 
{ 



 TFileDirectory folderRate = fs->mkdir("Rate");

 _hPlots["DTTFPtDistr"] = folderRate.make<TH1F>("hDTTFPtDistr", 
						"events Pt reco > pt Gen for;pt gen ",  
						 60,-1.,119.);
 
 _hPlots["GMTPtRate"] = folderRate.make<TH1F>("hGMTPtRate", 
					      "events Pt reco > pt Gen for;pt gen ",  
					       60,-1.,119.);
 
 _hPlots["DeltaPhiPtRate"] = folderRate.make<TH1F>("hDeltaPhiPtRate", 
						   "events Pt reco > pt Gen for;pt gen ",  
						    60,-1.,119.);

 _hPlots["BestPtRate"] = folderRate.make<TH1F>("hBestPtRate", 
					       "events Pt reco > pt Gen for;pt gen ",  
					        60,-1.,119.);
 
 _hPlots["MinPtRate"] = folderRate.make<TH1F>("hMinPtRate", 
					      "events Pt reco > pt Gen for;pt gen ",  
					       60,-1.,119.);
 
 _hPlots["MinPt2Rate"] = folderRate.make<TH1F>("hMinPt2Rate", 
					      "events Pt reco > pt Gen for;pt gen ",  
					       60,-1.,119.);
 
 _hPlots["PhiBendPtRate"] = folderRate.make<TH1F>("hPhiBPtRate", 
						    "events Pt reco > pt Gen for;pt gen ",  
						     60,-1.,119.);
 }


void ChPairPlotterRate::draw() const
{

  gStyle->SetPalette(52);
  
  
  std::map<std::string, TH1 *>::const_iterator hPlotsIt = _hPlots.begin();
  std::map<std::string, TH1 *>::const_iterator hPlotsEnd = _hPlots.end();

  for (; hPlotsIt != hPlotsEnd ; ++hPlotsIt)
    {

      TH1 * hHisto = hPlotsIt->second;
 
      std::string histoName = hPlotsIt->first;
      
      std::string cName = "c" + histoName + "_2step_";
      
      system("mkdir -p plots/Rate/");
      
      TCanvas * c = new TCanvas(cName.c_str(),cName.c_str(),500,500);
      
      c->cd();
      c->SetGrid();
      
      hHisto->Draw("");

      c->SaveAs(("plots/Rate/" + cName+".pdf").c_str());
    }
}


// --------------------------------------------------
// Analyzer class, holds all ChPairPlotterRates 
// performs GMT in to GEN matching and plot booking,
// filling, drawing for all of them 
// --------------------------------------------------

class L1ITMBPtRate : public edm::EDAnalyzer {
  
public:
  
  L1ITMBPtRate( const edm::ParameterSet & );
  ~L1ITMBPtRate();

  void endJob();
  void analyze( const edm::Event &, const edm::EventSetup & );  

  float GetDeltaPhi(const L1ITMu::TriggerPrimitive * in, 
		    const L1ITMu::TriggerPrimitive * out); 
  
  float GetPhiBend(const L1ITMu::TriggerPrimitive * Prim);

  int objFromPrim(const L1ITMu::TriggerPrimitive * prim);
  
  ChPairPlotterRate* getPlotter(); 


private:

  TH1F *_hDist; 
  TH1F *_hRateDPhi; 
  TH1F *_hRateDttf; 

  edm::InputTag _mbTracksTag;  
  edm::InputTag _genParticlesTag;
  
  edm::Service<TFileService> _fs;

  std::vector<int> _mbPtChambObjects;   

  std::map< int, std::map<int, ChPairPlotterRate*> > histos; // <DTTF rawId, <ChambPair rawId, Plotter *> > 

  ChPairPlotterRate* histo;
 
  std::map<int, int> NGen;
};


L1ITMBPtRate::L1ITMBPtRate(const edm::ParameterSet& p) : 
  _mbTracksTag( p.getParameter<edm::InputTag>("MBTracksCollection") ),
  _genParticlesTag( p.getParameter<edm::InputTag>("GenParticlesCollection"))

{

  // CB: combinations to be plotted go here 
  // for now just booking DT IN, OUT CORR
  // and ch 1 and 2 of DTTF wh 1 sec 1
  
  _mbPtChambObjects.push_back(ChambPairId::DTIN);
  _mbPtChambObjects.push_back(ChambPairId::DTOUT);
  _mbPtChambObjects.push_back(ChambPairId::DTCORR);

  ChPairPlotterRate *plotter = new ChPairPlotterRate( _fs.operator->());  
  histo=plotter;

}
 


float L1ITMBPtRate::GetDeltaPhi(const L1ITMu::TriggerPrimitive * in,   
				    const L1ITMu::TriggerPrimitive * out){
 float inPhiValue  = in->getCMSGlobalPhi();  
 float outPhiValue = out->getCMSGlobalPhi();

 return  fabs(inPhiValue-outPhiValue);
}


float L1ITMBPtRate::GetPhiBend(const L1ITMu::TriggerPrimitive * Prim){
  
  return fabs(Prim->getDTData().bendingAngle);  
}
 
L1ITMBPtRate::~L1ITMBPtRate()
{

}

void L1ITMBPtRate::endJob()
{
  histo->draw();
}

void L1ITMBPtRate::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  edm::Handle<L1ITMu::MBTrackCollection> mbTracks;
  iEvent.getByLabel( _mbTracksTag, mbTracks);

  edm::Handle<reco::GenParticleCollection> gens;
  iEvent.getByLabel(_genParticlesTag, gens);

  if (!gens.isValid() || !mbTracks.isValid()) {
    std::cout << "[L1ITMBPtRate]::analyze : GEN Muons or MBTrack collections are not valid. skipping event!)\n";
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
    if (abs(wheel)==3) continue; // CB hack for test just using a given sector
    
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
	std::cout << "[L1ITMBPtRate]::analyze : WRONG quality for DTTF matched primitive. skipping.\n";
	continue;
      }

      dtBestPrims[iSt] = dtMatch.get();
      
      //  std:: cout<<"dtbestPrimis "<<dtBestPrims[iSt]<<std::endl;  //GP Add all chambers combinations
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
      
     
      ChPairPlotterRate *plotter = getPlotter();
      
      NGen[1]=1;
      
      if (plotter && gmtInPt >= 0.){ 	
	
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
	
	plotter->fillRate(gmtInPt, ptBendIn, ptBendOut, ptDeltaPhi,mb1Obj,mb2Obj,PlaceName.str()); 	          
      }	  
    } 
  }
}

					  
int L1ITMBPtRate::objFromPrim(const L1ITMu::TriggerPrimitive * prim) 
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

ChPairPlotterRate* L1ITMBPtRate::getPlotter() {

  return histo;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPtRate);
