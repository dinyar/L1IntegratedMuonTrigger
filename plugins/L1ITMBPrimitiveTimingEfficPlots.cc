/*
 * \file L1ITMBPrimitiveTimingEffPlots.cc
 * 
 * $Date: 2013/12/18 13:23:26 $
 * $Revision: 1.1 $
 * \author C.F. Bedoya, J.M. Cela - CIEMAT
 *
 */

// Framework
#include "FWCore/Framework/interface/EventSetup.h"

// DT trigger
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DQM/DTMonitorModule/interface/DTTrigGeomUtils.h"

// Geometry
#include "DataFormats/GeometryVector/interface/Pi.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <FWCore/Framework/interface/LuminosityBlock.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include <vector>
#include <string>
#include <map>

//Root
#include"TH1.h"
#include"TEfficiency.h"
#include"TAxis.h"

#include <sstream>
#include <iostream>
#include <math.h>

class L1ITMBPrimitiveTimingEffPlots: public edm::EDAnalyzer{
  
 public:
  
  /// Constructor
  L1ITMBPrimitiveTimingEffPlots(const edm::ParameterSet& ps );
  
  /// Destructor
  virtual ~L1ITMBPrimitiveTimingEffPlots();
  
 protected:
  
  // BeginJob
  void beginJob();

  ///BeginRun
  void beginRun(const edm::Run& , const edm::EventSetup&);

  /// Find best (highest qual) DCC trigger segments
  void searchDccBest(std::vector<L1MuDTChambPhDigi>* trigs);
  
  /// Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c);

  /// To reset the MEs
  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, const edm::EventSetup& context) ;
  
  /// EndJob
  void endJob(void);  

 private:  

  /// Book histos
  void bookHistos(DTChamberId chId);

 private :

  int nEvents;
  int nLumis;
  int nPhiBins, nPhibBins,nBXBins, nTPBins, nEffBins;
  double rangePhi, rangePhiB,rangeBX,rangeTP, rangeEff;
    
  int nPosBins;
  double rangePos;
  
  bool overUnderIn;

  edm::InputTag dccInputTag;
  edm::InputTag segInputTag;
 
  int trigQualBest[6][5][13];
  const L1MuDTChambPhDigi* trigBest[6][5][13];
  bool track_ok[6][5][15]; // CB controlla se serve

  edm::Service<TFileService> fs;

  edm::ParameterSet parameters;
  edm::ESHandle<DTGeometry> muonGeom;
  std::string theGeomLabel;
  DTTrigGeomUtils* trigGeomUtils;

  std::map<uint32_t, std::map<std::string, TH1F*> > chHistos;
  std::map<uint32_t, std::map<std::string, TEfficiency*> > EffchHistos;
  std::map<std::string, TH1F*> sHistos;
  std::map<std::string, TH1F*> sHistos2;
  std::map<std::string, TH1F*> sHistos3;
   
};


using namespace edm;
using namespace std;

L1ITMBPrimitiveTimingEffPlots::L1ITMBPrimitiveTimingEffPlots(const edm::ParameterSet& ps) : trigGeomUtils(0) {
	
  LogTrace("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: Constructor"<<endl;
	
  dccInputTag  = ps.getUntrackedParameter<InputTag>("inputTagDCC");
  segInputTag  = ps.getUntrackedParameter<InputTag>("inputTagSEG");

  overUnderIn      = ps.getUntrackedParameter<bool>("rebinOutFlowsInGraph");
  theGeomLabel     = ps.getUntrackedParameter<string>("geomLabel");
  
  
  nBXBins  = 16;
  rangeBX  = 8;

  nTPBins  = 15;
  rangeTP  = 15;
  
  nEffBins  = 200;
  rangeEff  = 200;
  
  parameters = ps;
	
}


L1ITMBPrimitiveTimingEffPlots::~L1ITMBPrimitiveTimingEffPlots() {
  
  LogTrace("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: analyzed " << nEvents << " events" << endl;
  if (trigGeomUtils) { delete trigGeomUtils; }
	
}


void L1ITMBPrimitiveTimingEffPlots::beginJob(){
	
  LogTrace("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: BeginJob" << endl;
  nEvents = 0;
  nLumis  = 0;
	
}

void L1ITMBPrimitiveTimingEffPlots::bookHistos(DTChamberId chId) {

  stringstream wheel; wheel << chId.wheel();	
  stringstream sector; sector << chId.sector();
  stringstream station; station << chId.station();	

  TFileDirectory folder  = fs->mkdir(("Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str()).c_str());
  
  string chTag = "_W" + wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
  std::map<std::string, TH1F*> &chambMap = chHistos[chId.rawId()];
  std::map<std::string, TEfficiency*> &effchambMap = EffchHistos[chId.rawId()];
  
  
  string hName = "Timing_DTprim_BX_ALLQ";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for ALL qualities",nBXBins,-rangeBX,rangeBX);
/*
  hName = "Timing_DTprim_BX_Q03";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q<4",nBXBins,-rangeBX,rangeBX);
   
  hName = "Timing_DTprim_BX_Q46";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q>=4",nBXBins,-rangeBX,rangeBX);
*/
  hName = "Timing_DTprim_BX_Q0";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=0",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q1";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=1",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q2";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=2",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q3";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=3",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q4";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=4",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q5";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=5",nBXBins,-rangeBX,rangeBX);

    hName = "Timing_DTprim_BX_Q6";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"BX distribution of DT trigger primitives for Q=6",nBXBins,-rangeBX,rangeBX);

  
  hName = "Clasific_OOT";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Number of entries for different clasification",5,0,4);

   hName = "Timing_DTprim_BX_outoftime";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Number of entries of out of time DT trigger primitives ALL qualities",nBXBins,-rangeBX,rangeBX);

    hName = "Number_TPrimitives";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Number of trigger primitives per event",nTPBins,0,rangeTP);

   hName = "TP_Efficiency_posPhi";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_posPhi",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_posEta";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_posEta",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_posPhi_BX0";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_posPhi_BX0",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_posEta_BX0";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_posEta_BX0",nEffBins,-rangeEff,rangeEff);
  
     hName = "TP_Efficiency_dirPhi";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_dirPhi",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_dirEta";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_dirEta",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_dirPhi_BX0";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_dirPhi_BX0",nEffBins,-rangeEff,rangeEff);

   hName = "TP_Efficiency_dirEta_BX0";
  effchambMap[hName] = folder.make<TEfficiency>((hName+chTag).c_str(),"TP_Efficiency_dirEta_BX0",nEffBins,-rangeEff,rangeEff);
  
}

  

void L1ITMBPrimitiveTimingEffPlots::beginRun(const edm::Run& run, const edm::EventSetup& context) {
	
  LogTrace("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: BeginRun" << endl;   
	
  context.get<MuonGeometryRecord>().get(theGeomLabel,muonGeom);
  trigGeomUtils = new DTTrigGeomUtils(muonGeom);
	
  std::vector<DTChamber*>::const_iterator chambIt  = muonGeom->chambers().begin();
  std::vector<DTChamber*>::const_iterator chambEnd = muonGeom->chambers().end();
  
  for (; chambIt!=chambEnd; ++chambIt)
    bookHistos((*chambIt)->id());

   TFileDirectory crisfolder  = fs->mkdir("Overall");
   std::map<std::string, TH1F*> &crischambMap = sHistos;
   crischambMap["ALL_Timing_DTprim_BX_ALLQ"]=crisfolder.make<TH1F>("ALL_Timing_DTprim_BX_ALLQ","Summary BX",nBXBins,-rangeBX,rangeBX);
   /*
   std::map<std::string, TH1F*> &crischambMap2 = sHistos2;
   crischambMap2["ALL_Timing_DTprim_BX_Q03"]=crisfolder.make<TH1F>("ALL_Timing_DTprim_BX_Q03","Summary BX Q<4",nBXBins,-rangeBX,rangeBX);
   std::map<std::string, TH1F*> &crischambMap3 = sHistos3;
   crischambMap3["ALL_Timing_DTprim_BX_Q46"]=crisfolder.make<TH1F>("ALL_Timing_DTprim_BX_Q46","Summary BX Q>=4",nBXBins,-rangeBX,rangeBX);
   */
}


void L1ITMBPrimitiveTimingEffPlots::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {

  nLumis++;
  LogTrace("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: Begin of LS transition" << endl;
  
}


void L1ITMBPrimitiveTimingEffPlots::endJob(){
	
  LogVerbatim("L1ITMBPrimitiveTimingEffPlots") << "[L1ITMBPrimitiveTimingEffPlots]: analyzed " << nEvents << " events" << endl;
	
}


void L1ITMBPrimitiveTimingEffPlots::analyze(const edm::Event& e, const edm::EventSetup& c){
	
  nEvents++;
    
  edm::Handle<L1MuDTChambPhContainer> trigHandle;
  e.getByLabel(dccInputTag,trigHandle);
  vector<L1MuDTChambPhDigi>* trigs = trigHandle->getContainer();
  searchDccBest(trigs);

  Handle<DTRecSegment4DCollection> segments4D;
  e.getByLabel(segInputTag,segments4D);  		
  DTRecSegment4DCollection::id_iterator chamberId;

  // Preliminary loop finds best 4D Segment and high quality ones
  vector<const DTRecSegment4D*> best4DSegments;
  for (chamberId = segments4D->id_begin(); chamberId != segments4D->id_end(); ++chamberId){
		
    DTRecSegment4DCollection::range  rangeInCh = segments4D->get(*chamberId);
    DTRecSegment4DCollection::const_iterator trackIt  = rangeInCh.first;
    DTRecSegment4DCollection::const_iterator trackEnd = rangeInCh.second;

    const DTRecSegment4D* tmpBest = 0;
    int tmpdof = 0;
    int dof = 0;
    
    for (; trackIt!=trackEnd; ++trackIt){
      if(trackIt->hasPhi()) {				
	dof = trackIt->phiSegment()->degreesOfFreedom();
	if (dof>tmpdof) {
	  tmpBest = &(*trackIt);
	  tmpdof = dof;	
	}
      }
 
    }
    if (tmpBest) best4DSegments.push_back(tmpBest);
  } //end loop search best 4D segments

  vector<const DTRecSegment4D*>::const_iterator bestTrackIt  = best4DSegments.begin();
  vector<const DTRecSegment4D*>::const_iterator bestTrackEnd = best4DSegments.end();
  
  //loop on best quality segments
  for (; bestTrackIt!=bestTrackEnd; ++bestTrackIt) { 
    if((*bestTrackIt)->hasPhi()) {
      
      DTChamberId chId = (*bestTrackIt)->chamberId();
      int nHitsPhi = (*bestTrackIt)->phiSegment()->degreesOfFreedom()+2;
      int wheel    = chId.wheel();
      int station  = chId.station();
      int scsector = 0; //conputeSCCoordinates updates it
      float trackPosPhi, trackPosEta, trackDirPhi, trackDirEta;
      trigGeomUtils->computeSCCoordinates((*bestTrackIt),scsector,trackPosPhi,trackDirPhi,trackPosEta,trackDirEta);
      //trackPosPhi is position inside the chamber obtained from the DT segment, from -150 to +150 in an MB3
      //trigGeomUtils->trigToSeg(station,trigPos,trackDirPhi); La contraria es la q hace falta para phase 2
	
      map<string, TH1F*> &chMap = chHistos[chId.rawId()];
      map<string, TEfficiency*> &effchMap = EffchHistos[chId.rawId()];

      map<string, TH1F*> &crischMap = sHistos;
 //     map<string, TH1F*> &crischMap2 = sHistos2;
  //    map<string, TH1F*> &crischMap3 = sHistos3;
      bool TP_good=false;
      bool TP_good_BX0=false;
      
      int BXbest=-1;
      if (trigQualBest[wheel+3][station][scsector] > -1 &&  //fills up the BX histogram for all qualities
	  trigQualBest[wheel+3][station][scsector] < 7  &&
	  nHitsPhi>=7 ) {
	BXbest=trigBest[wheel+3][station][scsector]->bxNum();
	chMap.find("Timing_DTprim_BX_ALLQ")->second->Fill(BXbest);
	crischMap.find("ALL_Timing_DTprim_BX_ALLQ")->second->Fill(BXbest);
	} //trigQualBest if

      BXbest=-1;
      if (trigQualBest[wheel+3][station][scsector] > -1 &&  //EFFICIENCY tests
	  trigQualBest[wheel+3][station][scsector] < 7 ) {
	   // && nHitsPhi>=7 ) {
	BXbest=trigBest[wheel+3][station][scsector]->bxNum();
	TP_good=true;
	if (BXbest==0) TP_good_BX0=true;
	  
	  } //trigQualBest if
	  
	effchMap.find("TP_Efficiency_posPhi")->second->Fill(TP_good,trackPosPhi);
	effchMap.find("TP_Efficiency_posEta")->second->Fill(TP_good,trackPosEta);
	effchMap.find("TP_Efficiency_posPhi_BX0")->second->Fill(TP_good_BX0,trackPosPhi);
	effchMap.find("TP_Efficiency_posEta_BX0")->second->Fill(TP_good_BX0,trackPosEta);
	effchMap.find("TP_Efficiency_dirPhi")->second->Fill(TP_good,trackDirPhi);
	effchMap.find("TP_Efficiency_dirEta")->second->Fill(TP_good,trackDirEta);
	effchMap.find("TP_Efficiency_dirPhi_BX0")->second->Fill(TP_good_BX0,trackDirPhi);
	effchMap.find("TP_Efficiency_dirEta_BX0")->second->Fill(TP_good_BX0,trackDirEta);
     
	//} //trigQualBest if

/*
      BXbest=-1;
      if (trigQualBest[wheel+3][station][scsector] > -1 &&   //fills up the BX histogram for qualities 0 to 3
	  trigQualBest[wheel+3][station][scsector] < 4  &&
	  nHitsPhi>=7 ) {
	BXbest=trigBest[wheel+3][station][scsector]->bxNum();
	chMap.find("Timing_DTprim_BX_Q03")->second->Fill(BXbest);
	crischMap2.find("ALL_Timing_DTprim_BX_Q03")->second->Fill(BXbest);
      } //trigQualBest if

      BXbest=-1;
      if (trigQualBest[wheel+3][station][scsector] > 3 &&   //fills up the BX histogram for qualities 4 to 6
	  trigQualBest[wheel+3][station][scsector] < 7  &&
	  nHitsPhi>=7 ) {
	BXbest=trigBest[wheel+3][station][scsector]->bxNum();
	chMap.find("Timing_DTprim_BX_Q46")->second->Fill(BXbest);
	crischMap3.find("ALL_Timing_DTprim_BX_Q46")->second->Fill(BXbest);
      } //trigQualBest if
*/
//individual qualities      
      BXbest=-1;
      if (nHitsPhi>=7 ) {
          if (trigQualBest[wheel+3][station][scsector] == 0) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q0")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 1) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q1")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 2) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q2")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 3) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q3")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 4) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q4")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 5) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q5")->second->Fill(BXbest);
          }
          if (trigQualBest[wheel+3][station][scsector] == 6) {
          	  BXbest=trigBest[wheel+3][station][scsector]->bxNum();
          	  chMap.find("Timing_DTprim_BX_Q6")->second->Fill(BXbest);
          }
      } //trigQualBest if individual qualities
 
  vector<L1MuDTChambPhDigi>::const_iterator trigIt  = trigs->begin();
  vector<L1MuDTChambPhDigi>::const_iterator trigEnd = trigs->end();
  int NumTP=0;
  int BX_i=-1;
  for(; trigIt!=trigEnd; ++trigIt) {
    int wh   = trigIt->whNum();
    int sec  = trigIt->scNum() + 1; // DTTF -> DT sector range transform
    int st   = trigIt->stNum();
     
  //  int qual = trigIt->code();
  if ((wh==wheel) && (sec==scsector) && (st==station) && (nHitsPhi>7)){
  	   BX_i=trigIt->bxNum();
     	   NumTP=NumTP+1;
	   chMap.find("Timing_DTprim_BX_outoftime")->second->Fill(BX_i);
   }//if matching same station
  }//for trigs
 
    chMap.find("Number_TPrimitives")->second->Fill(NumTP);
 
    //Clasification OOT TP
    if (NumTP==1) {
    	    //select if it has the BXin =0
    	    if (BX_i==0){    chMap.find("Clasific_OOT")->second->Fill(0);
	    }else{    chMap.find("Clasific_OOT")->second->Fill(1);
	    }
    }else {//nuntp=1
    	    //search for the best quality 
    	if (BXbest!=-1){
    	    //select if it has the BXin =0
    	    if (BXbest==0){  chMap.find("Clasific_OOT")->second->Fill(2);
	    }else{chMap.find("Clasific_OOT")->second->Fill(3);
	    }
	}
    }//nuntp=1
    
    
   
    } //if hasPhi
  }//for bestTracks 
  
}//::analyze

void L1ITMBPrimitiveTimingEffPlots::searchDccBest( std::vector<L1MuDTChambPhDigi>* trigs ){
  
  string histoType ;
  string histoTag ;
  
  // define best quality trigger segment
  // start from 1 and zero is kept empty
  for (int st=0;st<=4;++st)
    for (int wh=0;wh<=5;++wh)
      for (int sec=0;sec<=12;++sec)
	trigQualBest[wh][st][sec] = -1;    
	
  vector<L1MuDTChambPhDigi>::const_iterator trigIt  = trigs->begin();
  vector<L1MuDTChambPhDigi>::const_iterator trigEnd = trigs->end();
  for(; trigIt!=trigEnd; ++trigIt) {
    
    int wh   = trigIt->whNum();
    int sec  = trigIt->scNum() + 1; // DTTF -> DT sector range transform
    int st   = trigIt->stNum();
    int qual = trigIt->code();

    if(qual>trigQualBest[wh+3][st][sec] && qual<7) {
      trigQualBest[wh+3][st][sec]=qual; 
      trigBest[wh+3][st][sec] = &(*trigIt);
    }
    
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPrimitiveTimingEffPlots);


