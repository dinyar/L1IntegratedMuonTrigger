/*
 * \file L1ITMBPrimitiveTimingPlots.cc
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
#include"TAxis.h"

#include <sstream>
#include <iostream>
#include <math.h>

class L1ITMBPrimitiveTimingPlots: public edm::EDAnalyzer{
  
 public:
  
  /// Constructor
  L1ITMBPrimitiveTimingPlots(const edm::ParameterSet& ps );
  
  /// Destructor
  virtual ~L1ITMBPrimitiveTimingPlots();
  
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
  int nPhiBins, nPhibBins,nBXBins;
  double rangePhi, rangePhiB,rangeBX;
  
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
  std::map<std::string, TH1F*> sHistos;
   
};


using namespace edm;
using namespace std;

L1ITMBPrimitiveTimingPlots::L1ITMBPrimitiveTimingPlots(const edm::ParameterSet& ps) : trigGeomUtils(0) {
	
  LogTrace("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: Constructor"<<endl;
	
  dccInputTag  = ps.getUntrackedParameter<InputTag>("inputTagDCC");
  segInputTag  = ps.getUntrackedParameter<InputTag>("inputTagSEG");

  overUnderIn      = ps.getUntrackedParameter<bool>("rebinOutFlowsInGraph");
  theGeomLabel     = ps.getUntrackedParameter<string>("geomLabel");
  
  
  nBXBins  = 10;
  rangeBX  = 5;
  
  parameters = ps;
	
}


L1ITMBPrimitiveTimingPlots::~L1ITMBPrimitiveTimingPlots() {
  
  LogTrace("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: analyzed " << nEvents << " events" << endl;
  if (trigGeomUtils) { delete trigGeomUtils; }
	
}


void L1ITMBPrimitiveTimingPlots::beginJob(){
	
  LogTrace("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: BeginJob" << endl;
  nEvents = 0;
  nLumis  = 0;
	
}

void L1ITMBPrimitiveTimingPlots::bookHistos(DTChamberId chId) {

  stringstream wheel; wheel << chId.wheel();	
  stringstream sector; sector << chId.sector();
  stringstream station; station << chId.station();	

  TFileDirectory folder  = fs->mkdir(("Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str()).c_str());
  
  string chTag = "_W" + wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
  std::map<std::string, TH1F*> &chambMap = chHistos[chId.rawId()];

  
  string hName = "Timing_DTprim_BX";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Number of entries of DT trigger primitives versus BX",nBXBins,-rangeBX,rangeBX);

}

  

void L1ITMBPrimitiveTimingPlots::beginRun(const edm::Run& run, const edm::EventSetup& context) {
	
  LogTrace("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: BeginRun" << endl;   
	
  context.get<MuonGeometryRecord>().get(theGeomLabel,muonGeom);
  trigGeomUtils = new DTTrigGeomUtils(muonGeom);
	
  std::vector<DTChamber*>::const_iterator chambIt  = muonGeom->chambers().begin();
  std::vector<DTChamber*>::const_iterator chambEnd = muonGeom->chambers().end();
  
  for (; chambIt!=chambEnd; ++chambIt)
    bookHistos((*chambIt)->id());

   TFileDirectory crisfolder  = fs->mkdir("Overall");
   std::map<std::string, TH1F*> &crischambMap = sHistos;
   crischambMap["ALL_Timing_DTprim_BX"]=crisfolder.make<TH1F>("ALL_Timing_DTprim_BX","Summary BX",nBXBins,-rangeBX,rangeBX);

}


void L1ITMBPrimitiveTimingPlots::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {

  nLumis++;
  LogTrace("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: Begin of LS transition" << endl;
  
}


void L1ITMBPrimitiveTimingPlots::endJob(){
	
  LogVerbatim("L1ITMBPrimitiveTimingPlots") << "[L1ITMBPrimitiveTimingPlots]: analyzed " << nEvents << " events" << endl;
	
}


void L1ITMBPrimitiveTimingPlots::analyze(const edm::Event& e, const edm::EventSetup& c){
	
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
  
  }

  vector<const DTRecSegment4D*>::const_iterator bestTrackIt  = best4DSegments.begin();
  vector<const DTRecSegment4D*>::const_iterator bestTrackEnd = best4DSegments.end();
  
  for (; bestTrackIt!=bestTrackEnd; ++bestTrackIt) {
    
    if((*bestTrackIt)->hasPhi()) {
      
      DTChamberId chId = (*bestTrackIt)->chamberId();
      int nHitsPhi = (*bestTrackIt)->phiSegment()->degreesOfFreedom()+2;
      
      int wheel    = chId.wheel();
      int station  = chId.station();
      int scsector = 0;
      float trackPosPhi, trackPosEta, trackDirPhi, trackDirEta;
      trigGeomUtils->computeSCCoordinates((*bestTrackIt),scsector,trackPosPhi,trackDirPhi,trackPosEta,trackDirEta);
      
      map<string, TH1F*> &chMap = chHistos[chId.rawId()];

      if (trigQualBest[wheel+3][station][scsector] > -1 &&  // residuals only for correlate triggers
	  trigQualBest[wheel+3][station][scsector] < 7  &&
	  nHitsPhi>=7 ) {

/*					
	float trigPos = trigGeomUtils->trigPos(trigBest[wheel+3][station][scsector]);
	float trigDir = trigGeomUtils->trigDir(trigBest[wheel+3][station][scsector]);
	trigGeomUtils->trigToSeg(station,trigPos,trackDirPhi);
										
	double deltaPos = trigPos-trackPosPhi;
	deltaPos = overUnderIn ? max(min(deltaPos,rangePhi-0.01),-rangePhi+0.01) : deltaPos;
	double deltaDir = trigDir-trackDirPhi;
	deltaDir = overUnderIn ? max(min(deltaDir,rangePhiB-0.01),-rangePhiB+0.01) : deltaDir;
*/	
	
	int BXbest=0;
	BXbest=trigBest[wheel+3][station][scsector]->bxNum();
	chMap.find("Timing_DTprim_BX")->second->Fill(BXbest);

	map<string, TH1F*> &crischMap = sHistos;

	crischMap.find("ALL_Timing_DTprim_BX")->second->Fill(BXbest);
	
      }
      
    }
  } 
  
}

void L1ITMBPrimitiveTimingPlots::searchDccBest( std::vector<L1MuDTChambPhDigi>* trigs ){
  
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
DEFINE_FWK_MODULE(L1ITMBPrimitiveTimingPlots);


