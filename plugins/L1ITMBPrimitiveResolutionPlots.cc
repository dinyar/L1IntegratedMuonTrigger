/*
 * \file DTLocalLutTriggerTask.cc
 * 
 * $Date: 2011/06/10 13:23:26 $
 * $Revision: 1.1 $
 * \author D. Fasanella - INFN Bologna
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

class L1ITMBPrimitiveResolutionPlots: public edm::EDAnalyzer{
  
public:
  
  /// Constructor
  L1ITMBPrimitiveResolutionPlots(const edm::ParameterSet& ps );
  
  /// Destructor
  virtual ~L1ITMBPrimitiveResolutionPlots();
  
protected:
  
  // BeginJob
  void beginJob();

  ///BeginRun
  void beginRun(const edm::Run& , const edm::EventSetup&);

  /// Find best (highest qual) DCC trigger segments
  void searchDccBest(const std::vector<L1MuDTChambPhDigi>& trigs);
  
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
  int nPhiBins, nPhibBins;
  double rangePhi, rangePhiB;
  
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

};


using namespace edm;
using namespace std;

L1ITMBPrimitiveResolutionPlots::L1ITMBPrimitiveResolutionPlots(const edm::ParameterSet& ps) : trigGeomUtils(0) {
	
  LogTrace("L1ITMBPrimitiveResolutionPlots") << "[L1ITMBPrimitiveResolutionPlots]: Constructor"<<endl;
	
  dccInputTag  = ps.getUntrackedParameter<InputTag>("inputTagDCC");
  segInputTag  = ps.getUntrackedParameter<InputTag>("inputTagSEG");

  overUnderIn      = ps.getUntrackedParameter<bool>("rebinOutFlowsInGraph");
  theGeomLabel     = ps.getUntrackedParameter<string>("geomLabel");
  
  nPhiBins  = 401;
  rangePhi  = 10.025;
  nPhibBins = 401;
  rangePhiB = 10.025;

  parameters = ps;
	
}


L1ITMBPrimitiveResolutionPlots::~L1ITMBPrimitiveResolutionPlots() {
  
  LogTrace("L1ITMBPrimitiveResolutionPlots") << "[L1ITMBPrimitiveResolutionPlots]: analyzed " << nEvents << " events" << endl;
  if (trigGeomUtils) { delete trigGeomUtils; }
	
}


void L1ITMBPrimitiveResolutionPlots::beginJob(){
	
  LogTrace("L1ITMBPrimitiveResolutionPlots") << "[L1ITMBPrimitiveResolutionPlots]: BeginJob" << endl;
  nEvents = 0;
  nLumis  = 0;
	
}

void L1ITMBPrimitiveResolutionPlots::bookHistos(DTChamberId chId) {

  stringstream wheel; wheel << chId.wheel();	
  stringstream sector; sector << chId.sector();
  stringstream station; station << chId.station();	

  TFileDirectory folder  = fs->mkdir(("Wheel" + wheel.str() + "/Sector" + sector.str() +
				      "/Station" + station.str() + "/Segment").c_str());
  
  string chTag = "_W" + wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
  std::map<std::string, TH1F*> &chambMap = chHistos[chId.rawId()];

  string hName = "DCC_PhiResidual";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Trigger local position - Segment local position (correlated triggers)",nPhiBins,-rangePhi,rangePhi); 
  hName = "DCC_PhibResidual";
  chambMap[hName] = folder.make<TH1F>((hName+chTag).c_str(),"Trigger local direction - Segment local direction (correlated triggers)",nPhibBins,-rangePhiB,rangePhiB);

}

  

void L1ITMBPrimitiveResolutionPlots::beginRun(const edm::Run& run, const edm::EventSetup& context) {
	
  LogTrace("L1ITMBPrimitiveResolutionPlots" ) << "[L1ITMBPrimitiveResolutionPlots]: BeginRun" << endl;   
	
  context.get<MuonGeometryRecord>().get(theGeomLabel,muonGeom);
  trigGeomUtils = new DTTrigGeomUtils(muonGeom);
  const std::vector<const DTChamber*> & chambers = muonGeom->chambers();
  
  for ( const auto & chambIt : chambers ) {
    bookHistos( chambIt->id() );
  }
}


void L1ITMBPrimitiveResolutionPlots::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {

  nLumis++;
  LogTrace("L1ITMBPrimitiveResolutionPlots") << "[L1ITMBPrimitiveResolutionPlots]: Begin of LS transition" << endl;
  
}


void L1ITMBPrimitiveResolutionPlots::endJob(){
	
  LogVerbatim("L1ITMBPrimitiveResolutionPlots") << "[L1ITMBPrimitiveResolutionPlots]: analyzed " << nEvents << " events" << endl;
	
}


void L1ITMBPrimitiveResolutionPlots::analyze(const edm::Event& e, const edm::EventSetup& c){
	
  nEvents++;
    
  edm::Handle<L1MuDTChambPhContainer> trigHandle;
  e.getByLabel(dccInputTag,trigHandle);
  const vector<L1MuDTChambPhDigi> & trigs = *(trigHandle->getContainer());
  searchDccBest(trigs);

  Handle<DTRecSegment4DCollection> segments4D;
  e.getByLabel(segInputTag,segments4D);		
  DTRecSegment4DCollection::id_iterator chamberId;

  // Preliminary loop finds best 4D Segment and high quality ones
  vector<const DTRecSegment4D*> best4DSegments;
  for ( chamberId = segments4D->id_begin(); chamberId != segments4D->id_end();
	++chamberId ) {
		
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
  
  for ( const auto & bestTrackIt : best4DSegments ) {
    
    if ( bestTrackIt->hasPhi() ) {
      
      DTChamberId chId = bestTrackIt->chamberId();
      int nHitsPhi = bestTrackIt->phiSegment()->degreesOfFreedom()+2;
      
      int wheel    = chId.wheel();
      int station  = chId.station();
      int scsector = 0;
      float trackPosPhi, trackPosEta, trackDirPhi, trackDirEta;
      trigGeomUtils->computeSCCoordinates(bestTrackIt,scsector,trackPosPhi,trackDirPhi,trackPosEta,trackDirEta);
      
      map<string, TH1F*> &chMap = chHistos[chId.rawId()];

      if (trigQualBest[wheel+3][station][scsector] > -1 &&  // residuals only for correlate triggers
	  trigQualBest[wheel+3][station][scsector] < 7  &&
	  nHitsPhi>=7 ) {
					
	float trigPos = trigGeomUtils->trigPos(trigBest[wheel+3][station][scsector]);
	float trigDir = trigGeomUtils->trigDir(trigBest[wheel+3][station][scsector]);
	trigGeomUtils->trigToSeg(station,trigPos,trackDirPhi);
										
	double deltaPos = trigPos-trackPosPhi;
	deltaPos = overUnderIn ? max(min(deltaPos,rangePhi-0.01),-rangePhi+0.01) : deltaPos;
	double deltaDir = trigDir-trackDirPhi;
	deltaDir = overUnderIn ? max(min(deltaDir,rangePhiB-0.01),-rangePhiB+0.01) : deltaDir;
	chMap.find("DCC_PhiResidual")->second->Fill(deltaPos);
	chMap.find("DCC_PhibResidual")->second->Fill(deltaDir);

      }
      
    }
  } 
  
}

void L1ITMBPrimitiveResolutionPlots::searchDccBest( const std::vector<L1MuDTChambPhDigi> & trigs ){
  
  string histoType ;
  string histoTag ;
  
  // define best quality trigger segment
  // start from 1 and zero is kept empty
  for (int st=0;st<=4;++st)
    for (int wh=0;wh<=5;++wh)
      for (int sec=0;sec<=12;++sec)
	trigQualBest[wh][st][sec] = -1;    
	
  for ( const auto & trigIt : trigs ) {
    
    int wh   = trigIt.whNum();
    int sec  = trigIt.scNum() + 1; // DTTF -> DT sector range transform
    int st   = trigIt.stNum();
    int qual = trigIt.code();

    if(qual>trigQualBest[wh+3][st][sec] && qual<7) {
      trigQualBest[wh+3][st][sec]=qual; 
      trigBest[wh+3][st][sec] = &(trigIt);
    }
    
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1ITMBPrimitiveResolutionPlots);


