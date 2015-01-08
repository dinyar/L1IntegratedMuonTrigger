/*
 * EDAnalyzer L1ITMBPrimitiveResolutionComparePlots
 *
 * This analyzer is based on the L1ITMBPrimitiveResolutionPlots.
 * It creates the same plots but for two DCC input tags at once (e.g.
 * the old and the new trigger algorithms). It also creates the related
 * scatterplots.
 *
 * The code is a heavily modified version of the original class/config.
 * In addition to implementing the new features it has also been
 * cleaned up.
 *
 * The main differences are:
 *  - the analyzer now takes two DCC input tags (inputTagDCC_old, inputTagDCC_new)
 *  - it does not generate the unnecessary "Segment" directory
 */
#include "L1Trigger/L1IntegratedMuonTrigger/interface/L1DTTrigGeomUtils.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TH2.h>

#include <sstream>
#include <string>
#include <map>

//const int nquals = 16;
const int nquals = 7;

class L1ITMBPrimitiveResolutionComparePlots : public edm::EDAnalyzer {
public:
  explicit L1ITMBPrimitiveResolutionComparePlots(const edm::ParameterSet& ps);
  virtual ~L1ITMBPrimitiveResolutionComparePlots();

protected:
  void beginJob();
  void beginRun(const edm::Run& , const edm::EventSetup&);
  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, const edm::EventSetup& context);
  void endJob(void);
  void analyze(const edm::Event& e, const edm::EventSetup& c);

private:
  int nEvents, nLumis;
  int nPhiBins, nPhibBins;
  double rangePhi, rangePhiB;

  bool overUnderIn;
  edm::InputTag segInputTag;
  edm::Service<TFileService> fs;

  std::string theGeomLabel;
  L1DTTrigGeomUtils* trigGeomUtils;
  edm::ESHandle<DTGeometry> muonGeom;
  std::map<uint32_t, std::map<std::string, TH2F*> > scatterHistos;
  std::map<uint32_t, std::map<std::string, TH2F*> > scatterHistosQs[nquals];

  struct t_vars {
    edm::InputTag dccInputTag;
    std::string identifier;
    int trigQualBest[6][5][13];
    const L1MuDTChambPhDigi* trigBest[6][5][13];
    std::map<uint32_t, std::map<std::string, TH1F*>> chHistos;
    std::map<uint32_t, std::map<std::string, TH1F*>> chHistosQs[nquals];
  } oldVars, newVars;

  struct t_deltas {
    float deltaPos, deltaDir;
    int quality;
  };

  void bookHistos(const std::string &hName, const DTChamberId &chId, TFileDirectory &folder);
  bool analyze_tag( t_vars &vars, const DTRecSegment4D* bestTrack, t_deltas &deltas);
  void searchDccBest( t_vars &vars, const std::vector<L1MuDTChambPhDigi> &trigs);
};

using namespace std;
using namespace edm;

L1ITMBPrimitiveResolutionComparePlots::L1ITMBPrimitiveResolutionComparePlots(const edm::ParameterSet& ps) : trigGeomUtils(0) {
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: Constructor" << endl;

  oldVars.dccInputTag = ps.getUntrackedParameter<InputTag>("inputTagDCC_old");
  newVars.dccInputTag = ps.getUntrackedParameter<InputTag>("inputTagDCC_new");
  segInputTag = ps.getUntrackedParameter<InputTag>("inputTagSEG");
  overUnderIn = ps.getUntrackedParameter<bool>("rebinOutFlowsInGraph");
  theGeomLabel = ps.getUntrackedParameter<string>("geomLabel");

  nPhiBins  = 401;
  rangePhi  = 10.025;
  nPhibBins = 401;
  rangePhiB = 10.025;
}

L1ITMBPrimitiveResolutionComparePlots::~L1ITMBPrimitiveResolutionComparePlots() {
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: analyzed " << nEvents << " events" << endl;
  if(trigGeomUtils) delete trigGeomUtils;
}

void L1ITMBPrimitiveResolutionComparePlots::beginJob(){
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: BeginJob" << endl;
  nEvents = 0;
  nLumis  = 0;
}

void L1ITMBPrimitiveResolutionComparePlots::beginRun(const edm::Run& run, const edm::EventSetup& context) {
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: BeginRun" << endl;

  oldVars.identifier = "old";
  newVars.identifier = "new";

  context.get<MuonGeometryRecord>().get(theGeomLabel,muonGeom);
  trigGeomUtils = new L1DTTrigGeomUtils(muonGeom);

  for(const DTChamber * chamber : muonGeom->chambers()) {
    DTChamberId chId = chamber->id();

    ostringstream tdir_path;
    tdir_path << "Wheel" << chId.wheel() << "/Sector" << chId.sector() << "/Station" << chId.station();
    TFileDirectory folder = fs->mkdir(tdir_path.str());

    bookHistos("DCC_PhiResidual", chId, folder);
    bookHistos("DCC_PhibResidual", chId, folder);
  }
}

void L1ITMBPrimitiveResolutionComparePlots::bookHistos(const string &hName, const DTChamberId &chId, TFileDirectory &folder) {
  TFileDirectory subfolder = folder.mkdir(hName);

  ostringstream chTag_oss;
  chTag_oss << "_W" << chId.wheel() << "_Sec" << chId.sector() << "_St" << chId.station();
  string chTag = chTag_oss.str();

  this->newVars.chHistos[chId.rawId()][hName] = subfolder.make<TH1F>(
								     (hName+"_"+newVars.identifier+chTag).c_str(),
								     "Trigger local position - Segment local position (correlated triggers)",
								     nPhiBins,-rangePhi,rangePhi
								     );

  this->oldVars.chHistos[chId.rawId()][hName] = subfolder.make<TH1F>(
								     (hName+"_"+oldVars.identifier+chTag).c_str(),
								     "Trigger local position - Segment local position (correlated triggers)",
								     nPhiBins,-rangePhi,rangePhi
								     );

  this->scatterHistos[chId.rawId()][hName+"_scatter"] = subfolder.make<TH2F>(
									     (hName+"_scatter"+chTag).c_str(),
									     "Trigger local position - Segment local position (correlated triggers) - scatter plot new vs. old",
									     nPhiBins,-rangePhi,rangePhi,nPhiBins,-rangePhi,rangePhi
									     );

  for(int i=0; i<nquals; i++) {
    ostringstream qdir_oss; qdir_oss << "quality_" << i;
    TFileDirectory qdir = subfolder.mkdir(qdir_oss.str());

    string chTagQ = chTag + "_" + qdir_oss.str();
    this->newVars.chHistosQs[i][chId.rawId()][hName] = qdir.make<TH1F>(
								       (hName+"_"+newVars.identifier+chTagQ).c_str(),
								       "Trigger local position - Segment local position (correlated triggers)",
								       nPhiBins,-rangePhi,rangePhi
								       );

    this->oldVars.chHistosQs[i][chId.rawId()][hName] = qdir.make<TH1F>(
								       (hName+"_"+oldVars.identifier+chTagQ).c_str(),
								       "Trigger local position - Segment local position (correlated triggers)",
								       nPhiBins,-rangePhi,rangePhi
								       );

    this->scatterHistosQs[i][chId.rawId()][hName+"_scatter"] = qdir.make<TH2F>(
									       (hName+"_scatter"+chTagQ).c_str(),
									       "Trigger local position - Segment local position (correlated triggers) - scatter plot new vs. old",
									       nPhiBins,-rangePhi,rangePhi,nPhiBins,-rangePhi,rangePhi
									       );
  }
}

/*
 * To reset the MEs
 */
void L1ITMBPrimitiveResolutionComparePlots::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {
  nLumis++;
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: Begin of LS transition" << endl;
}


void L1ITMBPrimitiveResolutionComparePlots::endJob(){
  LogVerbatim("L1ITMBPrimitiveResolutionComparePlots") << "[L1ITMBPrimitiveResolutionComparePlots]: analyzed " << nEvents << " events" << endl;
}

void L1ITMBPrimitiveResolutionComparePlots::analyze(const edm::Event& e, const edm::EventSetup& c) {
  nEvents++;
  LogTrace("L1ITMBPrimitiveResolutionComparePlots") << "Analyzing event" << nEvents << endl;

  edm::Handle<L1MuDTChambPhContainer> trigHandle;

  e.getByLabel(oldVars.dccInputTag, trigHandle);
  searchDccBest(oldVars, *trigHandle->getContainer());

  e.getByLabel(newVars.dccInputTag, trigHandle);
  searchDccBest(newVars, *trigHandle->getContainer());

  Handle<DTRecSegment4DCollection> segments4D;
  e.getByLabel(segInputTag,segments4D);		
  DTRecSegment4DCollection::id_iterator chamberId;

  // Preliminary loop finds best 4D Segment and high quality ones
  vector<const DTRecSegment4D*> best4DSegments;
  for ( chamberId = segments4D->id_begin(); chamberId != segments4D->id_end();
	++chamberId ) {

    const DTRecSegment4D* tmpBest = 0;
    int tmpdof = 0;
    int dof = 0;

    DTRecSegment4DCollection::range rangeInCh = segments4D->get(*chamberId);
    DTRecSegment4DCollection::const_iterator trackIt  = rangeInCh.first;
    DTRecSegment4DCollection::const_iterator trackEnd = rangeInCh.second;

    for (; trackIt!=trackEnd; ++trackIt){

      if(trackIt->hasPhi()) {
	dof = trackIt->phiSegment()->degreesOfFreedom();
	if(dof>tmpdof) {
	  tmpBest = &(*trackIt);
	  tmpdof = dof;
	}
      }
    }

    if(tmpBest) best4DSegments.push_back(tmpBest);
  }

  for(const DTRecSegment4D* bestTrack : best4DSegments) {
    t_deltas deltas_new, deltas_old;
    if(analyze_tag(newVars, bestTrack, deltas_new) && analyze_tag(oldVars, bestTrack, deltas_old)) {
      map<string, TH2F*> &scatterMap = scatterHistos[bestTrack->chamberId().rawId()];
      scatterMap.find("DCC_PhiResidual_scatter")->second->Fill(deltas_old.deltaPos, deltas_new.deltaPos);
      scatterMap.find("DCC_PhibResidual_scatter")->second->Fill(deltas_old.deltaDir, deltas_new.deltaDir);

      map<string, TH2F*> &scatterMapQ = scatterHistosQs[deltas_new.quality][bestTrack->chamberId().rawId()];
      scatterMapQ.find("DCC_PhiResidual_scatter")->second->Fill(deltas_old.deltaPos, deltas_new.deltaPos);
      scatterMapQ.find("DCC_PhibResidual_scatter")->second->Fill(deltas_old.deltaDir, deltas_new.deltaDir);
    }
  }
}

bool L1ITMBPrimitiveResolutionComparePlots::analyze_tag( t_vars &vars, const DTRecSegment4D* bestTrack, t_deltas &deltas) {
  if(bestTrack->hasPhi()) {
    DTChamberId chId = bestTrack->chamberId();
    int nHitsPhi = bestTrack->phiSegment()->degreesOfFreedom()+2;

    int wheel    = chId.wheel();
    int station  = chId.station();
    int scsector = 0;
    float trackPosPhi, trackPosEta, trackDirPhi, trackDirEta;
    trigGeomUtils->computeSCCoordinates(bestTrack,scsector,trackPosPhi,trackDirPhi,trackPosEta,trackDirEta);

    int best_qualitycode = vars.trigQualBest[wheel+3][station][scsector];
    const L1MuDTChambPhDigi* best_trig = vars.trigBest[wheel+3][station][scsector];

    if(best_qualitycode > -1 && best_qualitycode < nquals && nHitsPhi >= 7) {
      float trigPos = trigGeomUtils->trigPos(best_trig);
      float trigDir = trigGeomUtils->trigDir(best_trig);
      trigGeomUtils->trigToSeg(station,trigPos,trackDirPhi);

      double deltaPos = trigPos-trackPosPhi;
      deltaPos = overUnderIn ? max(min(deltaPos,rangePhi-0.01),-rangePhi+0.01) : deltaPos;
      double deltaDir = trigDir-trackDirPhi;
      deltaDir = overUnderIn ? max(min(deltaDir,rangePhiB-0.01),-rangePhiB+0.01) : deltaDir;

      map<string, TH1F*> &chMap = vars.chHistos[chId.rawId()];
      chMap.find("DCC_PhiResidual")->second->Fill(deltaPos);
      chMap.find("DCC_PhibResidual")->second->Fill(deltaDir);

      map<string, TH1F*> &chMapQ = vars.chHistosQs[best_qualitycode][chId.rawId()];
      chMapQ.find("DCC_PhiResidual")->second->Fill(deltaPos);
      chMapQ.find("DCC_PhibResidual")->second->Fill(deltaDir);

      deltas.deltaPos = deltaPos;
      deltas.deltaDir = deltaDir;
      deltas.quality = best_qualitycode;
      return true;
    }
  }
  return false;
}

/*
 * Finds the best (i.e. highest quality) DCC trigger segments.
 */
void L1ITMBPrimitiveResolutionComparePlots::searchDccBest( t_vars &vars, const std::vector<L1MuDTChambPhDigi> &trigs) {
  // define best quality trigger segment
  // start from 1 and zero is kept empty
  for (int st=0;st<=4;++st) {
    for (int wh=0;wh<=5;++wh) {
      for (int sec=0;sec<=12;++sec) {
	vars.trigQualBest[wh][st][sec] = -1;
      }
    }
  }

  for ( const L1MuDTChambPhDigi& trig : trigs ) {
    int whId  = trig.whNum() + 3;
    int secId = trig.scNum() + 1; // DTTF -> DT sector range transform
    int stId  = trig.stNum();
    int qualitycode = trig.code();

    if(qualitycode > vars.trigQualBest[whId][stId][secId] && qualitycode < nquals) {
      vars.trigQualBest[whId][stId][secId] = qualitycode;
      vars.trigBest[whId][stId][secId] = &trig;
    }
  }
}

DEFINE_FWK_MODULE(L1ITMBPrimitiveResolutionComparePlots);
