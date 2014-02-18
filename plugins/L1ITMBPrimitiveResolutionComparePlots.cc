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
#include "DQM/DTMonitorModule/interface/DTTrigGeomUtils.h"

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
		DTTrigGeomUtils* trigGeomUtils;
		edm::ESHandle<DTGeometry> muonGeom;
		std::map<uint32_t, std::map<std::string, TH2F*> > scatterHistos;

		struct t_vars {
			edm::InputTag dccInputTag;
			std::string identifier;
			int trigQualBest[6][5][13];
			const L1MuDTChambPhDigi* trigBest[6][5][13];
			std::map<uint32_t, std::map<std::string, TH1F*>> chHistos;
		} oldVars, newVars;

		struct t_deltas {
			float deltaPos, deltaDir;
		};

		void bookHistos(t_vars &vars, TFileDirectory &folder, const std::string chTag, const DTChamberId &chId);
		bool analyze_tag(t_vars &vars, std::vector<const DTRecSegment4D*>::const_iterator &bestTrackIt, t_deltas &deltas);
		void searchDccBest(t_vars &vars, std::vector<L1MuDTChambPhDigi>* trigs);
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
	trigGeomUtils = new DTTrigGeomUtils(muonGeom);

	std::vector<DTChamber*>::const_iterator chambIt  = muonGeom->chambers().begin();
	std::vector<DTChamber*>::const_iterator chambEnd = muonGeom->chambers().end();

	for(; chambIt!=chambEnd; ++chambIt) {
		DTChamberId chId = (*chambIt)->id();

		ostringstream tdir_path, chTag;
		//tdir_path << "Wheel" << chId.wheel() << "/Sector" << chId.sector() << "/Station" << chId.station() << "/Segment";
		tdir_path << "Wheel" << chId.wheel() << "/Sector" << chId.sector() << "/Station" << chId.station();
		chTag << "_W" << chId.wheel() << "_Sec" << chId.sector() << "_St" << chId.station();

		TFileDirectory folder = fs->mkdir(tdir_path.str().c_str());

		bookHistos(newVars, folder, chTag.str(), chId);
		bookHistos(oldVars, folder, chTag.str(), chId);

		string hName;
		map<std::string, TH2F*> &chambMap2D = scatterHistos[chId.rawId()];

		hName = "DCC_PhiResidual_scatter";
		chambMap2D[hName] = folder.make<TH2F>(
			(hName+chTag.str()).c_str(),
			"Trigger local position - Segment local position (correlated triggers) - scatter plot new vs. old",
			nPhiBins,-rangePhi,rangePhi,nPhiBins,-rangePhi,rangePhi
		);

		hName = "DCC_PhibResidual_scatter";
		chambMap2D[hName] = folder.make<TH2F>(
			(hName+chTag.str()).c_str(),
			"Trigger local direction - Segment local direction (correlated triggers) - scatter plot new vs. old",
			nPhibBins,-rangePhiB,rangePhiB,nPhibBins,-rangePhiB,rangePhiB
		);
	}
}

void L1ITMBPrimitiveResolutionComparePlots::bookHistos(t_vars &vars, TFileDirectory &folder, const string chTag, const DTChamberId &chId) {
	string hName;
	map<std::string, TH1F*> &chambMap = vars.chHistos[chId.rawId()];

	hName = "DCC_PhiResidual";
	chambMap[hName] = folder.make<TH1F>(
		(hName+"_"+vars.identifier+chTag).c_str(),
		"Trigger local position - Segment local position (correlated triggers)",
		nPhiBins,-rangePhi,rangePhi\
	);

	hName = "DCC_PhibResidual";
	chambMap[hName] = folder.make<TH1F>(
		(hName+"_"+vars.identifier+chTag).c_str(),
		"Trigger local direction - Segment local direction (correlated triggers)",
		nPhibBins,-rangePhiB,rangePhiB
	);
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
	searchDccBest(oldVars, trigHandle->getContainer());

	e.getByLabel(newVars.dccInputTag, trigHandle);
	searchDccBest(newVars, trigHandle->getContainer());

	Handle<DTRecSegment4DCollection> segments4D;
	e.getByLabel(segInputTag,segments4D);

	// Preliminary loop finds best 4D Segment and high quality ones
	vector<const DTRecSegment4D*> best4DSegments;
	DTRecSegment4DCollection::id_iterator chamberId;
	for(chamberId = segments4D->id_begin(); chamberId != segments4D->id_end(); ++chamberId) {
		DTRecSegment4DCollection::range  rangeInCh = segments4D->get(*chamberId);
		DTRecSegment4DCollection::const_iterator trackIt  = rangeInCh.first;
		DTRecSegment4DCollection::const_iterator trackEnd = rangeInCh.second;

		const DTRecSegment4D* tmpBest = 0;
		int tmpdof = 0, dof = 0;

		for(; trackIt!=trackEnd; ++trackIt) {
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

	vector<const DTRecSegment4D*>::const_iterator bestTrackIt  = best4DSegments.begin();
	vector<const DTRecSegment4D*>::const_iterator bestTrackEnd = best4DSegments.end();
	for(; bestTrackIt!=bestTrackEnd; ++bestTrackIt) {
		t_deltas deltas_new, deltas_old;
		map<string, TH2F*> &scatterMap = scatterHistos[(*bestTrackIt)->chamberId().rawId()];
		if(analyze_tag(newVars, bestTrackIt, deltas_new) && analyze_tag(oldVars, bestTrackIt, deltas_old)) {
			scatterMap.find("DCC_PhiResidual_scatter")->second->Fill(deltas_old.deltaPos, deltas_new.deltaPos);
			scatterMap.find("DCC_PhibResidual_scatter")->second->Fill(deltas_old.deltaDir, deltas_new.deltaDir);
		}
	}
}

bool L1ITMBPrimitiveResolutionComparePlots::analyze_tag(t_vars &vars, vector<const DTRecSegment4D*>::const_iterator &bestTrackIt, t_deltas &deltas) {
	if((*bestTrackIt)->hasPhi()) {
		DTChamberId chId = (*bestTrackIt)->chamberId();
		int nHitsPhi = (*bestTrackIt)->phiSegment()->degreesOfFreedom()+2;

		int wheel    = chId.wheel();
		int station  = chId.station();
		int scsector = 0;
		float trackPosPhi, trackPosEta, trackDirPhi, trackDirEta;
		trigGeomUtils->computeSCCoordinates((*bestTrackIt),scsector,trackPosPhi,trackDirPhi,trackPosEta,trackDirEta);

		map<string, TH1F*> &chMap = vars.chHistos[chId.rawId()];

		if(vars.trigQualBest[wheel+3][station][scsector] > -1 // residuals only for correlate triggers
		  && vars.trigQualBest[wheel+3][station][scsector] < 7
		  && nHitsPhi>=7) {
			float trigPos = trigGeomUtils->trigPos(vars.trigBest[wheel+3][station][scsector]);
			float trigDir = trigGeomUtils->trigDir(vars.trigBest[wheel+3][station][scsector]);
			trigGeomUtils->trigToSeg(station,trigPos,trackDirPhi);

			double deltaPos = trigPos-trackPosPhi;
			deltaPos = overUnderIn ? max(min(deltaPos,rangePhi-0.01),-rangePhi+0.01) : deltaPos;
			double deltaDir = trigDir-trackDirPhi;
			deltaDir = overUnderIn ? max(min(deltaDir,rangePhiB-0.01),-rangePhiB+0.01) : deltaDir;

			chMap.find("DCC_PhiResidual")->second->Fill(deltaPos);
			chMap.find("DCC_PhibResidual")->second->Fill(deltaDir);

			deltas.deltaPos = deltaPos;
			deltas.deltaDir = deltaDir;
			return true;
		}
	}
	return false;
}

/*
 * Finds the best (i.e. highest quality) DCC trigger segments.
 */
void L1ITMBPrimitiveResolutionComparePlots::searchDccBest(t_vars &vars, std::vector<L1MuDTChambPhDigi>* trigs) {
	string histoType, histoTag;

	// define best quality trigger segment
	// start from 1 and zero is kept empty
	for (int st=0;st<=4;++st) {
		for (int wh=0;wh<=5;++wh) {
			for (int sec=0;sec<=12;++sec) {
				vars.trigQualBest[wh][st][sec] = -1;
			}
		}
	}

	vector<L1MuDTChambPhDigi>::const_iterator trigIt  = trigs->begin();
	vector<L1MuDTChambPhDigi>::const_iterator trigEnd = trigs->end();
	for(; trigIt!=trigEnd; ++trigIt) {
		int whId  = trigIt->whNum() + 3;
		int secId = trigIt->scNum() + 1; // DTTF -> DT sector range transform
		int stId  = trigIt->stNum();
		int qualitycode = trigIt->code();

		if(qualitycode > vars.trigQualBest[whId][stId][secId] && qualitycode < 7) {
			vars.trigQualBest[whId][stId][secId] = qualitycode;
			vars.trigBest[whId][stId][secId] = &(*trigIt);
		}
	}
}

DEFINE_FWK_MODULE(L1ITMBPrimitiveResolutionComparePlots);
