#include "L1Trigger/L1IntegratedMuonTrigger/interface/GeometryTranslator.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"

// event setup stuff / geometries
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "L1Trigger/DTUtilities/interface/DTTrigGeom.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include <cmath> // for pi

using namespace L1ITMu;

GeometryTranslator::GeometryTranslator():  
  _geom_cache_id(0ULL) {
}

GeometryTranslator::~GeometryTranslator() {  
}

double 
GeometryTranslator::calculateGlobalEta(const TriggerPrimitive& tp) const {
  switch(tp.subsystem()) {
  case TriggerPrimitive::kDT:
    return calcDTSpecificEta(tp);
    break;
  case TriggerPrimitive::kCSC:
    return calcCSCSpecificEta(tp);
    break;
  case TriggerPrimitive::kRPC:
    return calcRPCSpecificEta(tp);
    break;
  default:
    return std::nan("Invalid TP type!"); 
    break;
  }
}

double 
GeometryTranslator::calculateGlobalPhi(const TriggerPrimitive& tp) const {
  switch(tp.subsystem()) {
  case TriggerPrimitive::kDT:
    return calcDTSpecificPhi(tp);
    break;
  case TriggerPrimitive::kCSC:
    return calcCSCSpecificPhi(tp);
    break;
  case TriggerPrimitive::kRPC:
    return calcRPCSpecificPhi(tp);
    break;
  default:
    return std::nan("Invalid TP type!");
    break;
  }
}

double 
GeometryTranslator::calculateBendAngle(const TriggerPrimitive& tp) const {
  switch(tp.subsystem()) {
  case TriggerPrimitive::kDT:
    return calcDTSpecificBend(tp);
    break;
  case TriggerPrimitive::kCSC:
    return calcCSCSpecificBend(tp);
    break;
  case TriggerPrimitive::kRPC:
    return calcRPCSpecificBend(tp);
    break;
  default:
    return std::nan("Invalid TP type!");
    break;
  }
}

double 
GeometryTranslator::calculateGlobalRho(const TriggerPrimitive& tp) const {
  switch(tp.subsystem()) {
  case TriggerPrimitive::kDT:
    return calcDTSpecificRho(tp);
    break;
  case TriggerPrimitive::kCSC:
    return calcCSCSpecificRho(tp);
    break;
  case TriggerPrimitive::kRPC:
    return calcRPCSpecificRho(tp);
    break;
  default:
    return std::nan("Invalid TP type!");
    break;
  }
}

void GeometryTranslator::checkAndUpdateGeometry(const edm::EventSetup& es) {
  const MuonGeometryRecord& geom = es.get<MuonGeometryRecord>();
  unsigned long long geomid = geom.cacheIdentifier();
  if( _geom_cache_id != geomid ) {
    geom.get(_georpc);  
    geom.get(_geocsc);    
    geom.get(_geodt);
    _geom_cache_id = geomid;
  }  
}

GlobalPoint 
GeometryTranslator::getRPCSpecificPoint(const TriggerPrimitive& tp) const {
  const RPCDetId id(tp.detId<RPCDetId>());
  std::unique_ptr<const RPCRoll>  roll(_georpc->roll(id));
  const uint16_t strip = tp.getRPCData().strip;
  const LocalPoint lp = roll->centreOfStrip(strip);
  const GlobalPoint gp = roll->toGlobal(lp);
  
  roll.release();
 
  return gp;
}

double 
GeometryTranslator::calcRPCSpecificEta(const TriggerPrimitive& tp) const {  
  return getRPCSpecificPoint(tp).eta();
}

double 
GeometryTranslator::calcRPCSpecificPhi(const TriggerPrimitive& tp) const {  
  return getRPCSpecificPoint(tp).phi();
}

double 
GeometryTranslator::calcRPCSpecificRho(const TriggerPrimitive& tp) const {  
  return getRPCSpecificPoint(tp).perp();
}

// this function actually does nothing since RPC
// hits are point-like objects
double 
GeometryTranslator::calcRPCSpecificBend(const TriggerPrimitive& tp) const {
  return 0.0;
}


// alot of this is transcription and consolidation of the CSC
// global phi calculation code
// this works directly with the geometry 
// rather than using the old phi luts
GlobalPoint 
GeometryTranslator::getCSCSpecificPoint(const TriggerPrimitive& tp) const {
  const CSCDetId id(tp.detId<CSCDetId>()); 
  // we should change this to weak_ptrs at some point
  // requires introducing std::shared_ptrs to geometry
  std::unique_ptr<const CSCChamber> chamb(_geocsc->chamber(id));
  std::unique_ptr<const CSCLayerGeometry> layer_geom(
    chamb->layer(CSCConstants::KEY_ALCT_LAYER)->geometry()
    );
  std::unique_ptr<const CSCLayer> layer(
    chamb->layer(CSCConstants::KEY_ALCT_LAYER)
    );
  
  const uint16_t halfstrip = tp.getCSCData().strip;
  const uint16_t pattern = tp.getCSCData().pattern;
  const uint16_t keyWG = tp.getCSCData().keywire; 
  //const unsigned maxStrips = layer_geom->numberOfStrips();  

  // so we can extend this later 
  // assume TMB2007 half-strips only as baseline
  double offset = 0.0;
  switch(1) {
  case 1:
    offset = CSCPatternLUT::get2007Position(pattern);
  }
  const unsigned halfstrip_offs = unsigned(0.5 + halfstrip + offset);
  const unsigned strip = halfstrip_offs/2 + 1; // geom starts from 1

  // the rough location of the hit at the ALCT key layer
  // we will refine this using the half strip information
  const LocalPoint coarse_lp = 
    layer_geom->stripWireGroupIntersection(strip,keyWG);  
  const GlobalPoint coarse_gp = layer->surface().toGlobal(coarse_lp);  
  
  // the strip width/4.0 gives the offset of the half-strip
  // center with respect to the strip center
  const double hs_offset = layer_geom->stripPhiPitch()/4.0;
  
  // determine handedness of the chamber
  const bool ccw = isCSCCounterClockwise(layer);
  // we need to subtract the offset of even half strips and add the odd ones
  const double phi_offset = ( ( halfstrip_offs%2 ? 1 : -1)*
			      ( ccw ? -hs_offset : hs_offset ) );
  
  // the global eta calculation uses the middle of the strip
  // so no need to increment it
  const GlobalPoint final_gp( GlobalPoint::Polar( coarse_gp.theta(),
						  (coarse_gp.phi().value() + 
						   phi_offset),
						  coarse_gp.mag() ) );
    
  // We need to add in some notion of the 'error' on trigger primitives
  // like the width of the wire group by the width of the strip
  // or something similar      

  // release ownership of the pointers
  chamb.release();
  layer_geom.release();
  layer.release();
  
  return final_gp;
}

double 
GeometryTranslator::calcCSCSpecificEta(const TriggerPrimitive& tp) const {  
  return getCSCSpecificPoint(tp).eta();
}

double 
GeometryTranslator::calcCSCSpecificPhi(const TriggerPrimitive& tp) const {  
  return getCSCSpecificPoint(tp).phi();
}

double 
GeometryTranslator::calcCSCSpecificRho(const TriggerPrimitive& tp) const {  
  return getCSCSpecificPoint(tp).perp();
}

double 
GeometryTranslator::calcCSCSpecificBend(const TriggerPrimitive& tp) const {  
  return 0.0;
}

GlobalPoint
GeometryTranslator::calcDTSpecificPoint(const TriggerPrimitive& tp) const {  
  const DTChamberId baseid(tp.detId<DTChamberId>());
  // do not use this pointer for anything other than creating a trig geom
  std::unique_ptr<DTChamber> chamb(
    const_cast<DTChamber*>(_geodt->chamber(baseid)) 
    );
  std::unique_ptr<DTTrigGeom> trig_geom( new DTTrigGeom(chamb.get(),false) );
  chamb.release(); // release it here so no one gets funny ideas  
  // super layer one is the theta superlayer in a DT chamber
  // station 4 does not have a theta super layer
  // the BTI index from the theta trigger is an OR of some BTI outputs
  // so, we choose the BTI that's in the middle of the group
  // as the BTI that we get theta from
  // TODO:::::>>> need to make sure this ordering doesn't flip under wheel sign
  const int NBTI_theta = ( (baseid.station() != 4) ? 
			   trig_geom->nCell(2) : trig_geom->nCell(3) );
  const int bti_group = tp.getDTData().theta_bti_group;
  const unsigned bti_actual = bti_group*NBTI_theta/7 + NBTI_theta/14 + 1;  
  DTBtiId thetaBTI;  
  if ( baseid.station() != 4 && bti_group != -1) {
    thetaBTI = DTBtiId(baseid,2,bti_actual);
  } else {
    // since this is phi oriented it'll give us theta in the middle
    // of the chamber
    thetaBTI = DTBtiId(baseid,3,1); 
  }
  const GlobalPoint theta_gp = trig_geom->CMSPosition(thetaBTI);
  
  // local phi in sector -> global phi
  double phi_in_sector = ((double)tp.getDTData().radialAngle)/4096.0; 
  double sector_offset = tp.getDTData().sector*M_PI/6.0;
  double phi = phi_in_sector + sector_offset; // add sector offset  

  int qual = tp.getDTData().qualityCode;
  double z_chamber = (qual == 1 || qual == 3) ? trig_geom->ZSL(3) :
                     (qual == 0 || qual == 2) ? trig_geom->ZSL(1) :
                     trig_geom->ZcenterSL();

  GlobalPoint ch_coord = trig_geom->toGlobal( LocalPoint(0,0,z_chamber) );

  double rho_p = ch_coord.perp()/cos( acos(cos(ch_coord.phi()-sector_offset)) );
  double rho = rho_p/cos(phi_in_sector);
  
  return GlobalPoint( GlobalPoint::Cylindrical( rho, phi, theta_gp.z() ) );
}

double 
GeometryTranslator::calcDTSpecificEta(const TriggerPrimitive& tp) const {  
  return calcDTSpecificPoint(tp).eta();
}

double 
GeometryTranslator::calcDTSpecificPhi(const TriggerPrimitive& tp) const {
  return calcDTSpecificPoint(tp).phi();
}

double 
GeometryTranslator::calcDTSpecificRho(const TriggerPrimitive& tp) const {
  return calcDTSpecificPoint(tp).perp();
}

// we have the bend except for station 3
double 
GeometryTranslator::calcDTSpecificBend(const TriggerPrimitive& tp) const {
  int bend = tp.getDTData().bendingAngle;
  double bendf = bend/512.0;
  return bendf;
}

bool GeometryTranslator::
isCSCCounterClockwise(const std::unique_ptr<const CSCLayer>& layer) const {
  const int nStrips = layer->geometry()->numberOfStrips();
  const double phi1 = layer->centerOfStrip(1).phi();
  const double phiN = layer->centerOfStrip(nStrips).phi();
  return ( (std::abs(phi1 - phiN) < M_PI  && phi1 >= phiN) || 
	   (std::abs(phi1 - phiN) >= M_PI && phi1 < phiN)     );  
}
