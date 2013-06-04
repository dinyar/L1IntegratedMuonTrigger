

// #include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

// #include "DataFormats/MuonDetId/interface/CSCDetId.h"
// #include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBhelpers.h"

namespace {
  // DT TF relative segment address explanation
  //
  // a code of '15' (stations 2,3,4) or '3' (station 1 only)
  // means there was no valid track extrapolation in this DT station
  //
  // schematic diagram of DT codes with corresponding VHDL addresses in ():
  // the phi direction seems to be the direction with respect to the
  // track's sector processor. 
  // ( this is why station one only has addresses 1,2 )
  //         --------------------------------------
  //         |   4 (10)  5 (11) |   6 (2)  7 (3)  |   ( next sector )
  //      P  ------------+------------------------- 
  //      H  |   0 (8)   1 (9)  |   2 (0)  3 (1)  |   ( this sector )
  //      I  ------------+-------------------------
  //         |   8 (12)  9 (13) |  10 (4) 11 (5)  |   ( prev sector )
  //         ------------+-------------------------
  //               this Wheel       next Wheel
 
  bool isExtrapAcrossWheel(const int addr, const int station) {     
    if( station != 1 ) {
      switch(addr) {
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
      case 13:
	return false;
	break;
      default:
	return true;
      }
    } else {
      return !((bool)addr);
    }
    return false;
  } 
  
  int relativeSector(const int addr, const int station) {
    if( station != 1 ){
      switch(addr) {
      case 12:
      case 13:
      case 4:
      case 5:
	return -1;
	break;
      case 8:
      case 9:
      case 0:
      case 1:
	return 0;
	break;
      case 10:
      case 11:
      case 2:
      case 3:
	return 1;
	break;
      default:
	break;
      }
    }
    return 0;
  }
  
  

}

namespace L1ITMu {
  namespace MBhelpers {
    
    MBLTContainer 
    getPrimitivesByMBTriggerInfo(const int wheel,
				 const int sp_wheel,
				 const int sector,
                                const edm::Handle<MBLTContainer>& mbs,
				 const unsigned mode,
				 const std::vector<unsigned>& trkNmbs) {
      MBLTContainer result;
      auto mb = mbs->cbegin();
//       auto mbbeg = mbs->cbegin();
      auto mbend = mbs->cend();
      
      std::vector<unsigned>::const_iterator ista;
      auto sbeg = trkNmbs.cbegin();
      auto send = trkNmbs.cend();
      
      // the station and relative address
      int station, address;      
      // dt chamber identifiers
      int wheel_incr;
      int expectedwheel, dwheel, expectedsector, dsector;
/// JP : not sure if it's correct to match also the trkNmb (1)
//       unsigned expectedtrkNmb,dtrkNmb;

      for( ; mb != mbend; ++mb ) {                 
        
	for( ista = sbeg; ista != send; ++ista ) {
	  station = (ista - sbeg) + 1;
	  bool station_used = mode & ( 0x1 << (station-1) );
	  address = *ista;
          DTChamberId dtid = mb->first;
          if( !station_used || station != dtid.station() ) continue;
          wheel_incr = (isExtrapAcrossWheel(address,station) ? 1 : 0);
          expectedwheel = ( sp_wheel < 0 ? 
            wheel - wheel_incr :
            wheel + wheel_incr   );
	    
          dwheel = dtid.wheel();
          expectedsector = sector + relativeSector(address,station);
          expectedsector = ( expectedsector == 0 ? 12 : expectedsector);
          expectedsector = ( expectedsector == 13 ? 1 : expectedsector);
          dsector = dtid.sector();
          
/// JP : not sure if it's correct to match also the trkNmb (2)
//           expectedtrkNmb = address%2 + 1;
//           dtrkNmb = ;
          
          if( expectedsector == dsector &&
              expectedwheel  == dwheel            
/// JP : not sure if it's correct to match also the trkNmb (3)
//               && expectedtrkNmb == dtrkNmb
            ) {
               result[dtid] = L1ITMu::MBLTCollection(dtid);
          }
	}	
      }
      return result;
    }
  }
}
