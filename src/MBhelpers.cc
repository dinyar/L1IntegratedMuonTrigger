

// #include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

// #include "DataFormats/MuonDetId/interface/CSCDetId.h"
// #include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBhelpers.h"

namespace {
  // DTTF relative segment address explanation
  //
  // a code of '15' (stations 2,3,4) or '3' (station 1 only)
  // means there was no valid track extrapolation in this DT station
  //
  // schematic diagram of DT codes:
  //
  // STATION MB 1 
  // 
  // 3 => no primitive in this station
  //         -----------------------------
  //         |             |             |   ( next sector )
  //      P  -----------------------------
  //      H  |/// 2 (1) ///|             |   ( this sector )
  //      I  -----------------------------
  //         |             |             |   ( prev sector )
  //         -----------------------------
  //           this Wheel     next Wheel
  //           
  // STATION MB 2-3-4 
  // 
  // 15 => no primitive in this station
  //         -----------------------------
  //         |    12 (13)  |    4 (5)    |   ( next sector )
  //      P  -----------------------------
  //      H  |/// 8 (9) ///|    0 (1)    |   ( this sector )
  //      I  -----------------------------
  //         |    10 (11)  |    2 (3)    |   ( prev sector )
  //         -----------------------------
  //           this Wheel     next Wheel
  //  
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
    
    MBLTVectorRef
    getPrimitivesByMBTriggerInfo(const int wheel,
				  const int sp_wheel,
                                 const int sector,
                                 const int bx,
				  const edm::Handle<MBLTContainer>& mbs,
				  const unsigned mode,
				  const std::vector<unsigned>& trkNmbs) {
     MBLTVectorRef result;

      auto mb = mbs->cbegin();
      auto mbend = mbs->cend();
      
      std::vector<unsigned>::const_iterator ista;
      auto sbeg = trkNmbs.cbegin();
      auto send = trkNmbs.cend();
      
      // the station and relative address
      int station, address;      
      // dt chamber identifiers
      int wheel_incr;
      int expectedwheel, dwheel, expectedsector, dsector, expectedtr2tag, expectedbx;

      for ( unsigned int idx = 0; mb != mbend; ++mb, ++idx ) {
        
	for ( ista = sbeg; ista != send; ++ista ) {
	  station = (ista - sbeg) + 1;
	  bool station_used = mode & ( 0x1 << (station-1) );
	  address = *ista;
          DTChamberId dtid      = mb->first;
          MBLTCollection mbcoll = mb->second;
          if ( !station_used || station != dtid.station() ) continue;
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
          expectedtr2tag = address%2;
          expectedbx     = expectedtr2tag + bx; // If it is a 2nd track the trigger primitive bx is shifted (+ 1).
          
          if ( expectedsector == dsector &&
	       expectedwheel  == dwheel  ) {
            
            unsigned int tpmatched = 0;

//             result.push_back( MBLTContainerRef( mbs, idx ) );
//                 
            edm::LogWarning( "\tDTTF-MBLTCollection Matching:" ) 
                  << "\n\tMBLTCollection for ( Wheel / Station / Sector ) = ( " << dwheel << " / " << station << " / " << dsector << " )"
                  << "\n\t\tDTTriggerPrimitives matched : ";

            size_t dtListSize = mbcoll.getDtSegments().size();
            for ( size_t iDt = 0; iDt < dtListSize; ++iDt ) {
              
              const TriggerPrimitiveRef & dt = mbcoll.getDtSegments().at(iDt);
                            
              if ( dt->getDTData().qualityCode < 0) continue; // skipping theta chambers
              if ( expectedbx   != dt->getDTData().bx ) continue;  // match also in bx (correctly shifted according to the Tr2Tag)
              if ( expectedtr2tag != dt->getDTData().Ts2TagCode ) continue;  // match 1st / 2nd sub per Bx (if we have a 1st track at BX=0  and a 2nd track at BX=0, it means that the second track actually is at BX=-1)
              edm::LogWarning( "\t\tDTTriggerPrimitive" ) 
                  << "\n\t\tDTstub["<< iDt <<"]->wheel          = " << dt->getDTData().wheel
                  << "\n\t\tDTstub["<< iDt <<"]->station        = " << dt->getDTData().station
                  << "\n\t\tDTstub["<< iDt <<"]->sector         = " << dt->getDTData().sector
                  << "\n\t\tDTstub["<< iDt <<"]->qualityCode    = " << dt->getDTData().qualityCode
                  << "\n\t\tDTstub["<< iDt <<"]->bx             = " << dt->getDTData().bx
                  << "\n\t\tDTstub["<< iDt <<"]->Ts2TagCode     = " << dt->getDTData().Ts2TagCode;     
                  
              MBLTContainerRefPair mbdtmatch = std::make_pair(MBLTContainerRef( mbs, idx ), dt);              
                            
              result.push_back( mbdtmatch );
                  
              tpmatched++;
            }
            
            if (tpmatched != 1) {
              edm::LogError( "E R R O R ! ! !" )
                  << "\nFOUND " << tpmatched << " DT PRIMITIVES ASSOCIATED TO A SINGLE CHAMBER (via DTTF code)!\n";                 
            } 
             
          }
	}	
      }
      return result;
    }
  }
}
