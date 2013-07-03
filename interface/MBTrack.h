#ifndef __L1ITMU_MBTRACK_H__
#define __L1ITMU_MBTRACK_H__
// 
// Class: L1ITMu::MBTrack
//
// Info: This class represents (one of the) internal tracks processed
//       by L1ITMu before sending off to the GT or whatever comes after.
//       As such, they are similar to the L1MuRegionalCands in terms of
//       their meaning.
//       To exploit that we also allow these tracks to be built out of
//       old style regional cands so we can perform studies using them
//       in the new framework. So that debugging is easier in the case
//       of the new track not matching it's parent, a reference to the
//       parent track is stored.
//
// Author: 
//

#include <iostream>

#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitiveFwd.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/TriggerPrimitive.h"

#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollection.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/MBLTCollectionFwd.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/RegionalTracksFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTTrackCand.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1Track.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"

// Muons & Tracks Data Formats
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// GMT Data Formats
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutRecord.h"


namespace L1ITMu{
  class MBTrack : public L1MuRegionalCand {   
  public:

    // default constructor
    MBTrack():_wheel(0),_sector(0),_type(5) {}
    ~MBTrack() { _gmtout.clear(); _gmtin.clear(); }

    ///  construct starting from a DTTF candidate
    MBTrack(const L1MuDTTrackCand&);

// // //     // return the persistent pointer to the parent of this internal track
// // //     // may be null if this has no parent
// // //     RegionalCandBaseRef parent() const { return _parent; }
// // // 
// // //     /// set set parent in case you are starting from primitives
// // //     void setParent(const RegionalCandBaseRef& parent)
// // //        { _parent = parent; }

    // return the persistent pointer to the parent of this internal track
    // may be null if this has no parent
    const L1MuDTTrackCand& parent() const { return _parent; }

    /// set set parent in case you are starting from primitives
    void setParent(const L1MuDTTrackCand& parent)
       { _parent = parent; }       
       
    /// add primitives in case you are starting from tracks
    void addStub(const MBLTContainerRefPair& stub);
    /*     void addStub(const std::pair<const DTChamberId, L1ITMu::MBLTCollection>& stub); */

    /// return list of associated MBLT primitives
    const MBLTVectorRef& getStubs() const 
      { return _associatedStubs; }

    /// associate a GMTout 
    void associateGMTout(const L1MuGMTExtendedCand& gmtout);
   
    /// associate a GMTin 
    void associateGMTin(const L1MuRegionalCand& gmtin);
      
    /// return associated GMTout, if any
    const std::vector<L1MuGMTExtendedCand> getAssociatedGMTout() const 
      { return _gmtout; }
    
    /// return associated GMTin, if any
    const std::vector<L1MuRegionalCand> getAssociatedGMTin() const 
      { return _gmtin; }

  private:
    int _wheel, _sector;
    unsigned _type;
    //pointer to parent, if this was created from a CSC/DT/RPC track
    MBLTVectorRef _associatedStubs;
    L1MuDTTrackCand _parent;
    std::vector<L1MuGMTExtendedCand> _gmtout;
    std::vector<L1MuRegionalCand> _gmtin;    
  };
}

#endif
