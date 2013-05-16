#ifndef __L1ITMU_MBTrack_H__
#define __L1ITMU_MBTrack_H__

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

namespace L1ITMu {
  class MBTrack;
  typedef std::vector<MBTrack> MBTrackCollection;
  typedef edm::Ref<MBTrackCollection> MBTrackCollectionRef;

}

#endif
