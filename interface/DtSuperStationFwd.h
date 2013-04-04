#ifndef __L1ITMU_DTSUPERSTATION_H__
#define __L1ITMU_DTSUPERSTATION_H__

#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"

class DTChamberId;

namespace L1ITMu {
  class DtSuperStation;

  typedef std::pair<DTChamberId, DtSuperStation> DtSuperStationPair;
  typedef std::map<DTChamberId, DtSuperStation> DtSuperStationMap;
}

#endif
