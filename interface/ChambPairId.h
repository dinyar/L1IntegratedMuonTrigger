#ifndef __L1ITMU_CHAMBPAIRID_H__
#define __L1ITMU_CHAMBPAIRID_H__

// --------------------------------------------------
// Id class to handle chamber pair and the origin
// objects used to build them (e.g. inner RPC layer)
// --------------------------------------------------

#include "L1Trigger/L1IntegratedMuonTrigger/interface/DTTFId.h"

#include <string>

class ChambPairId {

public:

  enum chamb_objects {DTIN=0, DTCORR, DTDIR, DTOUT, NONE};

  ChambPairId(int wh, int sec, int inCh, int outCh, int inChObj, int outChObj);
  ChambPairId(DTTFId dttf, int inCh, int outCh, int inChObj, int outChObj);
  ChambPairId(const ChambPairId & id);

  ~ChambPairId() { };
  
  int rawId() const;
  std::string name() const;

  const DTTFId & dttfId() const { return _dttfId; };

  std::string inObjName()  const { return MBPtChambObjectName[_inChObj]; }; 
  std::string outObjName() const { return MBPtChambObjectName[_outChObj]; };

  static chamb_objects chambFromString( const std::string & str );
  
private:
  
  DTTFId _dttfId;
  int _inCh, _outCh;
  int _inChObj, _outChObj;

  std::string MBPtChambObjectName[5];
  
};

#endif
