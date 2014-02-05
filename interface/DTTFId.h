#ifndef __L1ITMU_DTTFID_H__
#define __L1ITMU_DTTFID_H__

#include<string>

// --------------------------------------------------
// Id class to handle DTTF identification
// --------------------------------------------------

class DTTFId {

public:

  DTTFId(int wh, int sec) : _wh(wh), _sec(sec) { };
  DTTFId(const DTTFId & id);

  ~DTTFId() { };

  int rawId() const { return (_wh+3) + 10*_sec; } ;
  std::string name() const;

private:
  
  int _wh, _sec;

};

#endif
