#include "L1Trigger/L1IntegratedMuonTrigger/interface/DTTFId.h"

#include<sstream>

DTTFId::DTTFId(const DTTFId& id) 
{

  _wh = id._wh;
  _sec = id._sec;
  
}

std::string DTTFId::name() const
{ 
  
  std::stringstream name;
  
  name << "Wh" << _wh << "Sc" << _sec;
  
  return name.str();
  
}


