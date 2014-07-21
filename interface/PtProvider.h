
#ifndef __L1ITMU_PTPROVIDER_H__
#define __L1ITMU_PTPROVIDER_H__

// 
// Class: PtProvider
//
// Info: Service to provide Pt from DTTF LUTs and weights for Barrel Upgrade Studies  
//
// Author: G. Pinna 

#include <map>
#include <string>

class PtProvider {

 public:

  PtProvider() {};
  ~PtProvider() {};

  static void SetWeight( float & whBIn, float &whBOut, float & whDPhi, int mb1Obj, int mb2Obj, std::string PlaceName,std::string WeightType, std::string MeanType);
  static float ptLut(std::string LutName,float PhiValue, int  QualIn, int QualOut);

  //private:


};

#endif
