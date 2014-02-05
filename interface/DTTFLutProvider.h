#ifndef __L1ITMU_DTTFLUTPROVIDER_H__
#define __L1ITMU_DTTFLUTPROVIDER_H__

// 
// Class: L1ITMBPtLutPlots
//
// Info: Service to provide DTTF LUTs for Barrel Upgrade Studies  
//
// Author: G. Codispoti, C. Battilana
//

#include <map>
#include <string>

class DTTFLutProvider {

 public:

  DTTFLutProvider() {};
  ~DTTFLutProvider() {};
 private:

  std::map< int, std::map< int, float > > phi_;
  std::map< int, std::map< int, float > > phib_;

  void parse( const std::string & inputdir );

  int parseLine( const std::string & line, int & inCh, int & outCh,
		 std::string & ref1, std::string & ref2,
		 float & pt, float & res, float & eff );

  void parse( const std::string & inputdir, const std::string & param,
	      std::map< int, std::map< int, float > > & effMap );

};

#endif

