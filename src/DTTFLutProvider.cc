// 
// Class: L1ITMBPtLutPlots
//
// Info: Service to provide DTTF LUTs for Barrel Upgrade Studies  
//
// Author: G. Codispoti, C. Battilana
//

#include "L1Trigger/L1IntegratedMuonTrigger/interface/DTTFLutProvider.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/ChambPairId.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/lexical_cast.hpp>

int DTTFLutProvider::parseLine( const std::string & line, int & inCh, int & outCh,
			      std::string & ref1, std::string & ref2,
			      float & pt, float & thr, float & eff )
{

  std::vector<std::string> elements;
  boost::algorithm::split(elements,line,boost::algorithm::is_any_of(" \t\n"));
  if (elements.size() != 7) {
    std::cerr << "[DTTPGParamsWriter] wrong number of entries in line : " << line << " pleas check your input file syntax!" << std::endl;
    return 0;
  } else {
    ref1 = elements[2];
    ref2 = elements[3];





    try {
      inCh = boost::lexical_cast<int>( elements[0] );
    } catch( boost::bad_lexical_cast const& ) {
      std::cerr << "Error: inCh " << elements[0] << " is not valid";
      return 0;
    }
    try {
      outCh = boost::lexical_cast<int>( elements[1] );
    } catch( boost::bad_lexical_cast const& ) {
      std::cerr << "Error: outCh " << elements[1] << " is not valid";
      return 0;
    }

    try {
      pt = boost::lexical_cast<float>( elements[4] );
    } catch( boost::bad_lexical_cast const& ) {
      std::cerr << "Error: pt " << elements[4] << " is not valid";
      return 0;
    }
    try {
      thr = boost::lexical_cast<float>( elements[5] );
    } catch( boost::bad_lexical_cast const& ) {
      std::cerr << "Error: phi/phib value " << elements[5] << " is not valid";
      return 0;
    }
    try {
      eff = boost::lexical_cast<float>( elements[6] );
    } catch( boost::bad_lexical_cast const& ) {
      std::cerr << "Error: efficiency " << elements[6] << " is not valid";
      return 0;
    }
  }
  return 1;
}




void DTTFLutProvider::parse( const std::string & inputdir, const std::string & param,
			     std::map< int, std::map< int, float > > & effMap )
{


  int wheels[6] = { -3, -2, -1, 1, 2, 3 };

  for ( int sector = 0; sector < 12; ++sector ) {
    for ( size_t w = 0; w < 6; ++w ) {

      DTTFId dttfId( wheels[w],  sector);

      std::ostringstream inputFileName ;
      inputFileName << inputdir << '/' << param << "Wh" << wheels[w] << "Sc" <<  sector;
      std::ifstream inputFile(inputFileName.str());
      int nLines=0;
      std::string line;

      while ( std::getline(inputFile, line) ) {

	int inCh;
	int outCh;
	std::string ref1;
	std::string ref2;
	float pt;
	float thr;
	float eff;

	if ( parseLine(line, inCh, outCh, ref1, ref2, pt, thr, eff ) ) {
	  int inChObj  = ChambPairId::chambFromString( ref1 );
	  int outChObj = ChambPairId::chambFromString( ref2 );
	  ChambPairId chId( wheels[w], sector, inCh, outCh, inChObj, outChObj );
	  effMap[dttfId.rawId()][chId.rawId()] = thr;
	  std::cout << wheels[w] << '\t' << sector << '\t' << inCh << '\t' << outCh
		    << '\t' << ref1 << " (" << inChObj << ")\t"
		    << '\t' << ref2 << " (" << outChObj << ")\t"
		    << pt << '\t' << thr << '\t' << eff << std::endl;
	}
	nLines++;
      }
      inputFile.close();
    }
  }
}


void DTTFLutProvider::parse( const std::string & inputdir )
{
  parse( inputdir, "Phi", phi_ );
  parse( inputdir, "PhiBend", phib_ );
}










