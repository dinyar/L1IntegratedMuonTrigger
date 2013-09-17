// 
// Class: L1ITMBPtLutPlots
//
// Info: Performs GEN-DTTF matching and computes LUT plots  
//
// Author: 
//


/* C++ Headers */
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/lexical_cast.hpp>



// --------------------------------------------------
// Id class to handle DTTF identification
// --------------------------------------------------

namespace sticazzi {

  class DTTFId {

  public:

    DTTFId(int wh, int sec);
    DTTFId(const DTTFId & id);

    ~DTTFId() { };

    int rawId() const { return (_wh+3) + 10*_sec; } ;
    std::string name() const;

  private:
  
    int _wh, _sec;

  };
}


sticazzi::DTTFId::DTTFId(int wh, int sec) : _wh(wh), _sec(sec) 
{
 
}


sticazzi::DTTFId::DTTFId(const DTTFId& id) 
{

  _wh = id._wh;
  _sec = id._sec;
  
}

std::string sticazzi::DTTFId::name() const
{ 
  
  std::stringstream name;
  
  name << "Wh" << _wh << "Sc" << _sec;
  
  return name.str();
  
}


// --------------------------------------------------
// Id class to handle chamber pair and the origin
// objects used to build them (e.g. inner RPC layer)
// --------------------------------------------------

namespace sticazzi {

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

}



sticazzi::ChambPairId::chamb_objects sticazzi::ChambPairId::chambFromString( const std::string & str )
{
  if ( str == "DTCORR" ) return  DTCORR;
  if ( str == "DTDIR" )  return  DTDIR;
  if ( str == "DTIN" )   return  DTIN;
  if ( str == "DTOUT" )  return  DTOUT;
  if ( str == "NONE" )   return  NONE;
  return NONE;
}




sticazzi::ChambPairId::ChambPairId(DTTFId dttf, int inCh, int outCh, int inChObj, int outChObj) :
  _dttfId(dttf), _inCh(inCh), _outCh(outCh), _inChObj(inChObj), _outChObj(outChObj) 
{

  MBPtChambObjectName[DTIN]   = "DTIN";
  MBPtChambObjectName[DTCORR] = "DTCORR";
  MBPtChambObjectName[DTDIR]  = "DTDIR";
  MBPtChambObjectName[DTOUT]  = "DTOUT";
  MBPtChambObjectName[NONE]   = "NONE";
    
}


sticazzi::ChambPairId::ChambPairId(int wh, int sec, int inCh, int outCh, int inChObj, int outChObj) :
  _dttfId(wh,sec), _inCh(inCh), _outCh(outCh), _inChObj(inChObj), _outChObj(outChObj) 
{

  MBPtChambObjectName[DTIN]   = "DTIN";
  MBPtChambObjectName[DTCORR] = "DTCORR";
  MBPtChambObjectName[DTDIR]  = "DTDIR";
  MBPtChambObjectName[DTOUT]  = "DTOUT";
  MBPtChambObjectName[NONE]   = "NONE";
    
}


sticazzi::ChambPairId::ChambPairId(const ChambPairId & id) :
  _dttfId(id._dttfId)
{
  
  _inCh = id._inCh;
  _outCh = id._outCh;
  _inChObj = id._inChObj;
  _outChObj = id._outChObj;

  MBPtChambObjectName[DTIN]   =   id.MBPtChambObjectName[DTIN];
  MBPtChambObjectName[DTCORR] =   id.MBPtChambObjectName[DTCORR];
  MBPtChambObjectName[DTDIR]  =   id.MBPtChambObjectName[DTDIR];
  MBPtChambObjectName[DTOUT]  =   id.MBPtChambObjectName[DTOUT];
  MBPtChambObjectName[NONE]   =   id.MBPtChambObjectName[NONE];
  
}


int sticazzi::ChambPairId::rawId() const
{ 

  int id = _dttfId.rawId() + 1000*_inCh + 10000*_outCh 
           + 100000*_inChObj + 10000000*_outChObj;
  
  return id;
  
}


std::string sticazzi::ChambPairId::name() const
{ 
  
  std::stringstream name;
  
  name << _dttfId.name() 
       << "inCh" << _inCh << "outCh" <<_outCh 
       << inObjName() << outObjName();
  
  return name.str();
  
}




////////////////////////////////////////

namespace sticazzi {

  class DTTFLutReader {

  public:

    DTTFLutReader() {};
    ~DTTFLutReader() {};
    void parse( const std::string & inputdir );

  private:

    std::map< int, std::map< int, float > > phi_;
    std::map< int, std::map< int, float > > phib_;
    int parseLine( const std::string & line, int & inCh, int & outCh,
		   std::string & ref1, std::string & ref2,
		   float & pt, float & res, float & eff );

    void parse( const std::string & inputdir, const std::string & param,
		std::map< int, std::map< int, float > > & effMap );

  };

}






int sticazzi::DTTFLutReader::parseLine( const std::string & line, int & inCh, int & outCh,
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




void sticazzi::DTTFLutReader::parse( const std::string & inputdir, const std::string & param,
				     std::map< int, std::map< int, float > > & effMap )
{


  int wheels[6] = { -3, -2, -1, 1, 2, 3 };

  for ( int sector = 0; sector < 12; ++sector ) {
    for ( size_t w = 0; w < 6; ++w ) {

      sticazzi::DTTFId dttfId( wheels[w],  sector);

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
	  int inChObj = sticazzi::ChambPairId::chambFromString( ref1 );
	  int outChObj = sticazzi::ChambPairId::chambFromString( ref2 );
	  sticazzi::ChambPairId chId( wheels[w], sector, inCh, outCh, inChObj, outChObj );
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


void sticazzi::DTTFLutReader::parse( const std::string & inputdir )
{
  parse( inputdir, "Phi", phi_ );
  parse( inputdir, "PhiBend", phib_ );
}








int main(int argc, char** argv) {


  if ( argc < 2 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <input>"
              << std::endl << std::endl;
    std::cout << "\tes. : " << argv[0] << " efficiency\n" << std::endl;
    return -1;
  }

  std::string inputdir = argv[1];
  sticazzi::DTTFLutReader reader;
  reader.parse( inputdir );

}

