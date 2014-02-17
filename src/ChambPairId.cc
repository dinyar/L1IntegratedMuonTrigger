#include "L1Trigger/L1IntegratedMuonTrigger/interface/ChambPairId.h"

#include <sstream>


ChambPairId::ChambPairId(DTTFId dttf, int inCh, int outCh, int inChObj, int outChObj) :
  _dttfId(dttf), _inCh(inCh), _outCh(outCh), _inChObj(inChObj), _outChObj(outChObj) 
{

  MBPtChambObjectName[DTIN]   = "DTIN";
  MBPtChambObjectName[DTCORR] = "DTCORR";
  MBPtChambObjectName[DTDIR]  = "DTDIR";
  MBPtChambObjectName[DTOUT]  = "DTOUT";
  MBPtChambObjectName[NONE]   = "NONE";
    
}


ChambPairId::ChambPairId(int wh, int sec, int inCh, int outCh, int inChObj, int outChObj) :
  _dttfId(wh,sec), _inCh(inCh), _outCh(outCh), _inChObj(inChObj), _outChObj(outChObj) 
{

  MBPtChambObjectName[DTIN]   = "DTIN";
  MBPtChambObjectName[DTCORR] = "DTCORR";
  MBPtChambObjectName[DTDIR]  = "DTDIR";
  MBPtChambObjectName[DTOUT]  = "DTOUT";
  MBPtChambObjectName[NONE]   = "NONE";
    
}


ChambPairId::ChambPairId(const ChambPairId & id) :
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


int ChambPairId::rawId() const
{ 

  int id = _dttfId.rawId() + 1000*_inCh + 10000*_outCh 
           + 100000*_inChObj + 10000000*_outChObj;
  
  return id;
  
}


std::string ChambPairId::name() const
{ 
  
  std::stringstream name;
  
  name << _dttfId.name() 
       << "inCh" << _inCh << "outCh" <<_outCh 
       << inObjName() << outObjName();
  
  return name.str();
  
}

ChambPairId::chamb_objects ChambPairId::chambFromString( const std::string & str )
{
  if ( str == "DTCORR" ) return  DTCORR;
  if ( str == "DTDIR" )  return  DTDIR;
  if ( str == "DTIN" )   return  DTIN;
  if ( str == "DTOUT" )  return  DTOUT;
  if ( str == "NONE" )   return  NONE;
  return NONE;
}




