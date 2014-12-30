// 
// Class: PtProvider
//
// Info: Service to provide Pt from DTTF LUTs and weights for Barrel Upgrade Studies  
//
// Author: G. Pinna 

#include "L1Trigger/L1IntegratedMuonTrigger/interface/PtProvider.h"
#include "L1Trigger/L1IntegratedMuonTrigger/interface/ChambPairId.h"
#include <iostream>
#include <fstream>
#include <sstream>


 float PtProvider::ptLut(std::string LutName,float PhiValue, int QualIn, int QualOut)
{

  bool PhiFound = 0;
  LutName="../Luts/lut_data/"+LutName;
  
  std::ifstream datafileIn(LutName.c_str(),std::ifstream::in);
  
  if (datafileIn.fail())
    {
      std::cout<<("Error opening "+LutName).c_str()<<std::endl;
      abort();
    }
  
  if (datafileIn.peek() == std::ifstream::traits_type::eof())
    {
      std::cout<<("File "+LutName+" is empty").c_str()<<std::endl;
      abort();
    }
  
  float PtLut = -1;
  float Phi1 =1000;  //Just a value higher then any Lut value.
  float Phi2 =0;
  std::string line;

  std::string QualInString, QualOutString;

  switch(QualIn){

  case 0: QualInString = "DTIN";
    break;
  case 1: QualInString = "DTCORR";
    break;
  case 2: QualInString = "DTDIRR";
    break;
  case 3: QualInString = "DTOUT";
    break;
  case 4: QualInString = "NONE";
    break;
  }

   switch(QualOut){

  case 0: QualOutString = "DTIN";
    break;
  case 1: QualOutString = "DTCORR";
    break;
  case 2: QualOutString = "DTDIRR";
    break;
  case 3: QualOutString = "DTOUT";
    break;
  case 4: QualOutString = "NONE";
    break;
  }

 float  PtLutFirst=8;

while (getline(datafileIn,line)){
		  
		
		  std::stringstream readline; 
		  std::string A, D1, D2;
		  readline<<line;
		  readline>>A>>A>>D1>>D2>>PtLut>>Phi2>>A;
    
		  if(D1 != QualInString && D2 != QualOutString ) continue;
		  if(PtLut<8) continue;    //GP skip first LUT values
    
		  if(PhiValue>=Phi2&&PhiValue<Phi1){
		    PhiFound=1;
		    break;
		  }
		  PtLutFirst = PtLut;
		  Phi1=Phi2; 
 }

 if(PhiFound)return PtLutFirst; 
 else return 120.; 
}


void PtProvider::SetWeight(float & whBIn, float &whBOut, float & whDPhi, int mb1Obj, int mb2Obj, std::string PlaceName, std::string WeightType, std::string MeanType){
  std::string  FileName;
  
  //std::cout<<" "<<MeanType<<std::endl;
  
  FileName = "../test/Weights/"+MeanType+PlaceName+".txt"; 
  
  std::ifstream datafileIn(FileName.c_str(),std::ifstream::in);
  
  if (datafileIn.fail())
    {
      std::cout<<("Error opening "+FileName).c_str()<<std::endl;
      abort();
    }
  
  std::string line;
  
  std::string QualInString, QualOutString;
  
  switch(mb1Obj){
    
  case 0: QualInString = "DTIN";
    break;
  case 1: QualInString = "DTCORR";
    break;
  case 2: QualInString = "DTDIRR";
    break;
  case 3: QualInString = "DTOUT";
    break;
  case 4: QualInString = "NONE";
    break;
  default: std::cout<<"[L1ITMBPtEff]::SetWeight: WRONG value in weight file.\n";
    break;
}
  switch(mb2Obj){  
  case 0: QualOutString = "DTIN";
    break;
  case 1: QualOutString = "DTCORR";
    break;
  case 2: QualOutString = "DTDIRR";
    break;
  case 3: QualOutString = "DTOUT";
    break;
  case 4: QualOutString = "NONE";
    break;
  default: std::cout<<"[L1ITMBPtEff]::SetWeight: WRONG value in weight file.\n";
    break;
  }
   
  while (getline(datafileIn,line)){
    std::stringstream readline; 
    std::string D1, D2, Type;
    readline<<line;
    
    readline>>Type>>D1>>D2>>whDPhi>>whBIn>>whBOut;
   
    if(Type=="#") continue;          
    if(D2 != QualOutString) continue;
    if(D1 != QualInString ) continue;   
    if(Type != WeightType) continue; 
    
    break; 
  }
 
}
