#include <stdlib.h>

#include <memory>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <ostream>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

void RMS(	   Float_t &RMSPhiBIn,
		   Float_t &RMSPhiBOut,
		   Float_t &RMSDeltaPhi,
		   Float_t &RMSGMT, string DT1 , string DT2, string wh, string sc, string cut);



void WeightCreator(string wh, string sc, string cut){ 

  Float_t RMSPhiBIn;
  Float_t RMSPhiBOut;
  Float_t RMSDeltaPhi;
  Float_t RMSGMT;
 
  string outputfile = "Weight/Wh" +wh+"Sc"+sc+".txt";
 
  ofstream outputFile(outputfile.c_str());
  outputFile<<"# DT1 DT2 RMSGMT RMSDeltaPhi RMSPhiBIn RMSPhiBOut "<<endl;

  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTCORR","DTCORR",wh, sc, cut);
  outputFile<<"DTCORR"<<" "<<"DTCORR"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
 
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTCORR","DTIN",wh, sc, cut);
  outputFile<<"DTCORR"<<" "<<"DTIN"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTCORR","DTOUT",wh, sc, cut);
  outputFile<<"DTCORR"<<" "<<"DTOUT"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTIN","DTCORR",wh, sc, cut);
  outputFile<<"DTIN"<<" "<<"DTCORR"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTIN","DTIN",wh, sc, cut);
  outputFile<<"DTIN"<<" "<<"DTIN"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTIN","DTOUT",wh, sc, cut);
  outputFile<<"DTIN"<<" "<<"DTOUT"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTOUT","DTCORR",wh, sc, cut);
  outputFile<<"DTOUT"<<" "<<"DTCORR"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTOUT","DTIN",wh, sc, cut);
  outputFile<<"DTOUT"<<" "<<"DTIN"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;
  
  RMS(RMSPhiBIn,RMSPhiBOut,RMSDeltaPhi,RMSGMT,"DTOUT","DTOUT",wh, sc, cut);
  outputFile<<"DTOUT"<<" "<<"DTOUT"<<" "<<RMSGMT<<" "<<" "<<RMSDeltaPhi<<" "<<RMSPhiBIn<<" "<<RMSPhiBOut<<endl;

  outputFile.close();
}

void RMS(	   Float_t &oRMSPhiBIn,
		   Float_t &oRMSPhiBOut,
		   Float_t &oRMSDeltaPhi,
		   Float_t &oRMSGMT, string DT1 , string DT2, string wh, string sc, string cut){
  
  string typeGMT  ="GMTPtResol";
  string typePhi  ="DeltaPhiPtResol";
  string typeBIn  ="PhiBInPtResol";
  string typeBOut ="PhiBOutPtResol";
  
  string namefileGMT = "h"+typeGMT+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut;
  string namepathGMT = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/"+namefileGMT;  

  string namefileBIn = "h"+typeBIn+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut;
  string namepathBIn = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/"+namefileBIn;

  string namefileBOut = "h"+typeBOut+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut;
  string namepathBOut = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/"+namefileBOut;

  string namefilePhi = "h"+typePhi+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut;
  string namepathPhi = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/"+namefilePhi;
  
  string namefile = "../L1ITMuonBarrelPtStudies.root";
                      
  TFile *file1 = new TFile(namefile.c_str());

  TH1F* HresGMT = (TH1F*)file1->Get(namepathGMT.c_str());

  TH1F* HresPhiBIn = (TH1F*)file1->Get(namepathBIn.c_str());

  TH1F* HresPhiBOut = (TH1F*)file1->Get(namepathBOut.c_str());

  TH1F* HresDeltaPhi = (TH1F*)file1->Get(namepathPhi.c_str());

  oRMSPhiBIn = HresPhiBIn->GetRMS();
  oRMSPhiBOut = HresPhiBOut->GetRMS();
  oRMSDeltaPhi = HresDeltaPhi->GetRMS();
  oRMSGMT = HresGMT->GetRMS();  

}


