#include <stdlib.h>
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TLine.h"
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

//---------------------------------------------------------------//
//Plot rate vs pt cut for different pt assignment algorithms
//Author: Gian Luca Pinna
//---------------------------------------------------------------//


std::string WeightExt = "Weights"; 
 
void SetWeight(float&, float&, float&, int, int, std::string, std::string, std::string);

typedef struct {Float_t PtBIn,PtBOut,PtDPhi,GMTPt,Wheel_,Qual1_,Qual2_;} PTs_;

void RateStudies(string FileIn,string WeightFolder="Weights",Int_t Wheel = 0){

  WeightExt=WeightFolder;

  TFile *f = new TFile(FileIn.c_str());
  TTree *t1 = (TTree*)f->Get("L1ITMuPtPlotter/ElderTree");

  PTs_ pts;
  t1->SetBranchAddress("Pts",&pts.PtBIn);

  //TH1F *hRateBIn= new TH1F("hRateBIn"," ", 60,-1.,119.);
  //TH1F *hRateBOut= new TH1F("hRateBOut"," ", 60,-1.,119.);
  TH1F *hRateDTTF= new TH1F("hRateDTTF"," ", 60,-1.,119.);
  TH1F *hRateDPhi= new TH1F("hRateBDPhi"," ", 60,-1.,119.);
  TH1F *hRatePtMin= new TH1F("hRatePtMin"," ", 60,-1.,119.);
  TH1F *hRatePtMin2= new TH1F("hRatePtMin2"," ", 60,-1.,119.);
  TH1F *hRatePtMin3= new TH1F("hRatePtMin3"," ", 60,-1.,119.);
  TH1F *hRatePtWeighted= new TH1F("hRatePtWeighted"," ", 60,-1.,119.);
  TH1F *hRateNewPt= new TH1F("hRateNewPt"," ", 60,-1.,119.);


  Int_t nentries = (Int_t)t1->GetEntries();
  Int_t nevents = 0;
  for (Int_t i=0; i<nentries; i++) {
    
    t1->GetEntry(i);

    if(Wheel==0){
      if(abs(pts.Wheel_)==2) continue; 
    }
    else if(pts.Wheel_!=Wheel) continue;  

    nevents++; 

    Float_t whWMTailBIn = 1;
    Float_t whWMTailBOut = 1;
    Float_t whWMTailDPhi = 1;

    Float_t whWMInTailBIn = 1;
    Float_t whWMInTailBOut = 1;
    Float_t whWMInTailDPhi = 1;

    Float_t whWMOutTailBIn = 1;
    Float_t whWMOutTailBOut = 1;
    Float_t whWMOutTailDPhi = 1;

    
    //example for RMS Weight
    // Float_t whWMRMSBIn = 1;
    // Float_t whWMRMSBOut = 1;
    // Float_t whWMRMSDPhi = 1;

    int sector = 0;
    stringstream PlaceName;
    PlaceName << "Wh" << pts.Wheel_ << "Sc" << sector;

    //example for RMS Weight
    //SetWeight( whWMTailBIn, whWMTailBOut, whWMTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMean","RMS");
    //SetWeight( whWMInTailBIn, whWMInTailBOut, whWMInTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanIn","RMS");
    //SetWeight( whWMOutTailBIn, whWMOutTailBOut, whWMOutTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanOut","RMS");
      

    //Weight for "tail" method
    SetWeight( whWMTailBIn, whWMTailBOut, whWMTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMean","Tail");
    SetWeight( whWMInTailBIn, whWMInTailBOut, whWMInTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanIn","Tail");
    SetWeight( whWMOutTailBIn, whWMOutTailBOut, whWMOutTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanOut","Tail");
    
    Float_t  WeightedPt = pts.PtDPhi;

    //Choice 1
    //if((pts.Qual1_==1) & (pts.Qual2_==1)) WeightedPt = whWMTailBIn*pts.PtBIn+whWMTailBOut*pts.PtBOut+whWMTailDPhi*pts.PtDPhi;
    //else if((pts.Qual1_==1) & (pts.Qual2_!=1)) WeightedPt = whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;
    //else if((pts.Qual1_!=1) & (pts.Qual2_==1)) WeightedPt =  whWMOutTailBOut*pts.PtBOut+whWMOutTailDPhi*pts.PtDPhi; 
    
    //Choice 2 
    if((pts.Qual1_==1)) WeightedPt =  whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;

    //Choice 3
    //if((pts.Qual1_==1)) WeightedPt =  whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;
    //else if((pts.Qual1_!=1) & (pts.Qual2_==1)) WeightedPt =  whWMOutTailBOut*pts.PtBOut+whWMOutTailDPhi*pts.PtDPhi; 
    
    
    
    Float_t NewPt = pts.PtDPhi;
    if(pts.Qual1_==1) NewPt =  pts.PtBIn;
    //else if((pts.Qual1_!=1) & (pts.Qual2_==1)) WeightedPt = pts.PtBOut; // Get Worse       
    
    Float_t MinPt = pts.PtDPhi;
    
    if((pts.Qual1_==1) & (pts.Qual2_==1)){  
      if(pts.PtBIn<MinPt) MinPt = pts.PtBIn;
      if(pts.PtBOut<MinPt) MinPt=pts.PtBOut;
    }
    else if((pts.Qual1_==1) & (pts.Qual2_!=1)){
      if(pts.PtBIn<MinPt) MinPt = pts.PtBIn;
    }
    else if((pts.Qual1_==1) & (pts.Qual2_!=1)){
      
      if(pts.PtBOut<MinPt) MinPt=pts.PtBOut;
    }
    
    
    Float_t MinPt2=pts.PtDPhi;;
    
    if(pts.Qual1_==1 &&  pts.Qual2_==1){
      if(fabs(pts.PtDPhi-pts.PtBIn)<120 && (pts.PtDPhi-pts.PtBOut)>15){
	if(pts.PtBIn<MinPt2) MinPt2=pts.PtBIn;
      }
      
      else if(abs(pts.PtDPhi-pts.PtBOut)<120 &&(pts.PtDPhi-pts.PtBIn)>15 ){
	if(pts.PtBOut<MinPt2) MinPt2=pts.PtBOut;
      }
    }
    
    else if(pts.Qual1_==1&&pts.Qual2_!=1){
      if(pts.PtBIn<MinPt2) MinPt2=pts.PtBIn;
    }
    
    else if(pts.Qual1_!=1&&pts.Qual2_==1){
      if(pts.PtBOut<MinPt2) MinPt2=pts.PtBOut;
    }
    
    
  Float_t MinPt3=pts.PtDPhi;
  
  if(pts.Qual1_==1 &&  pts.Qual2_==1){
    if(fabs(pts.PtDPhi-pts.PtBIn)<60 && (pts.PtDPhi-pts.PtBOut)>30){
      if(pts.PtBIn<MinPt3) MinPt3=pts.PtBIn;
    } 
    else if(abs(pts.PtDPhi-pts.PtBOut)<60 &&(pts.PtDPhi-pts.PtBIn)>30 ){
      if(pts.PtBOut<MinPt3) MinPt3=pts.PtBOut;
    }
  }
  else if(pts.Qual1_==1&&pts.Qual2_!=1){
    if(pts.PtBIn<MinPt3) MinPt3=pts.PtBIn;
  }
  
  else if(pts.Qual1_!=1&&pts.Qual2_==1){
    if(pts.PtBOut<MinPt3) MinPt3=pts.PtBOut;
  }

  for(int j = 0.; j<=120; j+=2){
    
    if(pts.GMTPt>=j)   hRateDTTF->Fill(j); 
    if(pts.PtDPhi>=j)  hRateDPhi->Fill(j);
    if(WeightedPt>=j)  hRatePtWeighted->Fill(j);
    if(NewPt>=j)   hRateNewPt->Fill(j);
    if(MinPt>=j)   hRatePtMin->Fill(j);
    if(MinPt2>=j)  hRatePtMin2->Fill(j);
    if(MinPt3>=j)  hRatePtMin3->Fill(j);  
 }
  }
  
  
  hRateDPhi->SetLineColor(2);
  hRateDPhi->SetTitle("Rate vs p_{T}; pT_{Cut} ; Rate"); 
  hRatePtWeighted->SetLineColor(93);
  hRateDTTF->SetLineColor(51);
  // hRateBIn->SetLineColor(2);
  // hRateBOut->SetLineColor(4);
  hRateNewPt->SetLineColor(64);
  hRatePtMin->SetLineColor(4);
  hRatePtMin2->SetLineColor(6);
  hRatePtMin3->SetLineColor(8);
  
  TLegend *label = new TLegend(.65,.12,.85,.26);
  
  label->AddEntry(hRateDTTF," DTTF ","l");
  label->AddEntry(hRateDPhi," #Delta #phi ","l");
  label->AddEntry(hRatePtMin," min pt","l");
  label->AddEntry(hRatePtMin2," min pt 2 ","l");
  label->AddEntry(hRatePtMin3," min pt 3 ","l");
  label->AddEntry(hRatePtWeighted," Tail method ","l");
  label->AddEntry(hRateNewPt," New method ","l");
  // label->AddEntry(hRateBIn," BIn ","p");
  // label->AddEntry(hRateBOut," BOut pt 2 ","p");
  
  TCanvas *c1 = new TCanvas("C1","Rate plots",500,500);
  c1->SetLogy();
  
  hRateDPhi->Draw("same");
  hRateDTTF->Draw("same");
  hRatePtMin->Draw("same");
  hRatePtWeighted->Draw("same");
  hRatePtMin2->Draw("same");
  hRatePtMin3->Draw("same");
  hRateNewPt->Draw("same");
  label->Draw("same");
  
}

void SetWeight(float & whBIn, float &whBOut, float & whDPhi, int mb1Obj, int mb2Obj, std::string PlaceName, std::string WeightType, std::string MeanType){
 std::string FileName;
  FileName = WeightExt+"/"+MeanType+PlaceName+".txt";
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
