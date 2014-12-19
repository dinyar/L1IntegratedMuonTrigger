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

std::string WeightExt = "Weights"; 

void SetWeight(float&, float&, float&, int, int, std::string, std::string, std::string);

typedef struct {Float_t PtTrue,PtBIn,PtBOut,PtDPhi,PtGMT,Wheel_,Qual1_,Qual2_;} PTs_;

void EfficiencyStudies(string FileIn, Float_t dttfPtCut =25.,string WeightFolder="Weights", Int_t Wheel = 0 ){

WeightExt=WeightFolder;

  TFile *f = new TFile(FileIn.c_str());
  TTree *t1 = (TTree*)f->Get("L1ITMuPtPlotter/ElderTree");

  PTs_ pts;
  t1->SetBranchAddress("Pts",&pts.PtTrue);

  TEfficiency *effBIn= new TEfficiency("effBIn"," ", 60,-1.,119.);
  TEfficiency *effBOut= new TEfficiency("effBOut"," ", 60,-1.,119.);
  TEfficiency *effDPhi= new TEfficiency("effBDPhi"," ", 60,-1.,119.);
  TEfficiency *effPtMin= new TEfficiency("effPtMin"," ", 60,-1.,119.);
  TEfficiency *effPtMin2= new TEfficiency("effPtMin2"," ", 60,-1.,119.);
  TEfficiency *effPtMin3= new TEfficiency("effPtMin3"," ", 60,-1.,119.);
  TEfficiency *effWeightPt= new TEfficiency("effWeightPt"," ", 60,-1.,119.);
  TEfficiency *effNewPt= new TEfficiency("effNewPt"," ", 60,-1.,119.);

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


    //Weight for "tail" method
     SetWeight( whWMTailBIn, whWMTailBOut, whWMTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMean","Tail");
     SetWeight( whWMInTailBIn, whWMInTailBOut, whWMInTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanIn","Tail");
     SetWeight( whWMOutTailBIn, whWMOutTailBOut, whWMOutTailDPhi, pts.Qual1_, pts.Qual2_, PlaceName.str(),"WMeanOut","Tail");

     Float_t  WeightedPt = pts.PtDPhi;

     //Choice 1
     if((pts.Qual1_==1) & (pts.Qual2_==1)) WeightedPt =  whWMTailBIn*pts.PtBIn+whWMTailBOut*pts.PtBOut+whWMTailDPhi*pts.PtDPhi;
     else if((pts.Qual1_==1) & (pts.Qual2_!=1)) WeightedPt =  whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;
     else if((pts.Qual1_==1) & (pts.Qual2_!=1)) WeightedPt =  whWMOutTailBOut*pts.PtBOut+whWMOutTailDPhi*pts.PtDPhi;
     
     //Choice 2 //Best
     //if((pts.Qual1_==1)) WeightedPt =  whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;
     
     //Choice 3
     //if((pts.Qual1_==1)) WeightedPt =  whWMInTailBIn*pts.PtBIn+whWMInTailDPhi*pts.PtDPhi;
     //else if((pts.Qual1_!=1) & (pts.Qual2_==1)) WeightedPt =  whWMOutTailBOut*pts.PtBOut+whWMOutTailDPhi*pts.PtDPhi; 
    

     
     Float_t NewPt = pts.PtDPhi;
     if(pts.Qual1_==1) NewPt =  pts.PtBIn;
     //else if((pts.Qual1_!=1) & (pts.Qual2_==1)) WeightedPt = pts.PtBOut; // Rate get Worse       
     
     
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
    
    Float_t MinPt2=pts.PtDPhi;
    
    if(pts.Qual1_==1 &&  pts.Qual2_==1){
      if(fabs(pts.PtDPhi-pts.PtBIn)<80 && (pts.PtDPhi-pts.PtBOut)>20){
	if(pts.PtBIn<MinPt2) MinPt2=pts.PtBIn;
      }
      
      else if(abs(pts.PtDPhi-pts.PtBOut)<80 &&(pts.PtDPhi-pts.PtBIn)>20 ){
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

  //Exemple for RMS weight method 
  //  Float_t WeightedPt =  whWMRMSBIn*pts.PtBIn+whWMRMSBOut*pts.PtBOut+whWMRMSDPhi*pts.PtDPhi;


  effWeightPt->Fill(WeightedPt>=dttfPtCut,pts.PtTrue);
  effBIn->Fill(pts.PtBIn>=dttfPtCut,pts.PtTrue);
  effBOut->Fill(pts.PtBOut>=dttfPtCut,pts.PtTrue);
  effDPhi->Fill(pts.PtDPhi>=dttfPtCut,pts.PtTrue);
  effPtMin->Fill(MinPt>=dttfPtCut,pts.PtTrue);
  effPtMin2->Fill(MinPt2>=dttfPtCut,pts.PtTrue);
  effPtMin3->Fill(MinPt3>=dttfPtCut,pts.PtTrue);
  effNewPt->Fill(NewPt>=dttfPtCut,pts.PtTrue);
  }

  cout<<"entries "<<nentries<<" central events "<<nevents<<endl;
  
  TLine *lineCut = new TLine(dttfPtCut,0,dttfPtCut,1.05);
  lineCut->SetLineColor(4);
  lineCut->SetLineStyle(4);
  lineCut->SetLineWidth(4);  

  effDPhi->SetMarkerStyle(21);
  effDPhi->SetMarkerColor(2);
  effDPhi->SetMarkerSize(.8);
 
  effBIn->SetMarkerStyle(20);
  effBIn->SetMarkerColor(91);
  effBIn->SetMarkerSize(.8);
 
  effBOut->SetMarkerStyle(20);
  effBOut->SetMarkerColor(95);
  effBOut->SetMarkerSize(.8); 

  effPtMin->SetMarkerStyle(21);
  effPtMin->SetMarkerColor(4);
  effPtMin->SetMarkerSize(.8);

  effPtMin2->SetMarkerStyle(21);
  effPtMin2->SetMarkerColor(6);
  effPtMin2->SetMarkerSize(.8); 

  effPtMin3->SetMarkerStyle(21);
  effPtMin3->SetMarkerColor(8);
  effPtMin3->SetMarkerSize(.8); 
  
  effWeightPt->SetMarkerStyle(22);
  effWeightPt->SetMarkerColor(93);
  effWeightPt->SetMarkerSize(.8); 

  effNewPt->SetMarkerStyle(22);
  effNewPt->SetMarkerColor(64);
  effNewPt->SetMarkerSize(.8); 


  TLegend *label = new TLegend(.65,.12,.85,.26);
  label->AddEntry(lineCut,"p_{T} threshold","l"); 
  label->AddEntry(effDPhi," #Delta #phi ","p");
  label->AddEntry(effPtMin," min pt","p");
  label->AddEntry(effWeightPt," Tail method","p");
  label->AddEntry(effPtMin2," min pt 2 ","p");
  label->AddEntry(effPtMin3," min pt 3 ","p");
  label->AddEntry(effNewPt," new Pt","p");
  // label->AddEntry(effBIn," BIn ","p");
  //label->AddEntry(effBOut," BOut pt 2 ","p");
  
  effDPhi->Draw("AP");
  effWeightPt->Draw("sameP");
  effPtMin->Draw("sameP");
  effPtMin2->Draw("sameP");
  effPtMin3->Draw("sameP");
  effNewPt->Draw("same");
  label->Draw("same");
  lineCut->Draw("same");  
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
