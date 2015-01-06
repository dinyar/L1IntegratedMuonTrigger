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
#include "TH2F.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TLine.h"
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------//
//Plot a 2d histogram of ratio of minPt2 rate and delta phi rate at pt =20 for different combination of parameters of minpt2 algorithm 
//Author: Gian Luca Pinna
//-------------------------------------------------------------------------------------------------------------------------------------//


typedef struct {Float_t PtBIn,PtBOut,PtDPhi,GMTPt_,Wheel_,Qual1_,Qual2_;} PTs_;

void RateStudiesLoop(string FileIn){

  TFile *f = new TFile(FileIn.c_str());
  TTree *t1 = (TTree*)f->Get("L1ITMuPtPlotter/ElderTree");

  PTs_ pts;
  t1->SetBranchAddress("Pts",&pts.PtBIn);
 
  TH2F *hValues= new TH2F("hValues"," ", 55,10.,120.,55,10.,120.);  

  for (Int_t l=10; l<120; l+=2) {
    for (Int_t m=10; m<120; m+=2) {

  TH1F *hRateDPhi= new TH1F("hRateBDPhi"," ", 60,-1.,119.);
  TH1F *hRatePtMin2= new TH1F("hRatePtMin2"," ", 60,-1.,119.);
  
Int_t nentries = (Int_t)t1->GetEntries();
  Int_t nevents = 0;
  for (Int_t i=0; i<nentries; i++) {
 
   
    t1->GetEntry(i);
    if(abs(pts.Wheel_)==3) continue; 
    nevents++; 
        
    Float_t MinPt2=pts.PtDPhi;;
    
    if(pts.Qual1_==1 &&  pts.Qual2_==1){
      if(fabs(pts.PtDPhi-pts.PtBIn)<80 && (pts.PtDPhi-pts.PtBOut)>10){
	if(pts.PtBIn<MinPt2) MinPt2=pts.PtBIn;
      }
      
      else if(abs(pts.PtDPhi-pts.PtBOut)<l &&(pts.PtDPhi-pts.PtBIn)>m ){
	if(pts.PtBOut<MinPt2) MinPt2=pts.PtBOut;
      }
    }
 
    else if(pts.Qual1_==1&&pts.Qual2_!=1){
      if(pts.PtBIn<MinPt2) MinPt2=pts.PtBIn;
    }
    
    else if(pts.Qual1_!=1&&pts.Qual2_==1){
      if(pts.PtBOut<MinPt2) MinPt2=pts.PtBOut;
    }
    

    for(int j = 0.; j<=120; j+=2){
      
      if(pts.PtDPhi>=j)  hRateDPhi->Fill(j);
      if(MinPt2>=j)  hRatePtMin2->Fill(j);
      
    }
  }
  
  Float_t Value=1-hRatePtMin2->GetBinContent(12)/hRateDPhi->GetBinContent(12);
  hValues->Fill(l,m,Value);
  hRatePtMin2->Delete();
  hRateDPhi->Delete();

    }

  }
  hValues->SetXTitle("l");
  hValues->SetYTitle("m");
  hValues->Draw("colz");
}


