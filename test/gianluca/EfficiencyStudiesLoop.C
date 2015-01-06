#include <stdlib.h>
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "TF1.h"

//------------------------------------------------------------------------------------------------------------------------------------//
//Plot a 2d histogram of plateau efficiency values for different combination of parameters of minpt2 algorithm
//Author: Gian Luca Pinna
//-------------------------------------------------------------------------------------------------------------------------------------//


typedef struct {Float_t PtTrue,PtBIn,PtBOut,PtDPhi,PtGMT,Wheel_,Qual1_,Qual2_;} PTs_;

void EfficiencyStudiesLoop(string FileIn, Float_t dttfPtCut =25., Int_t Wheel = 0 ){
  
  TFile *f = new TFile(FileIn.c_str());
  TTree *t1 = (TTree*)f->Get("L1ITMuPtPlotter/ElderTree");
  
  
  TH2F *hValues= new TH2F("hValues"," ", 55,10.,120.,55,10.,120.);  
  
  PTs_ pts;
  t1->SetBranchAddress("Pts",&pts.PtTrue);
  TF1 *f1 = new TF1("f1","pol1",dttfPtCut+25,130);
  
  for (Int_t l=10; l<=120; l+=2) {
    for (Int_t m=10; m<=120; m+=2) {
      
      
      TEfficiency *effPtMin2= new TEfficiency("effPtMin2"," ", 60,-1.,119.);
      
      Int_t nentries = (Int_t)t1->GetEntries();
      Int_t nevents = 0;
      for (Int_t i=0; i<nentries; i++) {
	
	t1->GetEntry(i);
	
	if(Wheel==0){
	  if(abs(pts.Wheel_)==2) continue; 
	}
	else if(pts.Wheel_!=Wheel) continue;  
	nevents++; 
	
	
	Float_t MinPt2=pts.PtDPhi;
	
	if(pts.Qual1_==1 &&  pts.Qual2_==1){
	  if(fabs(pts.PtDPhi-pts.PtBIn)<l && (pts.PtDPhi-pts.PtBOut)>m){
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
	effPtMin2->Fill(MinPt2>=dttfPtCut,pts.PtTrue);
      }
      
      effPtMin2->Fit(f1,"R");
      Float_t Value = f1->GetParameter(0);
      std::cout<<"l "<<l<<" m "<<m<<" Value "<<Value<<std::endl;
      hValues->Fill(l,m,Value);
      effPtMin2->Delete();
    }  
 }
  hValues->Draw("colz");
}


