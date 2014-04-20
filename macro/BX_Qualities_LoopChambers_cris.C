#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <iostream>
using namespace std;

void BX_Qualities_LoopChambers_cris()
{

gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);

//creating the Canvas
TCanvas *c1 = new TCanvas("c1", "Data", 1000, 600);
//c1->cd(1);

//0 all qualities, 1 Q46_old - Q56_new, 2 Q03_old - Q04_new

//FIRST FILE************************+
TFile f("L1ITMuonBarrelTimingEffPlots_17aprilT0.root");

//SECOND FILE************************+
TFile fsp("L1ITMuonBarrelTimingEffPlots_SP_17aprilT0.root");


TH1F *MB1s_bx_All=new TH1F("mb1All","BX distribution for MB1 All Qualities",16,-8,8);
TH1F *MB1s_bx_High=new TH1F("mb1High","BX distribution for MB1 High Qualities",16,-8,8);
TH1F *MB1s_bx_Low=new TH1F("mb1Low","BX distribution for MB1 Low Qualities",16,-8,8);

TH1F *MB2s_bx_All=new TH1F("mb2All","BX distribution for MB2 All Qualities",16,-8,8);
TH1F *MB2s_bx_High=new TH1F("mb2High","BX distribution for MB2 High Qualities",16,-8,8);
TH1F *MB2s_bx_Low=new TH1F("mb2Low","BX distribution for MB2 Low Qualities",16,-8,8);

TH1F *MB3s_bx_All=new TH1F("mb3All","BX distribution for MB3 All Qualities",16,-8,8);
TH1F *MB3s_bx_High=new TH1F("mb3High","BX distribution for MB3 High Qualities",16,-8,8);
TH1F *MB3s_bx_Low=new TH1F("mb3Low","BX distribution for MB3 Low Qualities",16,-8,8);

TH1F *MB4s_bx_All=new TH1F("mb4All","BX distribution for MB4 All Qualities",16,-8,8);
TH1F *MB4s_bx_High=new TH1F("mb4High","BX distribution for MB4 High Qualities",16,-8,8);
TH1F *MB4s_bx_Low=new TH1F("mb4Low","BX distribution for MB4 Low Qualities",16,-8,8);

TH1F *ALL_bx=new TH1F("ALL","BX distribution for ALL chambers ALL qualities",16,-8,8);
TH1F *ALL_bx_Qhigh=new TH1F("ALLHigh","BX distribution for ALL chambers High Qualities",16,-8,8);
TH1F *ALL_bx_Qlow=new TH1F("ALLLow","BX distribution for ALL chambers Low Qualities",16,-8,8);

TH1F *ALL_bx_Q0=new TH1F("ALL_Q0","BX distribution for ALL chambers Q=0",16,-8,8);
TH1F *ALL_bx_Q1=new TH1F("ALL_Q1","BX distribution for ALL chambers Q=1",16,-8,8);
TH1F *ALL_bx_Q2=new TH1F("ALL_Q2","BX distribution for ALL chambers Q=2",16,-8,8);
TH1F *ALL_bx_Q3=new TH1F("ALL_Q3","BX distribution for ALL chambers Q=3",16,-8,8);
TH1F *ALL_bx_Q4=new TH1F("ALL_Q4","BX distribution for ALL chambers Q=4",16,-8,8);
TH1F *ALL_bx_Q5=new TH1F("ALL_Q5","BX distribution for ALL chambers Q=5",16,-8,8);
TH1F *ALL_bx_Q6=new TH1F("ALL_Q6","BX distribution for ALL chambers Q=6",16,-8,8);

TString name0, name1, name2, name3, name4, name5, name6;
stringstream wheelstr;
stringstream sectorstr;
stringstream chamberstr;

TText *thlow=new TText(-5, 10, "LowQ_SP=[0,4]-LowQ_TP=[0,3]");
thlow->SetTextAlign(11);
//thlow->SetTextSize(0.12);

TText *thhi=new TText(-5, 10, "HighQ_SP=[5,6]-HighQ_TP=[4,6]");
thhi->SetTextAlign(11);
//thhi->SetTextSize(0.12);


for(int k=1;k<5;k++){ //loop chamber

for(int i=-2;i<3;i++){ //loop wheels
  for(int j=1;j<15;j++){ //loop sector
  	  
if (j>12 && k<4){
}else{
	
	wheelstr << i;
	sectorstr << j;
	chamberstr << k;
	
	name0="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q0"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh0=(TH1F*)f.Get(name0.Data());

	name1="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q1"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh1=(TH1F*)f.Get(name1.Data());
	
	name2="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q2"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh2=(TH1F*)f.Get(name2.Data());

	name3="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q3"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh3=(TH1F*)f.Get(name3.Data());

	name4="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q4"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh4=(TH1F*)f.Get(name4.Data());

	name5="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q5"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh5=(TH1F*)f.Get(name5.Data());

	name6="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q6"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh6=(TH1F*)f.Get(name6.Data());
	
	name0.Clear();
	name1.Clear();
	name2.Clear();
	name3.Clear();
	name4.Clear();
	name5.Clear();
	name6.Clear();
	wheelstr.str("");
	sectorstr.str("");
	chamberstr.str("");

	if (k==1){
	  MB1s_bx_All->Add(MB1s_bx_All,hh0);
	  MB1s_bx_All->Add(MB1s_bx_All,hh1);
	  MB1s_bx_All->Add(MB1s_bx_All,hh2);
	  MB1s_bx_All->Add(MB1s_bx_All,hh3);
	  MB1s_bx_All->Add(MB1s_bx_All,hh4);
	  MB1s_bx_All->Add(MB1s_bx_All,hh5);
	  MB1s_bx_All->Add(MB1s_bx_All,hh6);
	  
	  MB1s_bx_High->Add(MB1s_bx_High,hh4); //this not for SP
	  MB1s_bx_High->Add(MB1s_bx_High,hh5);
	  MB1s_bx_High->Add(MB1s_bx_High,hh6);
	  
	  MB1s_bx_Low->Add(MB1s_bx_Low,hh0);
	  MB1s_bx_Low->Add(MB1s_bx_Low,hh1);
	  MB1s_bx_Low->Add(MB1s_bx_Low,hh2);
	  MB1s_bx_Low->Add(MB1s_bx_Low,hh3);
	};
	if (k==2){
	  MB2s_bx_All->Add(MB2s_bx_All,hh0);
	  MB2s_bx_All->Add(MB2s_bx_All,hh1);
	  MB2s_bx_All->Add(MB2s_bx_All,hh2);
	  MB2s_bx_All->Add(MB2s_bx_All,hh3);
	  MB2s_bx_All->Add(MB2s_bx_All,hh4);
	  MB2s_bx_All->Add(MB2s_bx_All,hh5);
	  MB2s_bx_All->Add(MB2s_bx_All,hh6);
	  
	  MB2s_bx_High->Add(MB2s_bx_High,hh4); //this not for SP
	  MB2s_bx_High->Add(MB2s_bx_High,hh5);
	  MB2s_bx_High->Add(MB2s_bx_High,hh6);
	  
	  MB2s_bx_Low->Add(MB2s_bx_Low,hh0);
	  MB2s_bx_Low->Add(MB2s_bx_Low,hh1);
	  MB2s_bx_Low->Add(MB2s_bx_Low,hh2);
	  MB2s_bx_Low->Add(MB2s_bx_Low,hh3);
	};
	if (k==3){
	  MB3s_bx_All->Add(MB3s_bx_All,hh0);
	  MB3s_bx_All->Add(MB3s_bx_All,hh1);
	  MB3s_bx_All->Add(MB3s_bx_All,hh2);
	  MB3s_bx_All->Add(MB3s_bx_All,hh3);
	  MB3s_bx_All->Add(MB3s_bx_All,hh4);
	  MB3s_bx_All->Add(MB3s_bx_All,hh5);
	  MB3s_bx_All->Add(MB3s_bx_All,hh6);
	  
	  MB3s_bx_High->Add(MB3s_bx_High,hh4); //this not for SP
	  MB3s_bx_High->Add(MB3s_bx_High,hh5);
	  MB3s_bx_High->Add(MB3s_bx_High,hh6);
	  
	  MB3s_bx_Low->Add(MB3s_bx_Low,hh0);
	  MB3s_bx_Low->Add(MB3s_bx_Low,hh1);
	  MB3s_bx_Low->Add(MB3s_bx_Low,hh2);
	  MB3s_bx_Low->Add(MB3s_bx_Low,hh3);
	};
	if (k==4){
	  MB4s_bx_All->Add(MB4s_bx_All,hh0);
	  MB4s_bx_All->Add(MB4s_bx_All,hh1);
	  MB4s_bx_All->Add(MB4s_bx_All,hh2);
	  MB4s_bx_All->Add(MB4s_bx_All,hh3);
	  MB4s_bx_All->Add(MB4s_bx_All,hh4);
	  MB4s_bx_All->Add(MB4s_bx_All,hh5);
	  MB4s_bx_All->Add(MB4s_bx_All,hh6);
	  
	  MB4s_bx_High->Add(MB4s_bx_High,hh4); //this not for SP
	  MB4s_bx_High->Add(MB4s_bx_High,hh5);
	  MB4s_bx_High->Add(MB4s_bx_High,hh6);
	  
	  MB4s_bx_Low->Add(MB4s_bx_Low,hh0);
	  MB4s_bx_Low->Add(MB4s_bx_Low,hh1);
	  MB4s_bx_Low->Add(MB4s_bx_Low,hh2);
	  MB4s_bx_Low->Add(MB4s_bx_Low,hh3);
	};
	
	  ALL_bx->Add(ALL_bx,hh0);
	  ALL_bx->Add(ALL_bx,hh1);
	  ALL_bx->Add(ALL_bx,hh2);
	  ALL_bx->Add(ALL_bx,hh3);
	  ALL_bx->Add(ALL_bx,hh4);
	  ALL_bx->Add(ALL_bx,hh5);
	  ALL_bx->Add(ALL_bx,hh6);

  	  ALL_bx_Qhigh->Add(ALL_bx_Qhigh,hh4);//this not for SP
  	  ALL_bx_Qhigh->Add(ALL_bx_Qhigh,hh5);
  	  ALL_bx_Qhigh->Add(ALL_bx_Qhigh,hh6);

  	  ALL_bx_Qlow->Add(ALL_bx_Qlow,hh0);
  	  ALL_bx_Qlow->Add(ALL_bx_Qlow,hh1);
  	  ALL_bx_Qlow->Add(ALL_bx_Qlow,hh2);
  	  ALL_bx_Qlow->Add(ALL_bx_Qlow,hh3);

  	  ALL_bx_Q0->Add(ALL_bx_Q0,hh0);
  	  ALL_bx_Q1->Add(ALL_bx_Q1,hh1);
  	  ALL_bx_Q2->Add(ALL_bx_Q2,hh2);
  	  ALL_bx_Q3->Add(ALL_bx_Q3,hh3);
  	  ALL_bx_Q4->Add(ALL_bx_Q4,hh4);
  	  ALL_bx_Q5->Add(ALL_bx_Q5,hh5);
  	  ALL_bx_Q6->Add(ALL_bx_Q6,hh6);
    };
};//sector
};//wheel

};//chamber

MB1s_bx_All->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_All->GetXaxis()->CenterLabels();
MB1s_bx_All->GetXaxis()->SetTitle("BX");
MB1s_bx_All->GetYaxis()->SetTitle("# entries");
MB1s_bx_All->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_All->GetYaxis()->CenterTitle(1);
MB1s_bx_All->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_All->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_All->SetLineColor(kBlue);
MB1s_bx_All->SetLineStyle(1);
MB1s_bx_All->Draw();
c1->Update();
c1->SetLogy();

MB1s_bx_High->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_High->GetXaxis()->CenterLabels();
MB1s_bx_High->GetXaxis()->SetTitle("BX");
MB1s_bx_High->GetYaxis()->SetTitle("# entries");
MB1s_bx_High->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_High->GetYaxis()->CenterTitle(1);
MB1s_bx_High->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_High->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_High->SetLineColor(kBlue);
MB1s_bx_High->SetLineStyle(1);
MB1s_bx_High->Draw();
c1->Update();
c1->SetLogy();

MB1s_bx_Low->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_Low->GetXaxis()->CenterLabels();
MB1s_bx_Low->GetXaxis()->SetTitle("BX");
MB1s_bx_Low->GetYaxis()->SetTitle("# entries");
MB1s_bx_Low->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_Low->GetYaxis()->CenterTitle(1);
MB1s_bx_Low->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_Low->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_Low->SetLineColor(kBlue);
MB1s_bx_Low->SetLineStyle(1);
MB1s_bx_Low->Draw();
c1->Update();
c1->SetLogy();


MB2s_bx_All->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_All->GetXaxis()->CenterLabels();
MB2s_bx_All->GetXaxis()->SetTitle("BX");
MB2s_bx_All->GetYaxis()->SetTitle("# entries");
MB2s_bx_All->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_All->GetYaxis()->CenterTitle(1);
MB2s_bx_All->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_All->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_All->SetLineColor(kBlue);
MB2s_bx_All->SetLineStyle(1);
MB2s_bx_All->Draw();
c1->Update();
c1->SetLogy();

MB2s_bx_High->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_High->GetXaxis()->CenterLabels();
MB2s_bx_High->GetXaxis()->SetTitle("BX");
MB2s_bx_High->GetYaxis()->SetTitle("# entries");
MB2s_bx_High->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_High->GetYaxis()->CenterTitle(1);
MB2s_bx_High->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_High->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_High->SetLineColor(kBlue);
MB2s_bx_High->SetLineStyle(1);
MB2s_bx_High->Draw();
c1->Update();
c1->SetLogy();

MB2s_bx_Low->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_Low->GetXaxis()->CenterLabels();
MB2s_bx_Low->GetXaxis()->SetTitle("BX");
MB2s_bx_Low->GetYaxis()->SetTitle("# entries");
MB2s_bx_Low->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_Low->GetYaxis()->CenterTitle(1);
MB2s_bx_Low->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_Low->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_Low->SetLineColor(kBlue);
MB2s_bx_Low->SetLineStyle(1);
MB2s_bx_Low->Draw();
c1->Update();
c1->SetLogy();

MB3s_bx_All->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_All->GetXaxis()->CenterLabels();
MB3s_bx_All->GetXaxis()->SetTitle("BX");
MB3s_bx_All->GetYaxis()->SetTitle("# entries");
MB3s_bx_All->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_All->GetYaxis()->CenterTitle(1);
MB3s_bx_All->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_All->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_All->SetLineColor(kBlue);
MB3s_bx_All->SetLineStyle(1);
MB3s_bx_All->Draw();
c1->Update();
c1->SetLogy();

MB3s_bx_High->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_High->GetXaxis()->CenterLabels();
MB3s_bx_High->GetXaxis()->SetTitle("BX");
MB3s_bx_High->GetYaxis()->SetTitle("# entries");
MB3s_bx_High->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_High->GetYaxis()->CenterTitle(1);
MB3s_bx_High->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_High->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_High->SetLineColor(kBlue);
MB3s_bx_High->SetLineStyle(1);
MB3s_bx_High->Draw();
c1->Update();
c1->SetLogy();

MB3s_bx_Low->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_Low->GetXaxis()->CenterLabels();
MB3s_bx_Low->GetXaxis()->SetTitle("BX");
MB3s_bx_Low->GetYaxis()->SetTitle("# entries");
MB3s_bx_Low->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_Low->GetYaxis()->CenterTitle(1);
MB3s_bx_Low->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_Low->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_Low->SetLineColor(kBlue);
MB3s_bx_Low->SetLineStyle(1);
MB3s_bx_Low->Draw();
c1->Update();
c1->SetLogy();

MB4s_bx_All->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_All->GetXaxis()->CenterLabels();
MB4s_bx_All->GetXaxis()->SetTitle("BX");
MB4s_bx_All->GetYaxis()->SetTitle("# entries");
MB4s_bx_All->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_All->GetYaxis()->CenterTitle(1);
MB4s_bx_All->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_All->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_All->SetLineColor(kBlue);
MB4s_bx_All->SetLineStyle(1);
MB4s_bx_All->Draw();
c1->Update();
c1->SetLogy();

MB4s_bx_High->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_High->GetXaxis()->CenterLabels();
MB4s_bx_High->GetXaxis()->SetTitle("BX");
MB4s_bx_High->GetYaxis()->SetTitle("# entries");
MB4s_bx_High->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_High->GetYaxis()->CenterTitle(1);
MB4s_bx_High->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_High->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_High->SetLineColor(kBlue);
MB4s_bx_High->SetLineStyle(1);
MB4s_bx_High->Draw();
c1->Update();
c1->SetLogy();

MB4s_bx_Low->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_Low->GetXaxis()->CenterLabels();
MB4s_bx_Low->GetXaxis()->SetTitle("BX");
MB4s_bx_Low->GetYaxis()->SetTitle("# entries");
MB4s_bx_Low->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_Low->GetYaxis()->CenterTitle(1);
MB4s_bx_Low->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_Low->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_Low->SetLineColor(kBlue);
MB4s_bx_Low->SetLineStyle(1);
MB4s_bx_Low->Draw();
c1->Update();
c1->SetLogy();



ALL_bx->GetXaxis()->SetRangeUser(-5,5);
ALL_bx->GetXaxis()->CenterLabels();
ALL_bx->GetXaxis()->SetTitle("BX");
ALL_bx->GetYaxis()->SetTitle("# entries");
ALL_bx->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx->GetYaxis()->CenterTitle(1);
ALL_bx->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx->GetXaxis()->SetTitleOffset(1.2);
ALL_bx->SetLineColor(kBlue);
ALL_bx->SetLineStyle(1);
ALL_bx->Draw();
c1->Update();
c1->SetLogy();


ALL_bx_Qhigh->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Qhigh->GetXaxis()->CenterLabels();
ALL_bx_Qhigh->GetXaxis()->SetTitle("BX");
ALL_bx_Qhigh->GetYaxis()->SetTitle("# entries");
ALL_bx_Qhigh->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Qhigh->GetYaxis()->CenterTitle(1);
ALL_bx_Qhigh->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Qhigh->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Qhigh->SetLineColor(kBlue);
ALL_bx_Qhigh->SetLineStyle(1);
ALL_bx_Qhigh->Draw();
c1->Update();
c1->SetLogy();


ALL_bx_Qlow->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Qlow->GetXaxis()->CenterLabels();
ALL_bx_Qlow->GetXaxis()->SetTitle("BX");
ALL_bx_Qlow->GetYaxis()->SetTitle("# entries");
ALL_bx_Qlow->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Qlow->GetYaxis()->CenterTitle(1);
ALL_bx_Qlow->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Qlow->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Qlow->SetLineColor(kBlue);
ALL_bx_Qlow->SetLineStyle(1);
ALL_bx_Qlow->Draw();
c1->Update();
c1->SetLogy();

ALL_bx_Q0->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q0->GetXaxis()->CenterLabels();
ALL_bx_Q0->GetXaxis()->SetTitle("BX");
ALL_bx_Q0->GetYaxis()->SetTitle("# entries");
ALL_bx_Q0->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q0->GetYaxis()->CenterTitle(1);
ALL_bx_Q0->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q0->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q0->SetLineColor(kBlue);
ALL_bx_Q0->SetLineStyle(1);
ALL_bx_Q0->Draw();
c1->Update();
c1->SetLogy();

ALL_bx_Q1->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q1->GetXaxis()->CenterLabels();
ALL_bx_Q1->GetXaxis()->SetTitle("BX");
ALL_bx_Q1->GetYaxis()->SetTitle("# entries");
ALL_bx_Q1->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q1->GetYaxis()->CenterTitle(1);
ALL_bx_Q1->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q1->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q1->SetLineColor(kBlue);
ALL_bx_Q1->SetLineStyle(1);
ALL_bx_Q1->Draw();
c1->Update();
c1->SetLogy();


ALL_bx_Q2->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q2->GetXaxis()->CenterLabels();
ALL_bx_Q2->GetXaxis()->SetTitle("BX");
ALL_bx_Q2->GetYaxis()->SetTitle("# entries");
ALL_bx_Q2->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q2->GetYaxis()->CenterTitle(1);
ALL_bx_Q2->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q2->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q2->SetLineColor(kBlue);
ALL_bx_Q2->SetLineStyle(1);
ALL_bx_Q2->Draw();
c1->Update();
c1->SetLogy();


ALL_bx_Q3->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q3->GetXaxis()->CenterLabels();
ALL_bx_Q3->GetXaxis()->SetTitle("BX");
ALL_bx_Q3->GetYaxis()->SetTitle("# entries");
ALL_bx_Q3->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q3->GetYaxis()->CenterTitle(1);
ALL_bx_Q3->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q3->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q3->SetLineColor(kBlue);
ALL_bx_Q3->SetLineStyle(1);
ALL_bx_Q3->Draw();
c1->Update();
c1->SetLogy();


ALL_bx_Q4->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q4->GetXaxis()->CenterLabels();
ALL_bx_Q4->GetXaxis()->SetTitle("BX");
ALL_bx_Q4->GetYaxis()->SetTitle("# entries");
ALL_bx_Q4->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q4->GetYaxis()->CenterTitle(1);
ALL_bx_Q4->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q4->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q4->SetLineColor(kBlue);
ALL_bx_Q4->SetLineStyle(1);
ALL_bx_Q4->Draw();
c1->Update();
c1->SetLogy();

ALL_bx_Q5->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q5->GetXaxis()->CenterLabels();
ALL_bx_Q5->GetXaxis()->SetTitle("BX");
ALL_bx_Q5->GetYaxis()->SetTitle("# entries");
ALL_bx_Q5->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q5->GetYaxis()->CenterTitle(1);
ALL_bx_Q5->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q5->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q5->SetLineColor(kBlue);
ALL_bx_Q5->SetLineStyle(1);
ALL_bx_Q5->Draw();
c1->Update();
c1->SetLogy();

ALL_bx_Q6->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q6->GetXaxis()->CenterLabels();
ALL_bx_Q6->GetXaxis()->SetTitle("BX");
ALL_bx_Q6->GetYaxis()->SetTitle("# entries");
ALL_bx_Q6->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q6->GetYaxis()->CenterTitle(1);
ALL_bx_Q6->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q6->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q6->SetLineColor(kBlue);
ALL_bx_Q6->SetLineStyle(1);
ALL_bx_Q6->Draw();
c1->Update();
c1->SetLogy();

TH1F *MB1s_bx_All_SP=new TH1F("mb1All_SP","BX distribution for MB1 All Qualities",16,-8,8);
TH1F *MB1s_bx_High_SP=new TH1F("mb1High_SP","BX distribution for MB1 High Qualities",16,-8,8);
TH1F *MB1s_bx_Low_SP=new TH1F("mb1Low_SP","BX distribution for MB1 Low Qualities",16,-8,8);

TH1F *MB2s_bx_All_SP=new TH1F("mb2All_SP","BX distribution for MB2 All Qualities",16,-8,8);
TH1F *MB2s_bx_High_SP=new TH1F("mb2High_SP","BX distribution for MB2 High Qualities",16,-8,8);
TH1F *MB2s_bx_Low_SP=new TH1F("mb2Low_SP","BX distribution for MB2 Low Qualities",16,-8,8);

TH1F *MB3s_bx_All_SP=new TH1F("mb3All_SP","BX distribution for MB3 All Qualities",16,-8,8);
TH1F *MB3s_bx_High_SP=new TH1F("mb3High_SP","BX distribution for MB3 High Qualities",16,-8,8);
TH1F *MB3s_bx_Low_SP=new TH1F("mb3Low_SP","BX distribution for MB3 Low Qualities",16,-8,8);

TH1F *MB4s_bx_All_SP=new TH1F("mb4All_SP","BX distribution for MB4 All Qualities",16,-8,8);
TH1F *MB4s_bx_High_SP=new TH1F("mb4High_SP","BX distribution for MB4 High Qualities",16,-8,8);
TH1F *MB4s_bx_Low_SP=new TH1F("mb4Low_SP","BX distribution for MB4 Low Qualities",16,-8,8);

TH1F *ALL_bx_SP=new TH1F("ALL_SP","BX distribution for ALL chambers ALL qualities",16,-8,8);
TH1F *ALL_bx_Qhigh_SP=new TH1F("ALLHigh_SP","BX distribution for ALL chambers High Qualities",16,-8,8);
TH1F *ALL_bx_Qlow_SP=new TH1F("ALLLow_SP","BX distribution for ALL chambers Low Qualities",16,-8,8);

TH1F *ALL_bx_Q0_SP=new TH1F("ALL_Q0_SP","BX distribution for ALL chambers Q=0",16,-8,8);
TH1F *ALL_bx_Q1_SP=new TH1F("ALL_Q1_SP","BX distribution for ALL chambers Q=1",16,-8,8);
TH1F *ALL_bx_Q2_SP=new TH1F("ALL_Q2_SP","BX distribution for ALL chambers Q=2",16,-8,8);
TH1F *ALL_bx_Q3_SP=new TH1F("ALL_Q3_SP","BX distribution for ALL chambers Q=3",16,-8,8);
TH1F *ALL_bx_Q4_SP=new TH1F("ALL_Q4_SP","BX distribution for ALL chambers Q=4",16,-8,8);
TH1F *ALL_bx_Q5_SP=new TH1F("ALL_Q5_SP","BX distribution for ALL chambers Q=5",16,-8,8);
TH1F *ALL_bx_Q6_SP=new TH1F("ALL_Q6_SP","BX distribution for ALL chambers Q=6",16,-8,8);

for(int k=1;k<5;k++){ //loop chamber

for(int i=-2;i<3;i++){ //loop wheels
  for(int j=1;j<15;j++){ //loop sector
  	  
if (j>12 && k<4){
}else{
	
	wheelstr << i;
	sectorstr << j;
	chamberstr << k;
	
	name0="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q0"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh0_sp=(TH1F*)fsp.Get(name0.Data());

	name1="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q1"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh1_sp=(TH1F*)fsp.Get(name1.Data());
	
	name2="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q2"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh2_sp=(TH1F*)fsp.Get(name2.Data());

	name3="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q3"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh3_sp=(TH1F*)fsp.Get(name3.Data());

	name4="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q4"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh4_sp=(TH1F*)fsp.Get(name4.Data());

	name5="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q5"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh5_sp=(TH1F*)fsp.Get(name5.Data());

	name6="L1ITMuTimingPlotter/Wheel"+ wheelstr.str()+"/Sector"+sectorstr.str()+"/Station"+chamberstr.str()+"/Timing_DTprim_BX_Q6"+"_W"+wheelstr.str()+"_Sec"+sectorstr.str()+"_St"+chamberstr.str();
	TH1F *hh6_sp=(TH1F*)fsp.Get(name6.Data());
	
	name0.Clear();
	name1.Clear();
	name2.Clear();
	name3.Clear();
	name4.Clear();
	name5.Clear();
	name6.Clear();
	wheelstr.str("");
	sectorstr.str("");
	chamberstr.str("");

	if (k==1){
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh0_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh1_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh2_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh3_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh4_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh5_sp);
	  MB1s_bx_All_SP->Add(MB1s_bx_All_SP,hh6_sp);
	  
	  MB1s_bx_High_SP->Add(MB1s_bx_High_SP,hh5_sp);
	  MB1s_bx_High_SP->Add(MB1s_bx_High_SP,hh6_sp);
	  
	  MB1s_bx_Low_SP->Add(MB1s_bx_Low_SP,hh0_sp);
	  MB1s_bx_Low_SP->Add(MB1s_bx_Low_SP,hh1_sp);
	  MB1s_bx_Low_SP->Add(MB1s_bx_Low_SP,hh2_sp);
	  MB1s_bx_Low_SP->Add(MB1s_bx_Low_SP,hh3_sp);
	  MB1s_bx_Low_SP->Add(MB1s_bx_Low_SP,hh4_sp);
	};
	if (k==2){
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh0_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh1_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh2_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh3_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh4_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh5_sp);
	  MB2s_bx_All_SP->Add(MB2s_bx_All_SP,hh6_sp);
	  
	  MB2s_bx_High_SP->Add(MB2s_bx_High_SP,hh5_sp);
	  MB2s_bx_High_SP->Add(MB2s_bx_High_SP,hh6_sp);
	  
	  MB2s_bx_Low_SP->Add(MB2s_bx_Low_SP,hh0_sp);
	  MB2s_bx_Low_SP->Add(MB2s_bx_Low_SP,hh1_sp);
	  MB2s_bx_Low_SP->Add(MB2s_bx_Low_SP,hh2_sp);
	  MB2s_bx_Low_SP->Add(MB2s_bx_Low_SP,hh3_sp);
	  MB2s_bx_Low_SP->Add(MB2s_bx_Low_SP,hh4_sp);
	};
	if (k==3){
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh0_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh1_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh2_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh3_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh4_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh5_sp);
	  MB3s_bx_All_SP->Add(MB3s_bx_All_SP,hh6_sp);
	  
	  MB3s_bx_High_SP->Add(MB3s_bx_High_SP,hh5_sp);
	  MB3s_bx_High_SP->Add(MB3s_bx_High_SP,hh6_sp);
	  
	  MB3s_bx_Low_SP->Add(MB3s_bx_Low_SP,hh0_sp);
	  MB3s_bx_Low_SP->Add(MB3s_bx_Low_SP,hh1_sp);
	  MB3s_bx_Low_SP->Add(MB3s_bx_Low_SP,hh2_sp);
	  MB3s_bx_Low_SP->Add(MB3s_bx_Low_SP,hh3_sp);
	  MB3s_bx_Low_SP->Add(MB3s_bx_Low_SP,hh4_sp);
	};
	if (k==4){
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh0_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh1_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh2_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh3_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh4_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh5_sp);
	  MB4s_bx_All_SP->Add(MB4s_bx_All_SP,hh6_sp);
	  
	  MB4s_bx_High_SP->Add(MB4s_bx_High_SP,hh5_sp);
	  MB4s_bx_High_SP->Add(MB4s_bx_High_SP,hh6_sp);
	  
	  MB4s_bx_Low_SP->Add(MB4s_bx_Low_SP,hh0_sp);
	  MB4s_bx_Low_SP->Add(MB4s_bx_Low_SP,hh1_sp);
	  MB4s_bx_Low_SP->Add(MB4s_bx_Low_SP,hh2_sp);
	  MB4s_bx_Low_SP->Add(MB4s_bx_Low_SP,hh3_sp);
	  MB4s_bx_Low_SP->Add(MB4s_bx_Low_SP,hh4_sp);
	};
	
	  ALL_bx_SP->Add(ALL_bx_SP,hh0_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh1_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh2_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh3_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh4_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh5_sp);
	  ALL_bx_SP->Add(ALL_bx_SP,hh6_sp);

  	  ALL_bx_Qhigh_SP->Add(ALL_bx_Qhigh_SP,hh5_sp);
  	  ALL_bx_Qhigh_SP->Add(ALL_bx_Qhigh_SP,hh6_sp);

  	  ALL_bx_Qlow_SP->Add(ALL_bx_Qlow_SP,hh0_sp);
  	  ALL_bx_Qlow_SP->Add(ALL_bx_Qlow_SP,hh1_sp);
  	  ALL_bx_Qlow_SP->Add(ALL_bx_Qlow_SP,hh2_sp);
  	  ALL_bx_Qlow_SP->Add(ALL_bx_Qlow_SP,hh3_sp);
  	  ALL_bx_Qlow_SP->Add(ALL_bx_Qlow_SP,hh4_sp);

  	  ALL_bx_Q0_SP->Add(ALL_bx_Q0_SP,hh0_sp);
  	  ALL_bx_Q1_SP->Add(ALL_bx_Q1_SP,hh1_sp);
  	  ALL_bx_Q2_SP->Add(ALL_bx_Q2_SP,hh2_sp);
  	  ALL_bx_Q3_SP->Add(ALL_bx_Q3_SP,hh3_sp);
  	  ALL_bx_Q4_SP->Add(ALL_bx_Q4_SP,hh4_sp);
  	  ALL_bx_Q5_SP->Add(ALL_bx_Q5_SP,hh5_sp);
  	  ALL_bx_Q6_SP->Add(ALL_bx_Q6_SP,hh6_sp);
    };
};//sector
};//wheel

};//chamber

MB1s_bx_All_SP->SetTitle("MB1 BX distrib. of DT_TP vs SP for all Q");
MB1s_bx_High_SP->SetTitle("MB1 BX distrib. of DT_TP vs SP for High Qualities");
MB1s_bx_Low_SP->SetTitle("MB1 BX distrib. of DT_TP vs with SP for Low Qualities");

MB2s_bx_All_SP->SetTitle("MB2 BX distrib. of DT_TP vs SP for all Q");
MB2s_bx_High_SP->SetTitle("MB2 BX distrib. of DT_TP vs SP for High Qualities");
MB2s_bx_Low_SP->SetTitle("MB2 BX distrib. of DT_TP vs with SP for Low Qualities");

MB3s_bx_All_SP->SetTitle("MB3 BX distrib. of DT_TP vs SP for all Q");
MB3s_bx_High_SP->SetTitle("MB3 BX distrib. of DT_TP vs SP for High Qualities");
MB3s_bx_Low_SP->SetTitle("MB3 BX distrib. of DT_TP vs with SP for Low Qualities");

MB4s_bx_All_SP->SetTitle("MB4 BX distrib. of DT_TP vs SP for all Q");
MB4s_bx_High_SP->SetTitle("MB4 BX distrib. of DT_TP vs SP for High Qualities");
MB4s_bx_Low_SP->SetTitle("MB4 BX distrib. of DT_TP vs with SP for Low Qualities");

ALL_bx_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for all Q");
ALL_bx_Qhigh_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for High Qualities");
ALL_bx_Qlow_SP->SetTitle("All Chambers BX distrib. of DT_TP vs with SP for Low Qualities");

ALL_bx_Q0_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=0");
ALL_bx_Q1_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=1");
ALL_bx_Q2_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=2");
ALL_bx_Q3_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=3");
ALL_bx_Q4_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=4");
ALL_bx_Q5_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=5");
ALL_bx_Q6_SP->SetTitle("All Chambers BX distrib. of DT_TP vs SP for Q=6");

c1->Clear();

MB1s_bx_All_SP->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_All_SP->GetXaxis()->CenterLabels();
MB1s_bx_All_SP->GetXaxis()->SetTitle("BX");
MB1s_bx_All_SP->GetYaxis()->SetTitle("# entries");
MB1s_bx_All_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_All_SP->GetYaxis()->CenterTitle(1);
MB1s_bx_All_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_All_SP->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_All_SP->SetLineWidth(3);
MB1s_bx_All_SP->SetLineColor(kRed);
MB1s_bx_All_SP->SetLineStyle(2);
MB1s_bx_All_SP->Draw();
MB1s_bx_All->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB1s_bx_All,"Blue= DT primitives");
leg->AddEntry(MB1s_bx_All_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("MB1_BXdist_SPvsLegacy_allQ.gif");
c1->Clear();


MB1s_bx_High_SP->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_High_SP->GetXaxis()->CenterLabels();
MB1s_bx_High_SP->GetXaxis()->SetTitle("BX");
MB1s_bx_High_SP->GetYaxis()->SetTitle("# entries");
MB1s_bx_High_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_High_SP->GetYaxis()->CenterTitle(1);
MB1s_bx_High_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_High_SP->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_High_SP->SetLineColor(kRed);
MB1s_bx_High_SP->SetLineWidth(3);
MB1s_bx_High_SP->SetLineStyle(2);
MB1s_bx_High_SP->Draw();
MB1s_bx_High->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB1s_bx_High,"Blue= DT primitives");
leg->AddEntry(MB1s_bx_High_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thhi->Draw();

c1->Update();
c1->SaveAs("MB1_BXdist_SPvsLegacy_HighQ.gif");
c1->Clear();


MB1s_bx_Low_SP->GetXaxis()->SetRangeUser(-5,5);
MB1s_bx_Low_SP->GetXaxis()->CenterLabels();
MB1s_bx_Low_SP->GetXaxis()->SetTitle("BX");
MB1s_bx_Low_SP->GetYaxis()->SetTitle("# entries");
MB1s_bx_Low_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB1s_bx_Low_SP->GetYaxis()->CenterTitle(1);
MB1s_bx_Low_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB1s_bx_Low_SP->GetXaxis()->SetTitleOffset(1.2);
MB1s_bx_Low_SP->SetLineColor(kRed);
MB1s_bx_Low_SP->SetLineWidth(3);
MB1s_bx_Low_SP->SetLineStyle(2);
MB1s_bx_Low_SP->Draw();
MB1s_bx_Low->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB1s_bx_Low,"Blue= DT primitives");
leg->AddEntry(MB1s_bx_Low_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thlow->Draw();

c1->Update();
c1->SaveAs("MB1_BXdist_SPvsLegacy_LowQ.gif");
c1->Clear();

MB2s_bx_All_SP->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_All_SP->GetXaxis()->CenterLabels();
MB2s_bx_All_SP->GetXaxis()->SetTitle("BX");
MB2s_bx_All_SP->GetYaxis()->SetTitle("# entries");
MB2s_bx_All_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_All_SP->GetYaxis()->CenterTitle(1);
MB2s_bx_All_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_All_SP->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_All_SP->SetLineColor(kRed);
MB2s_bx_All_SP->SetLineWidth(3);
MB2s_bx_All_SP->SetLineStyle(2);
MB2s_bx_All_SP->Draw();
MB2s_bx_All->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB2s_bx_All,"Blue= DT primitives");
leg->AddEntry(MB2s_bx_All_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("MB2_BXdist_SPvsLegacy_allQ.gif");
c1->Clear();


MB2s_bx_High_SP->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_High_SP->GetXaxis()->CenterLabels();
MB2s_bx_High_SP->GetXaxis()->SetTitle("BX");
MB2s_bx_High_SP->GetYaxis()->SetTitle("# entries");
MB2s_bx_High_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_High_SP->GetYaxis()->CenterTitle(1);
MB2s_bx_High_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_High_SP->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_High_SP->SetLineColor(kRed);
MB2s_bx_High_SP->SetLineWidth(3);
MB2s_bx_High_SP->SetLineStyle(2);
MB2s_bx_High_SP->Draw();
MB2s_bx_High->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB2s_bx_High,"Blue= DT primitives");
leg->AddEntry(MB2s_bx_High_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thhi->Draw();

c1->Update();
c1->SaveAs("MB2_BXdist_SPvsLegacy_HighQ.gif");
c1->Clear();

MB2s_bx_Low_SP->GetXaxis()->SetRangeUser(-5,5);
MB2s_bx_Low_SP->GetXaxis()->CenterLabels();
MB2s_bx_Low_SP->GetXaxis()->SetTitle("BX");
MB2s_bx_Low_SP->GetYaxis()->SetTitle("# entries");
MB2s_bx_Low_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB2s_bx_Low_SP->GetYaxis()->CenterTitle(1);
MB2s_bx_Low_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB2s_bx_Low_SP->GetXaxis()->SetTitleOffset(1.2);
MB2s_bx_Low_SP->SetLineColor(kRed);
MB2s_bx_Low_SP->SetLineWidth(3);
MB2s_bx_Low_SP->SetLineStyle(2);
MB2s_bx_Low_SP->Draw();
MB2s_bx_Low->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB2s_bx_Low,"Blue= DT primitives");
leg->AddEntry(MB2s_bx_Low_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thlow->Draw();

c1->Update();
c1->SaveAs("MB2_BXdist_SPvsLegacy_LowQ.gif");
c1->Clear();

MB3s_bx_All_SP->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_All_SP->GetXaxis()->CenterLabels();
MB3s_bx_All_SP->GetXaxis()->SetTitle("BX");
MB3s_bx_All_SP->GetYaxis()->SetTitle("# entries");
MB3s_bx_All_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_All_SP->GetYaxis()->CenterTitle(1);
MB3s_bx_All_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_All_SP->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_All_SP->SetLineColor(kRed);
MB3s_bx_All_SP->SetLineWidth(3);
MB3s_bx_All_SP->SetLineStyle(2);
MB3s_bx_All_SP->Draw();
MB3s_bx_All->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB3s_bx_All,"Blue= DT primitives");
leg->AddEntry(MB3s_bx_All_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("MB3_BXdist_SPvsLegacy_allQ.gif");
c1->Clear();

MB3s_bx_High_SP->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_High_SP->GetXaxis()->CenterLabels();
MB3s_bx_High_SP->GetXaxis()->SetTitle("BX");
MB3s_bx_High_SP->GetYaxis()->SetTitle("# entries");
MB3s_bx_High_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_High_SP->GetYaxis()->CenterTitle(1);
MB3s_bx_High_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_High_SP->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_High_SP->SetLineColor(kRed);
MB3s_bx_High_SP->SetLineWidth(3);
MB3s_bx_High_SP->SetLineStyle(2);
MB3s_bx_High_SP->Draw();
MB3s_bx_High->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB3s_bx_High,"Blue= DT primitives");
leg->AddEntry(MB3s_bx_High_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thhi->Draw();

c1->Update();
c1->SaveAs("MB3_BXdist_SPvsLegacy_HighQ.gif");
c1->Clear();

MB3s_bx_Low_SP->GetXaxis()->SetRangeUser(-5,5);
MB3s_bx_Low_SP->GetXaxis()->CenterLabels();
MB3s_bx_Low_SP->GetXaxis()->SetTitle("BX");
MB3s_bx_Low_SP->GetYaxis()->SetTitle("# entries");
MB3s_bx_Low_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB3s_bx_Low_SP->GetYaxis()->CenterTitle(1);
MB3s_bx_Low_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB3s_bx_Low_SP->GetXaxis()->SetTitleOffset(1.2);
MB3s_bx_Low_SP->SetLineColor(kRed);
MB3s_bx_Low_SP->SetLineWidth(3);
MB3s_bx_Low_SP->SetLineStyle(2);
MB3s_bx_Low_SP->Draw();
MB3s_bx_Low->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB3s_bx_Low,"Blue= DT primitives");
leg->AddEntry(MB3s_bx_Low_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thlow->Draw();

c1->Update();
c1->SaveAs("MB3_BXdist_SPvsLegacy_LowQ.gif");
c1->Clear();

MB4s_bx_All_SP->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_All_SP->GetXaxis()->CenterLabels();
MB4s_bx_All_SP->GetXaxis()->SetTitle("BX");
MB4s_bx_All_SP->GetYaxis()->SetTitle("# entries");
MB4s_bx_All_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_All_SP->GetYaxis()->CenterTitle(1);
MB4s_bx_All_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_All_SP->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_All_SP->SetLineColor(kRed);
MB4s_bx_All_SP->SetLineWidth(3);
MB4s_bx_All_SP->SetLineStyle(2);
MB4s_bx_All_SP->Draw();
MB4s_bx_All->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB4s_bx_All,"Blue= DT primitives");
leg->AddEntry(MB4s_bx_All_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("MB4_BXdist_SPvsLegacy_allQ.gif");
c1->Clear();


MB4s_bx_High_SP->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_High_SP->GetXaxis()->CenterLabels();
MB4s_bx_High_SP->GetXaxis()->SetTitle("BX");
MB4s_bx_High_SP->GetYaxis()->SetTitle("# entries");
MB4s_bx_High_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_High_SP->GetYaxis()->CenterTitle(1);
MB4s_bx_High_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_High_SP->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_High_SP->SetLineColor(kRed);
MB4s_bx_High_SP->SetLineStyle(2);
MB4s_bx_High_SP->SetLineWidth(3);
MB4s_bx_High_SP->Draw();
MB4s_bx_High->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB4s_bx_High,"Blue= DT primitives");
leg->AddEntry(MB4s_bx_High_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thhi->Draw();

c1->Update();
c1->SaveAs("MB4_BXdist_SPvsLegacy_HighQ.gif");
c1->Clear();

MB4s_bx_Low_SP->GetXaxis()->SetRangeUser(-5,5);
MB4s_bx_Low_SP->GetXaxis()->CenterLabels();
MB4s_bx_Low_SP->GetXaxis()->SetTitle("BX");
MB4s_bx_Low_SP->GetYaxis()->SetTitle("# entries");
MB4s_bx_Low_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
MB4s_bx_Low_SP->GetYaxis()->CenterTitle(1);
MB4s_bx_Low_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
MB4s_bx_Low_SP->GetXaxis()->SetTitleOffset(1.2);
MB4s_bx_Low_SP->SetLineColor(kRed);
MB4s_bx_Low_SP->SetLineStyle(2);
MB4s_bx_Low_SP->SetLineWidth(3);
MB4s_bx_Low_SP->Draw();
MB4s_bx_Low->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(MB4s_bx_Low,"Blue= DT primitives");
leg->AddEntry(MB4s_bx_Low_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thlow->Draw();

c1->Update();
c1->SaveAs("MB4_BXdist_SPvsLegacy_LowQ.gif");
c1->Clear();

ALL_bx_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_SP->GetXaxis()->CenterLabels();
ALL_bx_SP->GetXaxis()->SetTitle("BX");
ALL_bx_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_SP->GetYaxis()->CenterTitle(1);
ALL_bx_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_SP->SetLineColor(kRed);
ALL_bx_SP->SetLineWidth(3);
ALL_bx_SP->SetLineStyle(2);
ALL_bx_SP->Draw();
ALL_bx->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx,"Blue= DT primitives");
leg->AddEntry(ALL_bx_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_allQ.gif");
c1->Clear();

ALL_bx_Qhigh_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Qhigh_SP->GetXaxis()->CenterLabels();
ALL_bx_Qhigh_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Qhigh_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Qhigh_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Qhigh_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Qhigh_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Qhigh_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Qhigh_SP->SetLineColor(kRed);
ALL_bx_Qhigh_SP->SetLineWidth(3);
ALL_bx_Qhigh_SP->SetLineStyle(2);
ALL_bx_Qhigh_SP->Draw();
ALL_bx_Qhigh->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Qhigh,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Qhigh_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thhi->Draw();

c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_HighQ.gif");
c1->Clear();

ALL_bx_Qlow_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Qlow_SP->GetXaxis()->CenterLabels();
ALL_bx_Qlow_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Qlow_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Qlow_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Qlow_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Qlow_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Qlow_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Qlow_SP->SetLineColor(kRed);
ALL_bx_Qlow_SP->SetLineWidth(3);
ALL_bx_Qlow_SP->SetLineStyle(2);
ALL_bx_Qlow_SP->Draw();
ALL_bx_Qlow->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Qlow,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Qlow_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();

thlow->Draw();


c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_LowQ.gif");
c1->Clear();





ALL_bx_Q0_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q0_SP->GetXaxis()->CenterLabels();
ALL_bx_Q0_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q0_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q0_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q0_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q0_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q0_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q0_SP->SetLineColor(kRed);
ALL_bx_Q0_SP->SetLineWidth(3);
ALL_bx_Q0_SP->SetLineStyle(2);
ALL_bx_Q0_SP->Draw();
ALL_bx_Q0->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q0,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q0_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q0.gif");
c1->Clear();


ALL_bx_Q1_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q1_SP->GetXaxis()->CenterLabels();
ALL_bx_Q1_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q1_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q1_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q1_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q1_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q1_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q1_SP->SetLineColor(kRed);
ALL_bx_Q1_SP->SetLineWidth(3);
ALL_bx_Q1_SP->SetLineStyle(2);
ALL_bx_Q1_SP->Draw();
ALL_bx_Q1->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q1,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q1_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q1.gif");
c1->Clear();

ALL_bx_Q2_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q2_SP->GetXaxis()->CenterLabels();
ALL_bx_Q2_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q2_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q2_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q2_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q2_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q2_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q2_SP->SetLineColor(kRed);
ALL_bx_Q2_SP->SetLineWidth(3);
ALL_bx_Q2_SP->SetLineStyle(2);
ALL_bx_Q2_SP->Draw();
ALL_bx_Q2->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q2,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q2_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q2.gif");
c1->Clear();

ALL_bx_Q3_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q3_SP->GetXaxis()->CenterLabels();
ALL_bx_Q3_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q3_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q3_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q3_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q3_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q3_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q3_SP->SetLineColor(kRed);
ALL_bx_Q3_SP->SetLineWidth(3);
ALL_bx_Q3_SP->SetLineStyle(2);
ALL_bx_Q3_SP->Draw();
ALL_bx_Q3->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q3,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q3_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q3.gif");
c1->Clear();

ALL_bx_Q4_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q4_SP->GetXaxis()->CenterLabels();
ALL_bx_Q4_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q4_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q4_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q4_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q4_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q4_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q4_SP->SetLineColor(kRed);
ALL_bx_Q4_SP->SetLineWidth(3);
ALL_bx_Q4_SP->SetLineStyle(2);
ALL_bx_Q4_SP->Draw();
ALL_bx_Q4->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q4,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q4_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q4.gif");
c1->Clear();

ALL_bx_Q5_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q5_SP->GetXaxis()->CenterLabels();
ALL_bx_Q5_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q5_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q5_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q5_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q5_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q5_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q5_SP->SetLineColor(kRed);
ALL_bx_Q5_SP->SetLineWidth(3);
ALL_bx_Q5_SP->SetLineStyle(2);
ALL_bx_Q5_SP->Draw();
ALL_bx_Q5->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q5,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q5_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q5.gif");
c1->Clear();

ALL_bx_Q6_SP->GetXaxis()->SetRangeUser(-5,5);
ALL_bx_Q6_SP->GetXaxis()->CenterLabels();
ALL_bx_Q6_SP->GetXaxis()->SetTitle("BX");
ALL_bx_Q6_SP->GetYaxis()->SetTitle("# entries");
ALL_bx_Q6_SP->GetXaxis()->CenterTitle(1); //centro el cartel del eje X
ALL_bx_Q6_SP->GetYaxis()->CenterTitle(1);
ALL_bx_Q6_SP->GetYaxis()->SetTitleOffset(1.4); //muevo el cartel del eje y a la izq
ALL_bx_Q6_SP->GetXaxis()->SetTitleOffset(1.2);
ALL_bx_Q6_SP->SetLineColor(kRed);
ALL_bx_Q6_SP->SetLineWidth(3);
ALL_bx_Q6_SP->SetLineStyle(2);
ALL_bx_Q6_SP->Draw();
ALL_bx_Q6->Draw("same");
leg = new TLegend(0.1,0.7,0.48,0.9);
leg->AddEntry(ALL_bx_Q6,"Blue= DT primitives");
leg->AddEntry(ALL_bx_Q6_SP,"Red = Superprimitives");
leg->Draw();
c1->SetLogy();
c1->Update();
c1->SaveAs("AllCh_BXdist_SPvsLegacy_Q6.gif");
c1->Clear();

}
