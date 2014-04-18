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

void Efficiency_tool_cris()
{

gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);

//creating the Canvas
TCanvas *c1 = new TCanvas("c1", "Data", 1000, 600);
//c1->cd(1);

//FIRST FILE************************+
TFile f("L1ITMuonBarrelTimingEffPlots_17aprilT0.root");
TEfficiency *eff=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_posPhi_W0_Sec1_St1");
eff->SetLineColor(kBlue);
eff->SetLineStyle(1);
eff->Draw();

//SECOND FILE************************+
TFile fsp("L1ITMuonBarrelTimingEffPlots_SP_17aprilT0.root");
TEfficiency *eff_SP=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_posPhi_W0_Sec1_St1");

eff_SP->SetLineColor(kRed);
eff_SP->SetLineStyle(2);
eff_SP->Draw("same");

leg1 = new TLegend(0.3,0.2,0.78,0.4);
leg1->AddEntry(eff,"Blue= DT primitives");
leg1->AddEntry(eff_SP,"Red = Superprimitives");
leg1->Draw();

c1->Update();

gPad->Update();
eff->GetPaintedGraph()->GetXaxis()->SetTitle("Track Local Phi Position");
eff->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at any BX");

TText *thlow=new TText(-10, 0.5, "W0 S1 MB1");
thlow->SetTextAlign(11);
thlow->Draw();

c1->SaveAs("Efficiency_phiPos_W0_Sec1_St1.gif");
c1->Clear();


//BX=0************************+
//FIRST FILE************************+
TEfficiency *effbx0=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_posPhi_BX0_W0_Sec1_St1");
effbx0->SetLineColor(kBlue);
effbx0->SetLineStyle(1);
effbx0->Draw();

//SECOND FILE************************+
TEfficiency *effbx0_SP=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_posPhi_BX0_W0_Sec1_St1");
effbx0_SP->SetLineColor(kRed);
effbx0_SP->SetLineStyle(2);

effbx0_SP->Draw("same");


leg2 = new TLegend(0.3,0.2,0.78,0.4);
leg2->AddEntry(effbx0,"Blue= DT primitives");
leg2->AddEntry(effbx0_SP,"Red = Superprimitives");
leg2->Draw();

c1->Update();

gPad->Update();
effbx0->GetPaintedGraph()->GetXaxis()->SetTitle("Track Local Phi Position");
effbx0->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at BX=0");

TText *thlow2=new TText(-10, 0.5, "W0 S1 MB1");
thlow2->SetTextAlign(11);
thlow2->Draw();

c1->SaveAs("Efficiency_atBX0_phiPos_W0_Sec1_St1.gif");
c1->Clear();


//PhiDIRECTION ************************+
//FIRST FILE************************+
TEfficiency *effdir=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_dirPhi_W0_Sec1_St1");
effdir->SetLineColor(kBlue);
effdir->SetLineStyle(1);
effdir->Draw();

//SECOND FILE************************+
TEfficiency *eff_SPdir=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_dirPhi_W0_Sec1_St1");
eff_SPdir->SetLineColor(kRed);
eff_SPdir->SetLineStyle(2);
eff_SPdir->Draw("same");


leg3 = new TLegend(0.3,0.2,0.78,0.4);
leg3->AddEntry(effdir,"Blue= DT primitives");
leg3->AddEntry(eff_SPdir,"Red = Superprimitives");
leg3->Draw();

c1->Update();

gPad->Update();
effdir->GetPaintedGraph()->GetXaxis()->SetTitle("Track Local Phi Direction");
effdir->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at any BX");

effdir->GetPaintedGraph()->SetTitle("W0 S1 MB1");


effdir->GetPaintedGraph()->GetXaxis()->SetRangeUser(-30,30);


TText *thlow3=new TText(-10, 0.5, "W0 S1 MB1");
thlow3->SetTextAlign(11);
thlow3->Draw();

c1->SaveAs("Efficiency_phiDir_W0_Sec1_St1.gif");
c1->Clear();


//PhiDIRECTION at BX=0 ************************+
//FIRST FILE************************+
TEfficiency *effbx0dir=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_dirPhi_BX0_W0_Sec1_St1");
effbx0dir->SetLineColor(kBlue);
effbx0dir->SetLineStyle(1);
effbx0dir->Draw();

//SECOND FILE************************+
TEfficiency *effbx0_SPdir=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_dirPhi_BX0_W0_Sec1_St1");
effbx0_SPdir->SetLineColor(kRed);
effbx0_SPdir->SetLineStyle(2);
effbx0_SPdir->Draw("same");

leg4 = new TLegend(0.3,0.2,0.78,0.4);
leg4->AddEntry(effbx0dir,"Blue= DT primitives");
leg4->AddEntry(effbx0_SPdir,"Red = Superprimitives");
leg4->Draw();

c1->Update();

gPad->Update();
effbx0dir->GetPaintedGraph()->GetXaxis()->SetTitle("Track Local Phi Direction");
effbx0dir->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at BX=0");

effbx0dir->GetPaintedGraph()->GetXaxis()->SetRangeUser(-30,30);

TText *thlow4=new TText(-10, 0.5, "W0 S1 MB1");
thlow4->SetTextAlign(11);
thlow4->Draw();

c1->SaveAs("Efficiency_atBX0_phiDir_W0_Sec1_St1.gif");
c1->Clear();



}
