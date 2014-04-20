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

void Efficiency_tool_cris_T0()
{

gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);

//creating the Canvas
TCanvas *c1 = new TCanvas("c1", "Data", 1000, 600);
//c1->cd(1);

//FIRST FILE************************+
TFile f("L1ITMuonBarrelTimingEffPlots_17aprilT0.root");
TEfficiency *effT0=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_T0_W0_Sec1_St1");
effT0->SetLineColor(kBlue);
effT0->SetLineStyle(1);
effT0->Draw();

//SECOND FILE************************+
TFile fsp("L1ITMuonBarrelTimingEffPlots_SP_17aprilT0.root");
TEfficiency *eff_SPT0=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_T0_W0_Sec1_St1");

eff_SPT0->SetLineColor(kRed);
eff_SPT0->SetLineStyle(2);
eff_SPT0->Draw("same");

leg1 = new TLegend(0.1,0.75,0.3,0.9);
leg1->AddEntry(effT0,"Blue= DT primitives");
leg1->AddEntry(eff_SPT0,"Red = Superprimitives");
leg1->Draw();

c1->Update();

gPad->Update();
effT0->GetPaintedGraph()->GetXaxis()->SetTitle("Segment T0");
effT0->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at any BX");

TText *thlow=new TText(-10, 0.1, "W0 S1 MB1");
thlow->SetTextAlign(11);
thlow->Draw();

c1->SaveAs("Efficiency_T0_W0_Sec1_St1.gif");
c1->Clear();


//BX=0************************+
//FIRST FILE************************+
TEfficiency *effbx0T0=(TEfficiency*)f.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_T0_BX0_W0_Sec1_St1");
effbx0T0->SetLineColor(kBlue);
effbx0T0->SetLineStyle(1);
effbx0T0->Draw();

//SECOND FILE************************+
TEfficiency *effbx0_SPT0=(TEfficiency*)fsp.Get("L1ITMuTimingPlotter/Wheel0/Sector1/Station1/TP_Efficiency_T0_BX0_W0_Sec1_St1");
effbx0_SPT0->SetLineColor(kRed);
effbx0_SPT0->SetLineStyle(2);

effbx0_SPT0->Draw("same");


leg2 = new TLegend(0.1,0.75,0.3,0.9);
leg2->AddEntry(effbx0T0,"Blue= DT primitives");
leg2->AddEntry(effbx0_SPT0,"Red = Superprimitives");
leg2->Draw();

c1->Update();

gPad->Update();
effbx0T0->GetPaintedGraph()->GetXaxis()->SetTitle("Segment T0");
effbx0T0->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency at BX=0");

TText *thlow2=new TText(-10, 0.1, "W0 S1 MB1");
thlow2->SetTextAlign(11);
thlow2->Draw();

c1->SaveAs("Efficiency_atBX0_T0_W0_Sec1_St1.gif");
c1->Clear();

}
