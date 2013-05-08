
TDirectory * d2;
TDirectory * d02;
TDirectory * d005;
TDirectory * d002;

void drawConfirmed(const TString & name )
{
  TCanvas * confirmed =  new TCanvas(name, name);
  confirmed->cd();
  // confirmed->SetLogy();
  confirmed->SetGridy();
  TH1* st1_2 = d2->Get(name);
  TH1* st1_02 = d02->Get(name);
  TH1* st1_002 = d002->Get(name);
  TH1* st1_005 = d005->Get(name);

  st1_2->SetLineWidth(2);
  st1_02->SetLineWidth(2);
  st1_002->SetLineWidth(2);
  st1_005->SetLineWidth(2);

  st1_2->SetStats(0);
  st1_02->SetStats(0);
  st1_002->SetStats(0);
  st1_005->SetStats(0);

  st1_02->SetLineColor( kGreen );
  st1_002->SetLineColor( kOrange );
  st1_005->SetLineColor( kRed );

  st1_2->Draw();
  // st1_2->Draw("sametext");
  st1_02->Draw("same");
  st1_002->Draw("same");
  st1_005->Draw("same");

  TLegend * leg = new TLegend(.7, 0.75, 0.9, 0.9  );
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
//   leg->SetTextSize( 0.035 );
//   leg->SetBorderSize(0);
  leg->AddEntry(st1_2, "|#Delta#phi| < 2", "pl");
  leg->AddEntry(st1_02, "|#Delta#phi| < 0.2", "pl");
  leg->AddEntry(st1_005, "|#Delta#phi| < 0.05", "pl");
  leg->AddEntry(st1_002, "|#Delta#phi| < 0.02", "pl");
  leg->Draw();

  confirmed->SaveAs( name + ".png" );
}



void makePlots()
{

  TFile * f2 = TFile::Open( "L1ITMuonBarrelPlots2.root");
  TFile * f02 = TFile::Open( "L1ITMuonBarrelPlots0.2.root");
  TFile * f005 = TFile::Open( "L1ITMuonBarrelPlots0.05.root");
  TFile * f002 = TFile::Open( "L1ITMuonBarrelPlots0.02.root");

  d2 = (TDirectory *) f2->Get( "L1ITMuPlotter");
  d02 = (TDirectory *) f02->Get( "L1ITMuPlotter");
  d005 = (TDirectory *) f005->Get( "L1ITMuPlotter");
  d002 = (TDirectory *) f002->Get( "L1ITMuPlotter");
  TString name;


  //////////////////
  /// Confirmed
  drawConfirmed( "confirmed_st1" );
  drawConfirmed( "confirmed_st2" );
  drawConfirmed( "confirmed_st3" );
  drawConfirmed( "confirmed_st4" );

  TH1* conf = d005->Get( "confirmed_st1" )->Clone();
  TH1* conf_2 = d005->Get( "confirmed_st2" );
  TH1* conf_3 = d005->Get( "confirmed_st3" );
  TH1* conf_4 = d005->Get( "confirmed_st4" );
  conf->Add( conf_2 );
  conf->Add( conf_3 );
  conf->Add( conf_4 );

  name = "confirmed_all";
  TCanvas * confirmed =  new TCanvas(name, name);
  confirmed->cd();
  // confirmed->SetLogy();
  confirmed->SetGridy();
  conf->SetLineWidth(2);
  conf->Draw();
  conf->Draw("sametext");
  conf->SetTitle("confirmed in all stations");
  conf->SetName("confirmed in all stations");
  confirmed->SaveAs( name + "_dphi0.05.png" );

  ////////////////////
  /// deltaPhiBin
  name = "deltaPhiBin";
  TCanvas * deltaPhiBin = new TCanvas(name, name);
  deltaPhiBin->cd();
  TH1* dphibin = d2->Get(name);
  dphibin->SetLineWidth(2);
  dphibin->Draw("text");
  dphibin->Draw("same");
  deltaPhiBin->SaveAs( name + ".png" );

  ////////////////////
  /// deltaPhi
  name = "deltaPhi";
  TCanvas * deltaPhi = new TCanvas(name, name);
  deltaPhi->cd();
  TH1* dphi = d2->Get(name);
  dphi->SetLineWidth(2);
  dphi->Draw();
  deltaPhi->SaveAs( name + ".png" );

  ////////////////////
  /// deltaPhiDt
  name = "deltaPhiDt";
  TCanvas * deltaPhiDt = new TCanvas(name, name);
  deltaPhiDt->cd();
  TH1* dphiDt = d005->Get(name);
  dphiDt->SetLineWidth(2);
  dphiDt->Draw();
  deltaPhiDt->SaveAs( name + ".png" );

  /////////////////////////
  /// rpcInHitsPerDtseg
  name = "rpcInHitsPerDtseg";
  TCanvas * rpcInHitsPerDtseg = new TCanvas(name, name);
  rpcInHitsPerDtseg->SetGridy();
  rpcInHitsPerDtseg->SetLogy();
  rpcInHitsPerDtseg->cd();
  TH1* rpcin = d005->Get(name);
  rpcin->SetLineWidth(2);
  rpcin->Draw();
  rpcInHitsPerDtseg->SaveAs( name + "_dphi0.05.png" );

  /////////////////////////
  /// rpcOutHitsPerDtseg
  name = "rpcOutHitsPerDtseg";
  TCanvas * rpcOutHitsPerDtseg =  new TCanvas(name, name);
  rpcOutHitsPerDtseg->SetGridy();
  rpcOutHitsPerDtseg->SetLogy();
  rpcOutHitsPerDtseg->cd();
  TH1* rpcout = d005->Get(name);
  rpcout->SetLineWidth(2);
  rpcout->Draw();
  rpcOutHitsPerDtseg->SaveAs( name + "_dphi0.05.png" );

  ////////////////////
  /// dtQuality_st1
  name = "dtQuality_st1";
  TCanvas * dtQuality_st1 = new TCanvas(name, name);
  dtQuality_st1->SetGridx();
  dtQuality_st1->SetGridy();
  dtQuality_st1->cd();
  // dtQuality_st1->SetLogz();
  TH1* qual_st1 = d005->Get(name);
  qual_st1->GetZaxis()->SetRangeUser( 0., 1. );
  qual_st1->SetStats( 0 );
  qual_st1->Draw("colztext");
  dtQuality_st1->SaveAs( name + "_dphi0.05.png" );


  ////////////////////
  /// dtQuality_st2
  name = "dtQuality_st2";
  TCanvas * dtQuality_st2 = new TCanvas(name, name);
  dtQuality_st2->SetGridx();
  dtQuality_st2->SetGridy();
  dtQuality_st2->cd();
  // dtQuality_st2->SetLogz();
  TH1* qual_st2 = d005->Get(name);
  qual_st2->GetZaxis()->SetRangeUser( 0., 1. );
  qual_st2->SetStats( 0 );
  qual_st2->Draw("colztext");
  dtQuality_st2->SaveAs( name + "_dphi0.05.png" );


  ////////////////////
  /// dtQualityNew
  name = "dtQualityNew";
  TCanvas * dtQualityNew = new TCanvas(name, name);
  dtQualityNew->SetGridx();
  dtQualityNew->SetGridy();
  dtQualityNew->cd();
  // dtQualityNew->SetLogz();
  TH1* qualNew = d005->Get(name);
  qualNew->GetZaxis()->SetRangeUser( 0., 0.5 );
  qualNew->SetStats( 0 );
  qualNew->Draw("colztext");
  dtQualityNew->SaveAs( name + "_dphi0.05.png" );

}
