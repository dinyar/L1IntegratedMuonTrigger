
void makePrimitivesPlots()
{
  TFile * f2 = TFile::Open( "L1ITMuonBarrelPrimitivePlots.root");

  TDirectory * d005 = (TDirectory *) f2->Get( "L1ITMuPlotter");
  TString name;



  ////////////////////
  /// uncorrelated
  name = "uncorrelated";
  TCanvas * uncorrelated = new TCanvas(name, name);
  uncorrelated->cd();
  uncorrelated->SetGridy();
  TH1* unc = d005->Get(name);
  unc->SetLineWidth(2);
  // unc->SetStats( 0 );
  unc->Draw();
  unc->Draw("sametext");
  uncorrelated->SaveAs( name + ".png" );


  ////////////////////
  /// uncorrelatedNorm
  name = "uncorrelatedNorm";
  TCanvas * uncorrelatedNorm = new TCanvas(name, name);
  uncorrelatedNorm->cd();
  uncorrelatedNorm->SetGridy();
  TH1* uncNorm = d005->Get(name);
  uncNorm->GetYaxis()->SetRangeUser( 0., 1.1 );
  uncNorm->GetXaxis()->SetNdivisions( 11 );
  uncNorm->SetLineWidth(2);
  // uncNorm->SetStats( 0 );
  uncNorm->Draw();
  uncNorm->Draw("sametext");
  uncorrelatedNorm->SaveAs( name + ".png" );


  ////////////////////
  /// bendingAngleDiff
  name = "bendingAngleDiff";
  TCanvas * bendingAngleDiff = new TCanvas(name, name);
  bendingAngleDiff->cd();
  TH2* bendAngleDiff = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  bendAngleDiff->Draw();
  bendingAngleDiff->SaveAs( name + ".png" );


  ////////////////////
  /// radialAngleDiff
  name = "radialAngleDiff";
  TCanvas * radialAngleDiff = new TCanvas(name, name);
  radialAngleDiff->cd();
  TH1* radAngleDiff = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  radAngleDiff->Draw();
  radialAngleDiff->SaveAs( name + ".png" );


  ////////////////////
  /// bendingAngle
  name = "bendingAngle";
  TCanvas * bendingAngle = new TCanvas(name, name);
  bendingAngle->cd();
  TH2* bendAngle = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  bendAngle->Draw("box");
  bendingAngle->SaveAs( name + ".png" );


  ////////////////////
  /// radialAngle
  name = "radialAngle";
  TCanvas * radialAngle = new TCanvas(name, name);
  radialAngle->cd();
  TH1* radAngle = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  radAngle->Draw("box");
  radialAngle->SaveAs( name + ".png" );


  ////////////////////
  /// bendingAnglePhysDiff
  name = "bendingAnglePhysDiff";
  TCanvas * bendingAnglePhysDiff = new TCanvas(name, name);
  bendingAnglePhysDiff->cd();
  TH2* bendAnglePhysDiff = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  bendAnglePhysDiff->Draw();
  bendingAnglePhysDiff->SaveAs( name + ".png" );


  ////////////////////
  /// radialAnglePhysDiff
  name = "radialAnglePhysDiff";
  TCanvas * radialAnglePhysDiff = new TCanvas(name, name);
  radialAnglePhysDiff->cd();
  TH1* radAnglePhysDiff = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  radAnglePhysDiff->Draw();
  radialAnglePhysDiff->SaveAs( name + ".png" );


  ////////////////////
  /// bendingAnglePhys
  name = "bendingAnglePhys";
  TCanvas * bendingAnglePhys = new TCanvas(name, name);
  bendingAnglePhys->cd();
  TH2* bendAnglePhys = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  bendAnglePhys->Draw("box");
  bendingAnglePhys->SaveAs( name + ".png" );


  ////////////////////
  /// radialAnglePhys
  name = "radialAnglePhys";
  TCanvas * radialAnglePhys = new TCanvas(name, name);
  radialAnglePhys->cd();
  TH1* radAnglePhys = d005->Get(name);
  gStyle->SetOptStat( 111111 );
  radAnglePhys->Draw("box");
  radialAnglePhys->SaveAs( name + ".png" );

  return;

  ////////////////////
  /// dtQuality
  name = "dtQuality";
  TCanvas * dtQuality = new TCanvas(name, name);
  dtQuality->SetGridx();
  dtQuality->SetGridy();
  dtQuality->cd();
  // dtQuality->SetLogz();
  TH1* qual = d005->Get(name);
  qual->GetZaxis()->SetRangeUser( 0., 0.5 );
  qual->SetStats( 0 );
  qual->Draw("colztext");
  dtQuality->SaveAs( name + ".png" );

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
  dtQualityNew->SaveAs( name + ".png" );


}
