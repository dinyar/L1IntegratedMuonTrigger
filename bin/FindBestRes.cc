
// 
// Class: FindBestRes
//
// Info: Find the weighted average the provides best resolution and print results in a file.
//
// Author: Gian Luca Pinna Angioni
//


#include <TROOT.h>
#include <TSystem.h>
#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TF1.h>
#include<iostream>
#include<fstream>
#include<TCanvas.h>
#include<stdlib.h>
#include <map>



class FindBestRes {

private:
 
  TString outdir_;
  TFile * file_;
  std::string cut_;
public:
  FindBestRes(const TString & infile, const TString & outdir, std::string cut);
 
  void run();

  Float_t Round(Float_t f);

  bool MakePlotsRes(std::string type,  
		    Float_t & GainRMS,
		    Float_t & GainSigma,
		    Int_t Inch,Int_t OutCh,std::string DT1 ,std::string DT2,std::string wh,std::string sc);

  virtual ~FindBestRes() {};
};

FindBestRes::FindBestRes( const TString & infile, const TString & outdir, std::string cut)
:  outdir_( outdir ) , cut_(cut)
{
    file_ = TFile::Open( infile );
  TString dir = "mkdir -p " + outdir_;
  system( dir.Data() );

  TString dirPlot = "mkdir -p plots/CombResults";
  system( dirPlot.Data() );
}


Float_t FindBestRes::Round(Float_t f){
  Float_t y, d;
  y = f*100;
  d = y - (int)y;
  y = (float)(int)(f*100)/100;
  if (d > 0.5)
    y += 0.01;
  return y;
}


bool FindBestRes:: MakePlotsRes(std::string type, 
				Float_t & GainRMS,
				Float_t & GainSigma,
				Int_t Inch,Int_t OutCh,std::string DT1 ,std::string DT2,std::string wh,std::string sc){
  
  
  std::string namefileGMT = "hGMTPtResolWh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_;
  std::string namepathGMT = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileGMT;  
  
  std::string namefilehRes = "h"+type+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_;
  std::string namepathhRes = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefilehRes;
  
  std::cout<<"opening histos :"<<std::endl<<namepathGMT.c_str()<<std::endl<<namepathhRes.c_str()<<std::endl<<std::endl;
  
  TH1F* HresGMT = (TH1F*)file_->Get(namepathGMT.c_str());
  TH1F* HresRes = (TH1F*)file_->Get(namepathhRes.c_str());
  
  Int_t entries = HresGMT->GetEntries();
  if(entries>0){
    
    TF1 *f1 = new TF1("f1","gaus",-3,1);
    
    HresGMT->Fit("f1");  
    Float_t SigmaGMT = f1->GetParameter(2);
    HresGMT->GetFunction("f1")->SetLineColor(6);
    HresRes->Fit("f1");
   
    Float_t SigmahRes = f1->GetParameter(2);
    Float_t RMSGMT=     HresGMT->GetRMS(); 
    Float_t RMShRes=    HresRes->GetRMS(); 
    
    GainSigma = Round((SigmaGMT/SigmahRes - 1)*100);
    GainRMS = Round((RMSGMT/RMShRes - 1)*100);
   

    HresGMT->Delete();
    HresRes->Delete();

    return 1;  
  }
  else {
    GainSigma = 1;
    GainRMS =1;
    std::cout<<" No Entries in "<<namepathhRes.c_str()<<std::endl<<" gain  = 1"<<std::endl;
    return 0;
  } 
}




 void FindBestRes::run( )
{

  std::string wheels[6] = { "-3", "-2", "-1", "1","2","3" };
  std::string ref[3] = { "DTCORR", "DTIN", "DTOUT" };
  
  
  std::map<int,std::string> Type;
  Type[1]= "DeltaPhiPtResol";
  Type[2]= "PhiBInPtResol";
  Type[3]= "PhiBOutPtResol";
  Type[4]= "MeanPtResol";
  Type[5]= "RMSWMeanPtResol";
  Type[6]= "RMSWMeanInPtResol";
  Type[7]= "RMSWMeanOutPtResol";
  Type[8]= "SigmaWMeanPtResol";
  Type[9]= "SigmaWMeanInPtResol";
  Type[10]= "SigmaWMeanOutPtResol";

  Float_t GainRMS;
  Float_t GainSigma;

  
  int inCh[2] = { 1, 2 };
  int outCh[2] = { 1, 2 };

  
  
  for ( int sector = 0; sector < 12; ++sector ) {
   
    for ( size_t w = 0; w < 6; ++w ) {
      for ( size_t ch1 = 0; ch1 < 1; ++ch1 ) {
	for ( size_t ch2 = 1; ch2 < 2; ++ch2 ) {
	  
	 
	  TH1F *hResultSigma = new TH1F("hSigma","Results #sigma gain deff",9,1,10);
	  TH1F *hResultRMS = new TH1F("hRMS","Results RMS gain deff",9,1,10);	  	 
 
	  hResultRMS->SetStats(0);
	  hResultSigma->SetStats(0);

	  TAxis *axXRMS = hResultRMS->GetXaxis();
	  TAxis *axXSigma = hResultSigma->GetXaxis();
	  	 
	  for(int i = 1; i<=3;i++) {
	    for(int j = 1; j<=3;j++) {
	      axXRMS->SetBinLabel(3*i+j-3,(ref[i-1]+ref[j-1]).c_str());
	      axXSigma->SetBinLabel(3*i+j-3,(ref[i-1]+ref[j-1]).c_str());
	    }
	  }
 
	  TString Namefile = Form("/Wh%sSc%i.txt",wheels[w].c_str(),sector);
	  std::cout<<" name file "<<Namefile.Data()<<std::endl;
	  Namefile=outdir_+Namefile;
	  
	  ofstream outputFile(Namefile.Data());
	  
	  outputFile<<" Cat1 Cat2  RMSGain method  SigmaGain method"<<std::endl;  	           
	  
 	  for ( size_t p = 0; p <=2; ++p ) {
	    for ( size_t q = 0; q <=2; ++q ) {
	      
	      Float_t gainGoodRMS = -1000;
	      Float_t gainGoodSigma = -1000;

	      
	      std::string GoodTypeRMS;
	      std::string GoodTypeSigma;
	      
	      bool Result = 0;
	      for(int k = 1;k<=10;k++){

		if(q!=0&& (k==3||k==4||k==5||k==7||k==8||k==10)) continue;
		if(p!=0&& (k==2||k==4||k==5||k==6||k==8||k==9)) continue;
		Result = MakePlotsRes(Type[k],GainRMS,GainSigma,inCh[ch1],outCh[ch2],ref[p],ref[q],wheels[w], "0");
		
	       	       if(Result){ 
				   std::cout<<"Qual "<< (ref[p]+ref[q]).c_str()<<" type "<<Type[k].c_str()<<" gain RMS "<<GainRMS<<std::endl;
				   if(GainRMS>gainGoodRMS) {
				     gainGoodRMS=GainRMS; 
				     GoodTypeRMS=Type[k];
				   }
				   
				   if(GainSigma>gainGoodSigma) {
				     gainGoodSigma=GainSigma; 
				     GoodTypeSigma=Type[k];
				   }  
		       }	     
		       else{
			 gainGoodSigma = 0;
			 gainGoodRMS = 0;
			 GoodTypeRMS = Type[1];
			 GoodTypeSigma = Type[1];
			 std::cout<<" no possible comparison, best type is set automatically Delta Phi "<<std::endl;
		       }
	      }
	     
	      hResultRMS->Fill(3*p+q+1,gainGoodRMS);
	      hResultSigma->Fill(3*p+q+1,gainGoodSigma);
	     
	      hResultRMS->SetYTitle("gain"); 
	      hResultSigma->SetYTitle("gain"); 
	      
	      TCanvas *cRes =new TCanvas("Res","results",1000,600);
	      
	      cRes->Divide(2,1);
	      cRes->cd(1);
	      hResultRMS->Draw();
	      
	      cRes->cd(2);
	      hResultSigma->Draw(); 
	     
	      cRes->SaveAs(("plots/CombResults/ResultWh"+wheels[w]+".png").c_str());
	     
	      outputFile<<(ref[p]).c_str()<<" "<<(ref[q]).c_str()<<" "<<gainGoodRMS<<"  "<<GoodTypeRMS<<" "<<gainGoodSigma<<" "<<GoodTypeSigma<<std::endl;      
	      
	    }
	  }
	}
      }
    } 
  }
}


int main(int argc, char** argv) {
  
  if ( argc < 4 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <input> <outdir>  <usually cut[=30,35,40,45]>"
              << std::endl << std::endl;
   
    return -1;
  }
 
  std::string arg3 = argv[3];
  std::cout<<" arg 3 "<<arg3.c_str()<<std::endl;



  FindBestRes Find(argv[1],argv[2],argv[3]);
 
  Find.run();
  
  return 0;
}
