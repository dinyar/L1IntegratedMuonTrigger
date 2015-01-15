#include <TROOT.h>
#include <TSystem.h>
#include<TString.h>
#include<TFile.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TF1.h>
#include<iostream>
#include<fstream>

#include<stdlib.h>

class WeightMaker {

private:
 
  TString outdir_;
  TFile * file_;
  std::string cut_;
  TString TailCut_;
public:
  WeightMaker(const TString & infile, const TString & outdir, std::string cut, TString TailCut);
 
  void run();

  void GetValues(std::string Type,	  
		 Float_t &WSigPhiBIn,
		 Float_t &WSigPhiBOut,
		 Float_t &WSigDeltaPhi,
		 Float_t &WFitPhiBIn,
		 Float_t &WFitPhiBOut,
		 Float_t &WFitDeltaPhi,
		 Float_t &WTailPhiBIn,
		 Float_t &WTailPhiBOut,
		 Float_t &WTailDeltaPhi,
		 Int_t Inch,Int_t OutCh,std::string DT1 ,std::string DT2,std::string wh,std::string sc);
  
  virtual ~WeightMaker() {};
};

WeightMaker::WeightMaker( const TString & infile, const TString & outdir, std::string cut, TString TailCut)
  :  outdir_( outdir ) , cut_(cut) , TailCut_(TailCut)
{
    file_ = TFile::Open( infile );
  TString dir = "mkdir -p " + outdir_;
  system( dir.Data() );
}



void WeightMaker:: GetValues(std::string Type,	 
			     Float_t &WSigPhiBIn,
			     Float_t &WSigPhiBOut,
			     Float_t &WSigDeltaPhi,
			     Float_t &WFitPhiBIn,
			     Float_t &WFitPhiBOut,
			     Float_t &WFitDeltaPhi,
			     Float_t &WTailPhiBIn,
			     Float_t &WTailPhiBOut,
			     Float_t &WTailDeltaPhi,
			     Int_t Inch,Int_t OutCh,std::string DT1 ,std::string DT2,std::string wh,std::string sc){
    
 std::string typePhi  ="DeltaPhiPtResol";
 std::string typeBIn  ="PhiBInPtResol";
 std::string typeBOut ="PhiBOutPtResol";
 

 std::string namefileBIn = "h"+typeBIn+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_;
 std::string namepathBIn = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileBIn;

 std::string namefileBOut = "h"+typeBOut+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_;
 std::string namepathBOut = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileBOut;

 std::string namefilePhi = "h"+typePhi+"Wh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_;
 std::string namepathPhi = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefilePhi; 

 std::string namefileDeltaBIn = "DeltaBInPtResolWh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_; 
 std::string namepathDeltaBIn = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileDeltaBIn;  

 std::string namefileDeltaBOut = "DeltaBOutPtResolWh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_; 
 std::string namepathDeltaBOut = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileDeltaBOut;  
 
 std::string namefileBInBOut = "BInBOutPtResolWh"+wh+"Sc"+sc+"inCh1outCh2"+DT1+DT2+"Cut_"+cut_; 
 std::string namepathBInBOut = "L1ITMuPtPlotter/Wh"+wh+"Sc"+sc+"/Resolutions/Cut_"+cut_+"/"+namefileBInBOut;  


  TH1F* HresPhiBIn = (TH1F*)file_->Get(namepathBIn.c_str());
  TH1F* HresPhiBOut = (TH1F*)file_->Get(namepathBOut.c_str());
  TH1F* HresDeltaPhi = (TH1F*)file_->Get(namepathPhi.c_str());

  TH2F* HresDeltaBIn = (TH2F*)file_->Get(namepathDeltaBIn.c_str()); 
  TH2F* HresDeltaBOut = (TH2F*)file_->Get(namepathDeltaBOut.c_str()); 
  TH2F* HresBInBOut = (TH2F*)file_->Get(namepathBInBOut.c_str()); 


  TF1 *f1 = new TF1("f1","gaus",-3,1);
  
  TString SaveName = outdir_+"/SigmaWe"+wh+"Sc"+sc+"WeightPlots.root";

  TFile *FileOut = new TFile(SaveName.Data(),"UPDATE");
 
  HresPhiBIn->Fit("f1");
  Float_t FitSigBIn = f1->GetParameter(2);
  HresPhiBIn->Write();
  
  HresPhiBOut->Fit("f1");   
  Float_t FitSigBOut = f1->GetParameter(2);
  HresPhiBOut->Write();   
  
  HresDeltaPhi->Fit("f1");
  Float_t FitSigDPhi = f1->GetParameter(2);
  HresDeltaPhi->Write();
 

  Float_t VarDPhi = 1./HresDeltaBIn->GetCovariance(1,1);
  Float_t VarBOut = 1./HresDeltaBOut->GetCovariance(2,2);
  Float_t VarBIn = 1./HresBInBOut->GetCovariance(1,1);

  Float_t TailCut = TailCut_.Atof();

  TAxis *xaxisDPhi = HresDeltaPhi->GetXaxis();
  Int_t BinCutDPhi = xaxisDPhi->FindBin(TailCut);
  Float_t TailDPhi = 1./HresDeltaPhi->Integral(0,BinCutDPhi);
  TailDPhi=TailDPhi*TailDPhi;
  std::cout<<"Bin "<<BinCutDPhi<<" integral "<<HresDeltaPhi->Integral(0,BinCutDPhi)<<" weight "<<TailDPhi<<std::endl;


  TAxis *xaxisBIn = HresPhiBIn->GetXaxis();
  Int_t BinCutPhiBIn = xaxisBIn->FindBin(TailCut);
  Float_t TailBIn = 1./HresPhiBIn->Integral(0,BinCutPhiBIn);
  TailBIn=TailBIn*TailBIn;


  TAxis *xaxisBOut = HresPhiBOut->GetXaxis();
  Int_t BinCutPhiBOut = xaxisBOut->FindBin(TailCut);
  Float_t TailBOut = 1./HresPhiBOut->Integral(0,BinCutPhiBOut);
  TailBOut=TailBOut*TailBOut;

  // GP  weights for weighted average for correlated variables, not used yet

  //   Float_t VarDeltaBIn = 1./HresDeltaBIn->GetCovariance();
  //   Float_t VarDeltaBOut = 1./HresDeltaBOut->GetCovariance();
  //   Float_t VarBInBOut = 1./HresBInBOut->GetCovariance();
 
  
  //   RMSPhiBIn = (VarDeltaBIn+VarBInBOut+RMSBIn)/(VarDeltaBOut+VarDeltaBIn+VarBInBOut+RMSBIn+RMSBOut+RMSDelta);
  //   RMSPhiBOut = (VarDeltaBOut+VarBInBOut+RMSBOut)/(VarDeltaBOut+VarDeltaBIn+VarBInBOut+RMSBIn+RMSBOut+RMSDelta);
  //   RMSDeltaPhi = (VarDeltaBIn+VarDeltaBOut+RMSDelta)/(VarDeltaBOut+VarDeltaBIn+VarBInBOut+RMSBIn+RMSBOut+RMSDelta);
  
  //   RMSPhiBIn = HresPhiBIn->GetRMS();
  //   RMSPhiBOut = HresPhiBOut->GetRMS();
  //   RMSDeltaPhi = HresDeltaPhi->GetRMS();
  //   RMSGMT = HresGMT->GetRMS();  


  if(Type=="WMean"){
    
    Float_t SumWVar = VarDPhi+VarBIn+VarBOut;
    Float_t SumWFit = FitSigDPhi+FitSigBIn+FitSigBOut;
    Float_t SumWTail = TailDPhi+TailBIn+TailBOut;
   
    
    WSigPhiBIn=VarBIn/SumWVar;
    WSigPhiBOut=VarBOut/SumWVar;
    WSigDeltaPhi=VarDPhi/SumWVar;

    WFitPhiBIn=FitSigBIn/SumWFit;
    WFitPhiBOut=FitSigBOut/SumWFit;
    WFitDeltaPhi=FitSigDPhi/SumWFit;

    WTailPhiBIn=TailBIn/SumWTail;
    WTailPhiBOut=TailBOut/SumWTail;
    WTailDeltaPhi=TailDPhi/SumWTail;
}
  

  else if(Type=="WMeanIn"){
    
    Float_t SumWVarIn = VarDPhi+VarBIn;
    Float_t SumWFitIn = FitSigDPhi+FitSigBIn;
    Float_t SumWTailIn = TailDPhi+TailBIn;
    
    WSigPhiBIn=VarBIn/SumWVarIn;
    WSigPhiBOut=0;
    WSigDeltaPhi=VarDPhi/SumWVarIn;

    WFitPhiBIn=FitSigBIn/SumWFitIn;
    WFitPhiBOut=0;
    WFitDeltaPhi=FitSigDPhi/SumWFitIn;

    WTailPhiBIn=TailBIn/SumWTailIn;
    WTailPhiBOut=0;
    WTailDeltaPhi=TailDPhi/SumWTailIn;

  }
  

  else if(Type=="WMeanOut"){
    
    Float_t SumWVarOut= VarDPhi+VarBOut;
    Float_t SumWFitOut = FitSigDPhi+FitSigBOut;
    Float_t SumWTailOut = TailDPhi+TailBOut;
    
    WSigPhiBIn=0;
    WSigPhiBOut=VarBOut/SumWVarOut;
    WSigDeltaPhi=VarDPhi/SumWVarOut;

    WFitPhiBIn=0;
    WFitPhiBOut=FitSigBOut/SumWFitOut;
    WFitDeltaPhi=FitSigDPhi/SumWFitOut;

    WTailPhiBIn=0;
    WTailPhiBOut=TailBOut/SumWTailOut;
    WTailDeltaPhi=TailDPhi/SumWTailOut;

  }

  HresPhiBIn->Delete();
  HresPhiBOut->Delete();
  HresDeltaPhi->Delete();
 
  
  FileOut->Close();
}

void WeightMaker::run( )
{
  

  std::string wheels[6] = { "-3", "-2", "-1", "1","2","3" };
  std::string ref[3] = { "DTCORR", "DTIN", "DTOUT" };
  std::string Type[3] = { "WMean", "WMeanIn", "WMeanOut" };
  
  Float_t PhiBInRMS;
  Float_t PhiBOutRMS; 
  Float_t DeltaPhiRMS;
  
  Float_t PhiBInSigma;
  Float_t PhiBOutSigma; 
  Float_t DeltaPhiSigma;

  Float_t PhiBInTail;
  Float_t PhiBOutTail; 
  Float_t DeltaPhiTail;

  // int sc = 0;
  
  int inCh[2] = { 1, 2 };
  int outCh[2] = { 1, 2 };
  
  
  
   for ( int sector = 0; sector < 12; ++sector ) { 
    for ( size_t w = 0; w < 6; ++w ) {
      
      for ( size_t ch1 = 0; ch1 < 1; ++ch1 ) {
	for ( size_t ch2 = 1; ch2 < 2; ++ch2 ) {
	  
	  TString NamefileRMS = Form("/RMSWh%sSc%i.txt",wheels[w].c_str(),sector);
	  NamefileRMS=outdir_+NamefileRMS;

	  TString NamefileSigma = Form("/SigmaWh%sSc%i.txt",wheels[w].c_str(),sector);
	  NamefileSigma=outdir_+NamefileSigma;

	  TString NamefileTail = Form("/TailWh%sSc%i.txt",wheels[w].c_str(),sector);
	  NamefileTail=outdir_+NamefileTail;
	  
	  
	  
	  ofstream outputFileRMS(NamefileRMS.Data());
	  outputFileRMS<<"# Type DT1 DT2  Weight DeltaPhi   Weight PhiBIn   Weight PhiBOut "<<std::endl;
	  
	  ofstream outputFileSigma(NamefileSigma.Data());
	  outputFileSigma<<"# Type DT1 DT2  Weight DeltPhi   Weight PhiBIn   Weight PhiBOut "<<std::endl;

	  ofstream outputFileTail(NamefileTail.Data());
	  outputFileTail<<"# Type DT1 DT2  Weight DeltPhi   Weight PhiBIn   Weight PhiBOut "<<std::endl;
	  
	  
	  for ( size_t p = 0; p <=2; ++p ) {
	    for ( size_t q = 0; q <=2; ++q ) {
	      for( size_t t = 0; t<= 2; t++) {
		
		GetValues(Type[t],PhiBInRMS,PhiBOutRMS,DeltaPhiRMS,PhiBInSigma,PhiBOutSigma,DeltaPhiSigma,PhiBInTail,PhiBOutTail,DeltaPhiTail,inCh[ch1],outCh[ch2],ref[p],ref[q],wheels[w], "0"); //GP so far use sector 0 for every sectors, 
		
		std::cout<<ref[p]<<" "<<ref[q]<<" "<<" "<<DeltaPhiRMS<<" "<<PhiBInRMS<<" "<<PhiBOutRMS<<std::endl;
		
		outputFileRMS<<Type[t]<<" "<<ref[p]<<" "<<ref[q]<<" "<<" "<<DeltaPhiRMS<<" "<<PhiBInRMS<<" "<<PhiBOutRMS<<std::endl;
		outputFileSigma<<Type[t]<<" "<<ref[p]<<" "<<ref[q]<<" "<<" "<<DeltaPhiSigma<<" "<<PhiBInSigma<<" "<<PhiBOutSigma<<std::endl;
		outputFileTail<<Type[t]<<" "<<ref[p]<<" "<<ref[q]<<" "<<" "<<DeltaPhiTail<<" "<<PhiBInTail<<" "<<PhiBOutTail<<std::endl;
	      }
	    }
	  }
	}
      }
    } 
  }
}


int main(int argc, char** argv) {

  if ( argc < 5 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <input> <outdir>  <usualy cut[=30,35,40,45]> <tail cut>"
              << std::endl << std::endl;
   
    return -1;
  }
 
  std::string arg3 = argv[3];
  std::cout<<" arg 3 "<<arg3.c_str()<<std::endl;


  WeightMaker Weight(argv[1],argv[2],argv[3],argv[4]);
  Weight.run();
  
  return 0;
}
