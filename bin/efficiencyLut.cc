#include <TROOT.h>
#include <TSystem.h>

#include<TString.h>
#include<TFile.h>
#include<TH2F.h>
#include<iostream>
#include<fstream>

#include<stdlib.h>

class efficiency {

private:
  std::ofstream testfile_;
  double target_;
  TString outdir_;
  TFile * file_;
  
public:
  efficiency( const TString & infile, const TString & outfile, const TString & outdir, double target );
  void run(bool);
  void write( TH2F * plot, std::ofstream & outfile, int inCh, int outCh,
	      const TString & ref1, const TString & ref2 );
  void write2( TH2F * plot, std::ofstream & outfile, int inCh, int outCh,
	       const TString & ref1, const TString & ref2 );
  virtual ~efficiency() {};
};


efficiency::efficiency(const TString & infile, const TString & outfile, const TString & outdir, double target )
  : testfile_( outfile.Data() ), target_( target ), outdir_( outdir )
{
  file_ = TFile::Open( infile );
  TString dir = "mkdir -p ../Luts/" + outdir_;
  system( dir.Data() );
}


void efficiency::write2( TH2F * plot, std::ofstream & outfile, int inCh, int outCh,
			const TString & ref1, const TString & ref2 )
{

  TAxis * xAxis = plot->GetXaxis();
  TAxis * yAxis = plot->GetYaxis();

  size_t xBins = plot->GetNbinsX();
  size_t yBins = plot->GetNbinsY();
//   double binWidth = plot->GetYaxis()->GetBinWidth(1);
//   double binStart = plot->GetXaxis()->GetXmin();
  // std::cout << binWidth << '\t' << binStart << '\t' << xBins << '\t' <<yBins << std::endl;
  double * xSum = new double[xBins];

  for ( size_t x = 1; x <= xBins; ++x ) {

    xSum[x] = 0;
    for ( size_t z = x; z <= xBins; ++z ) {
      for ( size_t y = 1; y <= yBins; ++y ) {
	xSum[x] += plot->GetBinContent( x, y );
      }
    }

    double eff = 0;
    double effPrev = 0;
    size_t y = 1;
    for ( ; y <= yBins; ++y ) {
      /// fill the y bin content, starting from threshold
      for ( size_t z = x; z <= xBins; ++z ) {
	double binContent = plot->GetBinContent( x, y );
	if ( ! binContent ) continue;
	eff += binContent / xSum[x];
      }

      if ( eff >= target_ ) {
	if ( target_ - effPrev < eff - target_ ) {
	  eff = effPrev;
	  --y;
	}
	break;
      } 
      effPrev = eff;
    }
    // std::cout << binWidth << '\t' << binStart << '\t' << xBins << '\t' << binStart + binWidth*x  << '\t' << xSum << std::endl;
    testfile_ << plot->GetName() << '\t' << plot->GetBinCenter( x ) << '\t' << xSum[x] << '\t' << y << '\t' << effPrev << '\t' << eff << std::endl; 
    outfile << inCh << '\t' << outCh << '\t' << ref1 << '\t' << ref2 << '\t'
	    << xAxis->GetBinCenter( x ) << '\t' << yAxis->GetBinCenter( y ) << '\t' << eff << std::endl; 

  }

  delete [] xSum;
}



void efficiency::write( TH2F * plot, std::ofstream & outfile, int inCh, int outCh,
			const TString & ref1, const TString & ref2 )
{

  // TAxis * xAxis = plot->GetXaxis();
  TAxis * yAxis = plot->GetYaxis();

  size_t xBins = plot->GetNbinsX();
  size_t yBins = plot->GetNbinsY();
//   double binWidth = plot->GetYaxis()->GetBinWidth(1);
//   double binStart = plot->GetXaxis()->GetXmin();
  // std::cout << binWidth << '\t' << binStart << '\t' << xBins << '\t' <<yBins << std::endl;

  for ( size_t x = 1; x <= xBins; ++x ) {
    double xSum = 0;
    for ( size_t y = 1; y <= yBins; ++y ) {
      xSum += plot->GetBinContent( x, y );
    }

    double eff = 0;
    double effPrev = 0;
    size_t y = 1;
    for ( ; y <= yBins; ++y ) {
      double binContent = plot->GetBinContent( x, y );
      if ( ! binContent ) continue;
      eff += binContent / xSum;
      if ( eff >= target_ ) {
	if ( target_ - effPrev < eff - target_ ) {
	  eff = effPrev;
	  --y;
	}
	break;
      }
      effPrev = eff;
    }
    // std::cout << binWidth << '\t' << binStart << '\t' << xBins << '\t' << binStart + binWidth*x  << '\t' << xSum << std::endl;
    // std::cout << plot->GetName() << '\t' << plot->GetBinCenter( x ) << '\t' << xSum << '\t' << effPrev << '\t' << eff << std::endl; 
    testfile_ << plot->GetName() << '\t' << plot->GetBinCenter( x ) << '\t' << eff << std::endl; 

    outfile << inCh << '\t' << outCh << '\t' << ref1 << '\t' << ref2 << '\t'
      // << xAxis->GetBinCenter( x ) << '\t' << yAxis->GetBinCenter( y ) << '\t' << eff << std::endl; // use bin center instead of bin lowedge
	    <<plot->GetBinLowEdge(x) << '\t' << yAxis->GetBinCenter( y ) << '\t' << eff << std::endl; 
  }


}


void efficiency::run( bool byTh )
{


  int wheels[6] = { -3, -2, -1, 1, 2, 3 };
  TString params[3] = { "Phi", "PhiBendIn","PhiBendOut" };
  TString ref[4] = { "DTCORR", "DTDIR", "DTIN", "DTOUT" };
  int inCh[2] = { 1, 2 };
  int outCh[2] = { 1, 2 };


  for ( int sector = 0; sector < 12; ++sector ) {
    for ( size_t w = 0; w < 6; ++w ) {
      for ( size_t p = 0; p < 3; ++p ) {

	TString fname = Form( "../Luts/%s/%sWh%dSc%d", outdir_.Data(), params[p].Data(), wheels[w], sector );
	std::ofstream outfile( fname.Data() );
      
	for ( size_t rf1 = 0; rf1 < 4; ++rf1 ) {
	  for ( size_t rf2 = 0; rf2 < 4; ++rf2 ) {
	    for ( size_t ch1 = 0; ch1 < 2; ++ch1 ) {
	      for ( size_t ch2 = 0; ch2 < 2; ++ch2 ) {

	

		TString name = Form( "L1ITMuPtPlotter/Wh%dSc0/hd%svsPtWh%dSc0inCh%doutCh%d%s%s",
				     wheels[w], params[p].Data(), wheels[w], inCh[ch1], outCh[ch2], ref[rf1].Data(), ref[rf2].Data() );  //so far use sector 0 for every sectors
		
		//if(rf1==1) std::cout<<" name "<<name.Data()<<std::endl;
		TH2F * plot = (TH2F *)(file_->Get( name ));

		// std::cout <<  w << '\t' << wheels[w] << '\t' << name.Data() << std::endl;
		if ( plot && byTh ) write2( plot, outfile, inCh[ch1], outCh[ch2], ref[rf1].Data(), ref[rf2].Data() );
		else if ( plot ) write( plot, outfile, inCh[ch1], outCh[ch2], ref[rf1].Data(), ref[rf2].Data() );
	      }
	    }
	  }
	}
      }
    }
  }

}





int main(int argc, char** argv) {


  if ( argc < 3 ) {
    std::cout << "Error in number of arguments: " << argc << std::endl;
    std::cout << "Passed args: " << argc << std::endl;
    for ( int i = 1; i < argc; ++i ) {
      std::cout << "\t" << argv[i] << std::endl;
    }
    std::cout << "Usage: \n\t\t " <<  argv[0] << " <input> <outdir> <log_file[=efflog.log]> <cut[=0.9]> <byTh[=0]"
              << std::endl << std::endl;
    std::cout << "\tes. : " << argv[0] << " /data2/battilan/L1Trigger/BarrelDTTFUpgrade/Plots/L1ITMuonBarrelPtPlots.root efficiency efflog.log 0.9 0" << std::endl;
    return -1;
  }

  TString infile  = argv[1];
  TString outdir  = argv[2];
  TString outfile = "efflog.log";
  float threshold = 0.9;
  bool byTh = 0;
  if ( argc > 3 ) outfile  = argv[3];
  if ( argc > 4 ) threshold = atof( argv[4] );
  if ( argc > 5 ) byTh  = atoi( argv[5] );

  efficiency eff( infile, outfile, outdir, threshold );
  eff.run( byTh );

  std::cout << "\nYou run with: \n\tinfile:  " << infile << "\n\toutdir:   " << outdir << "\n\toutfile:  " << outfile
	    << "\n\tthreshold: " << threshold << "\n\tbyThr:    " << (byTh ? 1 : 0 ) << std::endl;
}
