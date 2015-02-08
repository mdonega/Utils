#include <TF1.h>
#include <TH1.h>
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"

// LOAD ROOFIT:
// {
//   gSystem->Load("libRooFit.so") ;
//   gSystem->Load("libRooFitCore.so") ;
// }

void fitCB(TH1 *h, double &mean, double &emean, double &chi2)
{
  bool dbg = 0;
  int nBinsErecEgen    = h->GetNbinsX();
  double ErecEgenMin   = h->GetBinLowEdge(1);
  double ErecEgenMax   = h->GetBinLowEdge(nBinsErecEgen+1);
  
  double startMean  = 0;
  double startSigma = 0;

  // temporary histogram to get CB sigma
  //
  TH1F *htmp = new TH1F("htmp","htmp",nBinsErecEgen,ErecEgenMin, ErecEgenMax);

  // protect the fits from small statistics rebinning
  //
  if (h->GetMaximum()< 50.) { h->Rebin(2); htmp->Rebin(2); }
  if (h->GetMaximum()< 25.) { h->Rebin(2); htmp->Rebin(2); }
  
  // get the right shoulder of the histogram
  double max = h->GetMaximum();
  for (int i=1; i<h->GetNbinsX(); ++i) {if (h->GetBinContent(i) > max*0.8) htmp->SetBinContent(i,h->GetBinContent(i));}
  for (int i=1; i<h->GetNbinsX(); ++i) {if (i > h->GetMaximumBin()) htmp->SetBinContent(i,h->GetBinContent(i)); }  
  
  // fit the right shoulder to a gaussian
  TF1 *gtmp = 0;  
  if (h->GetRMS() < 0.01) {
    // when it's very narrow there is no need to fit
    startMean  = htmp->GetMean();
    startSigma = h->GetRMS();      
  }
  else{      
    // fit to a gaussian
    gtmp = new TF1("gtmp","gaus(0)",0,2.0);    
    
    // use the mean and RMS of the Erec/Egen as initial parameters for the gaussian
    gtmp->SetParameters(10, htmp->GetMean(), htmp->GetRMS());

    // first symmetric fit in the 5-sigma range
    htmp->Fit("gtmp","","qsame",htmp->GetMean()-5*htmp->GetRMS(),htmp->GetMean()+5*htmp->GetRMS()); 

    // then refine range to catch the right shoulder (half sigma to the left, 5 sigmas to the right)
    htmp->Fit("gtmp","","qsame",
	      gtmp->GetParameter(1)-0.5*TMath::Abs(gtmp->GetParameter(2)),
	      gtmp->GetParameter(1)+5*TMath::Abs(gtmp->GetParameter(2)));

    startMean  =  gtmp->GetParameter(1);
    startSigma =  gtmp->GetParameter(2);

  }
  delete gtmp;
  delete htmp;
  
  // RooDataHist from input TH1F
  RooRealVar  x("x","x", ErecEgenMin, ErecEgenMax);  
  RooDataHist data("data","Ereco/Egen",x,h); 
 
  // Initialize CB parameters
  double  alphaStart =   0.5;
  double  alphaMin   =   0.1;
  double  alphaMax   =  10.0;

  double  nStart     =  50.0;
  double  nMin       =   0.0;
  double  nMax       = 200.0;
  
  // CrystalBall fit function
  RooRealVar alpha  ("alpha"  ,        "alpha" , alphaStart,      alphaMin,        alphaMax); 
  RooRealVar n      ("n"      ,            "n" ,     nStart,          nMin,            nMax); 
  RooRealVar cbmean ("cbmean" ,       "cbmean" , startMean ,   ErecEgenMin,     ErecEgenMax);
  RooRealVar cbsigma("cbsigma",      "cbsigma" , startSigma, startSigma*0.5, startSigma*1.5);  // constraint !

  RooCBShape cball  ("cball"  , "crystal ball" , x, cbmean, cbsigma, alpha, n);
  
  // Fit      
  RooFitResult* fitres =cball.fitTo(data,
				    RooFit::Range(0, 1.3),
				    RooFit::Minos(kFALSE),
				    RooFit::PrintEvalErrors(-1));    
  if (dbg) fitres->Print();
  mean  = cbmean.getVal();
  emean = cbmean.getError();
  
  RooPlot*    xframe = x.frame(RooFit::Name("xframe"),RooFit::Title("E^{reco} / E^{gen}")) ;
  data.plotOn(xframe);
  cball.plotOn(xframe);
  chi2 = xframe->chiSquare(4); // dof = 4 = number fof floating parameters

  return;
}
