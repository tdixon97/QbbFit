// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__QBBFIT__H
#define __BAT__QBBFIT__H

#include <BAT/BCModel.h>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include <string>
#include "TPad.h"
#include "TGraphErrors.h"
#include <vector>
#include "TCanvas.h"
// This is a QbbFit header file.
// Model source code is located in file QbbFit/QbbFit.cxx

// ---------------------------------------------------------
class QbbFit : public BCModel
{

public:

  // Constructor
  QbbFit(const std::string& name,TString path_shape,TString path_bkg,bool fix);
  void MakePlots(TCanvas *&c, TCanvas *&cK, std::vector<double>pars);
  void SetUpperLower(double low,double up){fUpper=up;fLower=low;};
  double GetUpper(){return fUpper;};
  double GetLower(){return fLower;};
  
  // Destructor
  ~QbbFit();
  
  // Overload LogLikelihood to implement model
  double LogLikelihood(const std::vector<double>& pars);
  bool fFix;
  // Overload LogAprioriProbability if not using built-in 1D priors
  // double LogAPrioriProbability(const std::vector<double> & pars);
  
  // Overload CalculateObservables if using observables
  // void CalculateObservables(const std::vector<double> & pars);
private:
  int fCounter;
  TH1D*fBkg;
  TH1D *fData;
  TF1 *fShape;
  TF1 *fEfficiencyBkg;
  TF1 *fEfficiency;
  double fUpper,fLower;
  int fBinning;
};
// ---------------------------------------------------------

#endif
