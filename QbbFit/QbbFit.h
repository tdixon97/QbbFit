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
#include "TH2D.h"
#include "TH3D.h"
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
  QbbFit(const std::string& name,TString path_shape,TString path_bkg,TString path_bias,bool fix,bool floatbias,bool fullefficiency,double effvalMC,bool is_imp);
  void MakePlots(TCanvas *&c, TCanvas *&cK, std::vector<double>pars);
  void SetUpperLower(double low,double up){fUpper=up;fLower=low;};
  double GetUpper(){return fUpper;};
  double GetLower(){return fLower;};
  double GetSpectrum(double &energy);
  TF1 * GetIntegralFunction(TF1 *&f,TString name);
  double Norm(double Q);
  void SetRatioPrior(double ratio,double error);

  // Destructor
  ~QbbFit();
  
  // Overload LogLikelihood to implement model
  double LogLikelihood(const std::vector<double>& pars);
  bool fFix;
  // Overload LogAprioriProbability if not using built-in 1D priors
  double LogAPrioriProbability(const std::vector<double> & pars);
  double GetEfficiency(double &energy,double & constant_eff,bool simple);
  void GetSignalBackground(double &energy,const std::vector<double> &pars,double &S_pred,double &B_pred);
  void SetConstantEff(double &mean,double &sigma){fConstantEffMean=mean; fConstantEffSigma=sigma;};
  void SetThreshold(double t){fLower=t;};
  void SetBinning(int b){fBinning=b;  fBkg->Rebin(fBinning);fData->Rebin(fBinning);};
  void SetLYPCA(TString LY,TString PCA)
  {
    TFile *fPCA = new TFile(PCA); fPCAEffHisto=(TH2D*)fPCA->Get("h_bias");fPCAEffHisto->RebinX(10);fPCAEffHisto->RebinY(10);
    
    TFile *fLY = new TFile(LY); fLYEffHisto=(TH2D*)fLY->Get("h_bias");fLYEffHisto->RebinX(10); fLYEffHisto->RebinY(10);
  }

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
  double fp1Sigma;
  int fBinning;
  bool fFloatBias;
  TH3D* fBias;
  bool fIsImp;
  std::map<int,TF1*>fShapes;
  std::map<int,TF1*>fIntegrals;
  TF1 *fIntegral;
  TF1 *fEnergyBias;
  double xi31,xi51;
  bool fIsRatioPrior;
  double fXiRatio;
  double fXiRatioError;
  std::map<std::string,int> fParameterMap;
  double fBiasp0,fBiasp1,fBiasp2;
  double fEffValMC;
  bool fFullEfficiency;
  TF1 *fLYEff;
  TF1 *fPCAEff;
  TH2D *fLYEffHisto;
  TH2D *fPCAEffHisto;
  double fConstantEffMean;
  double fConstantEffSigma;
};
// ---------------------------------------------------------

#endif
