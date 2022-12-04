
// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "QbbFit.h"
#include "TF1.h"
#include "TFile.h"
#include <BAT/BCMath.h>
void QbbFit::MakePlots(TCanvas *&c, TCanvas *&cK, std::vector<double>pars)
{
  
  double S  = pars[0];
  double B  = pars[1];
  double p1 = pars[2];
  double p2=pars[3];
  double Q  = pars[4];
  std::cout<<"best fit pars are : "<<std::endl;
  std::cout<<"S                 : "<<S<<std::endl;
  std::cout<<"B                 : "<<B<<std::endl;
  std::cout<<"p1                : "<<p1<<std::endl;
  std::cout<<"Q                 : "<<Q<<std::endl;


  TH1D *hKurieResid = (TH1D*)fData->Clone("hKurieResid");
  TH1D *hKurie =(TH1D*)fData->Clone("hKurie");
  TGraph *gK = new TGraph();
  TH1D *hPlotSig=(TH1D*)fData->Clone("hPlotSig");
  TH1D *hPlotBkg =(TH1D*)fData->Clone("hPlotBkg");
  TH1D *hPlotModel=(TH1D*)fData->Clone("hPlotModel");
  int counter=0;
  TGraphErrors *gResid = new TGraphErrors();
  hPlotModel->SetBinErrorOption(TH1::kPoisson);
  fData->SetBinErrorOption(TH1::kPoisson);
  hKurie->SetLineColor(kBlack);
  hKurieResid->SetLineColor(kBlack);
  
  

  for (int i=fData->FindBin(0);i<fData->FindBin(4000);i++)
    {

      //      p1=0;
      if (fFix)
	{
	  p1=0;
	}
      double energy =fData->GetBinCenter(i);
      fEfficiency->SetParameter(0,p1);
      fEfficiencyBkg->SetParameter(0,p1);
      double B_pred = fEfficiencyBkg->Eval(energy)*B*fBkg->GetBinContent(fBkg->FindBin(energy));
      hPlotBkg->SetBinContent(i,B_pred);
      // get the predicted signal                                                                                                                                                                          
      //fEfficiency->SetParameter(0,p1);
      double S_pred =fBinning*S*fShape->Eval(energy)*fEfficiency->Eval(energy)*pow(Q-energy,5);
      if (S_pred<0)
	S_pred=0;
      
      hPlotSig->SetBinContent(i,S_pred);
      
      
      hPlotModel->SetBinContent(i,S_pred+B_pred);


      // FILL THE KURIE PLOT

      double value = fData->GetBinContent(i)-B_pred;
      std::cout<<"value "<<pow(value/(fBinning*fShape->Eval(energy)),0.2)<<std::endl;
      std::cout<<fShape->Eval(energy)<<std::endl;
      std::cout<<"energy = "<<energy<<std::endl;
      if (value>0 &&fShape->Eval(energy)>0)
        hKurie->SetBinContent(i,pow(value/(fBinning*fShape->Eval(energy)),0.2));
      else if (value<0 &&fShape->Eval(energy)>0)
        hKurie->SetBinContent(i,-pow(-value/(fBinning*fShape->Eval(energy)),0.2));
      else
        hKurie->SetBinContent(i,0);

      if (value>0 &&fShape->Eval(energy)>0)
        hKurie->SetBinError(i,0,0.2*pow((value/(fBinning*fShape->Eval(energy))),-0.8)*fData->GetBinError(i)/(fBinning*fShape->Eval(energy)));
      else if (value<0 &&fShape->Eval(energy)>0)
        hKurie->SetBinError(i,0,0.2*pow(-(value/(fBinning*fShape->Eval(energy))),-0.8)*fData->GetBinError(i)/(fBinning*fShape->Eval(energy)));
      else
	hKurie->SetBinError(i,0,0);
      double resid = (hKurie->GetBinContent(i)-(Q-energy))/(hKurie->GetBinError(i));
      if (fabs(resid)<200)
	{
	  hKurieResid->SetBinContent(i,resid);
	  hKurieResid->SetBinError(i,0,1);
	}
      else
	{
	  hKurieResid->SetBinContent(i,0);
          hKurieResid->SetBinError(i,0,0);

	}
      // now get the errors up and down
      double eUp = fData->GetBinErrorUp(i);
      double eDown=fData->GetBinErrorLow(i);
      
      //      std::cout<<"Enrgy  = "<<energy<<" VALUE = "<<fData->GetBinContent(i)<<" + / -  "<<eUp<<" , "<<eDown<<std::endl;
      if (fData->GetBinContent(i)!=0)
	{
	  
	  if (fData->GetBinContent(i)>=S_pred+B_pred)
	    {
	      
	      gResid->SetPoint(counter,energy,(fData->GetBinContent(i)-S_pred-B_pred)/eUp);
	      gResid->SetPointError(counter,0,1);
	    }
	  else
	    {
	      gResid->SetPoint(counter,energy,(fData->GetBinContent(i)-S_pred-B_pred)/eDown);
	      gResid->SetPointError(counter,0,1);
	    }
	  counter++;
	}
      else
	{
	  gResid->SetPoint(counter,energy,0); //(fData->GetBinContent(i)-S_pred-B_pred)/eDown);
	  gResid->SetPointError(counter,0,0);
	  counter++;
	}
    }


  


  // lets make two plots a spectrum and a Kurie Plot

  c->cd();
  c->SetBottomMargin(0.2);
  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  fData->SetTitle(" ; Energy [keV] ; counts/bin ; ");
  //fData->GetYaxis()->SetRangeUser(0.0000000000001,10000000);
  fData->GetXaxis()->SetRangeUser(00,3200);
  fData->GetYaxis()->SetRangeUser(0.01,2000);
  fData->Draw("E");
  hPlotModel->SetLineColor(2);
  hPlotModel->Draw("sameHIST");
  hPlotSig->SetLineColor(8);
  hPlotBkg->SetLineColor(9);
  hPlotBkg->SetLineWidth(1);
  hPlotSig->SetLineWidth(1);
  fData->SetLineWidth(1);
  fData->SetLineColor(kBlack);
  hPlotModel->SetLineWidth(1);
  hPlotBkg->Draw("HISTsame");
  hPlotSig->Draw("HISTsame");
  l->AddEntry(fData,"Data");
  l->AddEntry(hPlotModel,"Reconstruction");
  l->AddEntry(hPlotBkg,"Background");
  l->AddEntry(hPlotSig,"Signal");
  l->Draw("same");
  c->SetLogy(1);


  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);

  pad->SetTopMargin(0.8);
  pad->Draw();
  pad->SetFillStyle(0);
  pad->cd();
  gResid->GetYaxis()->SetLabelSize(0.02);
  gResid->SetTitle("; Energy [keV] ; Residual ;  "); 
  TF1 *fConst = new TF1("fConst","0",0,3200);
  gResid->Draw("APE");
  gResid->GetXaxis()->SetRangeUser(000,3200);
  gResid->GetHistogram()->SetMaximum(5);
  gResid->GetHistogram()->SetMinimum(-5);
  fConst->Draw("same");

  cK->cd();
  //cK->SetBottomMargin(0.2);
  TLegend *l2 = new TLegend(0.7,0.6,0.9,0.9);
  hKurie->SetTitle(" ; Energy [keV] ; K(E) ; ");
  //fData->GetYaxis()->SetRangeUser(0.0000000000001,10000000);

  TF1 *fFit=  new TF1("fFit","(x<[0])*([0]-x)",0,4000);
  fFit->SetParameter(0,Q);
  fFit->SetNpx(1000);
  cK->SetLogy(0);
  hKurie->GetXaxis()->SetRangeUser(00,3200);
  hKurie->Draw("E");
  fFit->SetLineColor(2);
  fFit->Draw("sameC");
  hKurie->SetLineWidth(1);
  hKurie->SetLineColor(kBlack);
  l2->AddEntry(hKurie,"Data");

  l2->AddEntry(fFit,"Reconstruction");
  l2->Draw("same");
  cK->SetLogy(0);
  /*
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);

  pad2->SetTopMargin(0.8);
  pad2->Draw();
  hKurieResid->SetLineWidth(1);
  hKurieResid->SetLineColor(kBlack);
  pad2->SetFillStyle(0);
  pad2->cd();
  hKurieResid->GetYaxis()->SetLabelSize(0.02);
  hKurieResid->SetTitle("; Energy [keV] ; Residual ;  ");
  hKurieResid->Draw("E");
  hKurieResid->GetXaxis()->SetRangeUser(000,3200);
  hKurieResid->GetYaxis()->SetRangeUser(-5,5);

  fConst->Draw("same");
  */
  cK->cd();

}
  
  // ---------------------------------------------------------
QbbFit::QbbFit(const std::string& name,TString path_shape,TString path_bkg,bool fix)
    : BCModel(name)
{


  fFix=fix;
  // get the inputs

  TFile *f_shape = new TFile(path_shape);
  fShape =(TF1*)f_shape->Get("fu");
  f_shape->ls();
  TFile *f_bkg = new TFile(path_bkg);
  fBkg = (TH1D*)f_bkg->Get("Background_model");
  fData=(TH1D*)f_bkg->Get("Experimental_data");
  std::cout<<fData->GetNbinsX()<<std::endl;
  f_bkg->ls();
  fEfficiency =new TF1("fEfficiency","1+[0]*x/3034.",0,4000);
  fEfficiencyBkg =new TF1("fEfficiencyBkg","1+[0]*x/3034.",0,4000);

   
  // some defaults
  fUpper=3500.;
  fLower=100.;
  fCounter=0;
  fBinning=1;

  fBkg->Rebin(fBinning);
  fData->Rebin(fBinning);
  
  // Define parameters here in the constructor.
  // we have parameters
  // 1) Normalisation of signal from [0,2]
  // 2) Normalisation of bkg again [0,2]
  // 3) Efficiency fraction from [0.9,1.1]
  // 4) Qbb from [3034-100,3034+100]

  AddParameter("fS",0.9,1.1,"f_{S}"," ");
  AddParameter("fB",0.9,1.1,"f_{B}"," ");
  AddParameter("p1",-0.2,0.2,"p1"," ");
  AddParameter("p2",-0.2,0.2,"p2"," ");

  AddParameter("Qbb",3034-20,3034+20,"Q_{#beta#beta}","[keV]");


  GetParameters().SetPriorConstantAll();
  SetPriorGauss("p1",0.0,0.03,0.03);
  GetParameters().SetNBins(500);

  
}

// ---------------------------------------------------------
QbbFit::~QbbFit()
{
    // destructor
}

// ---------------------------------------------------------
double QbbFit::LogLikelihood(const std::vector<double>& pars)
{

  // the likelihood if a simple binned one
  int verbose;
  double logL=0;
  if (fCounter%10000000==0)
    {
      verbose=1;
    }
  else
    {
      verbose=0;
    }

  double S  = pars[0];
  double B  = pars[1];
  double p1 = pars[2];
  double p2=pars[3];
  double Q  = pars[4];
  for (int i=fData->FindBin(fLower);i<fData->FindBin(fUpper);i++)
    {;
      // get th
      double energy =fData->GetBinCenter(i);
      double N = fData->GetBinContent(i);
      //p1=0;
      if (fFix)
	{
	  p1=0;
	}
	  
      fEfficiency->SetParameter(0,p1);
      fEfficiencyBkg->SetParameter(0,p1);
      // get the predicted background
      double B_pred = fEfficiencyBkg->Eval(energy)*B*fBkg->GetBinContent(fBkg->FindBin(energy));
      
      // get the predicted signal
      //fEfficiency->SetParameter(0,p1);
      double S_pred =fBinning*S*fShape->Eval(energy)*fEfficiency->Eval(energy)*pow(Q-energy,5);
      if (energy>Q)
	{
	  S_pred=0;
	}
		    
      double lambda = B_pred+S_pred;
      if( lambda <= 0. ) continue;
      logL += -lambda + N * log( lambda ) - BCMath::LogFact( N );

      if (verbose &&i%10000==0)
	{
          std::cout<<"energy = "<<energy<<std::endl;
          std::cout<<"N      = "<<N<<std::endl;
          std::cout<<"b      = "<<B_pred<<std::endl;
          std::cout<<"sig    = "<<S_pred<<std::endl;
          std::cout<<"lambda = "<<lambda<<std::endl;
          std::cout<<"logL   = "<<logL<<std::endl;
          std::cout<<" "<<std::endl;
        }


    }
  fCounter++;
  return logL;

  
}

// ---------------------------------------------------------
// double QbbFit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
// void QbbFit::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }
