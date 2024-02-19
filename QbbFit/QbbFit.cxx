
// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "QbbFit.h"
#include "TF1.h"
#include "TH3D.h"
#include "TFile.h"
#include <BAT/BCMath.h>

double QbbFit::Norm(double Q)
{
  // get the normalisation of int A(E)*(Qnom-E)^5/ int A(E)*(Q-E)^5

  double QDef=3034.4;
  double NormDef;
  double Norm;
  if (fIsImp)
    {
      NormDef =fIntegrals[0]->Eval(QDef)+
        fIntegrals[2]->Eval(QDef)*xi31+
        xi31*xi31*fIntegrals[22]->Eval(QDef)/3.+
        (fIntegrals[4]->Eval(QDef)*(xi31*xi31/3.+xi51));
    }
  else
    {
      NormDef=fIntegral->Eval(QDef);
    }
  if (fIsImp)
    {
      Norm =fIntegrals[0]->Eval(Q)+
        fIntegrals[2]->Eval(Q)*xi31+
        xi31*xi31*fIntegrals[22]->Eval(Q)/3.+
        (fIntegrals[4]->Eval(Q)*(xi31*xi31/3.+xi51));
    }
  else
    {
      Norm=fIntegral->Eval(Q);
    }
  return NormDef/Norm;
}

  
void QbbFit::MakePlots(TCanvas *&c, TCanvas *&cK, std::vector<double>pars)
{
  
 

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
  
  double Q=pars[fParameterMap["Qbb"]];


  for (int i=fData->FindBin(0);i<fData->FindBin(4000);i++)
    {
      
      double energy =fData->GetBinCenter(i);

      double S_pred,B_pred;
      GetSignalBackground(energy,pars,S_pred,B_pred);
      hPlotSig->SetBinContent(i,S_pred);
      
      
      hPlotModel->SetBinContent(i,S_pred+B_pred);


      // FILL THE KURIE PLOT
      double value = fData->GetBinContent(i)-B_pred;
      if (value>0 &&GetSpectrum(energy)>0)
        hKurie->SetBinContent(i,pow(value/(fBinning*GetSpectrum(energy)),0.2));
      else if (value<0 &&GetSpectrum(energy)>0)
        hKurie->SetBinContent(i,-pow(-value/(fBinning*GetSpectrum(energy)),0.2));
      else
        hKurie->SetBinContent(i,0);

      if (value>0 &&GetSpectrum(energy)>0)
        hKurie->SetBinError(i,0,0.2*pow((value/(fBinning*GetSpectrum(energy))),-0.8)*fData->GetBinError(i)/(fBinning*GetSpectrum(energy)));
      else if (value<0 &&GetSpectrum(energy)>0)
        hKurie->SetBinError(i,0,0.2*pow(-(value/(fBinning*GetSpectrum(energy))),-0.8)*fData->GetBinError(i)/(fBinning*GetSpectrum(energy)));
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

  TLegend *l2 = new TLegend(0.7,0.6,0.9,0.9);
  hKurie->SetTitle(" ; Energy [keV] ; K(E) ; ");
  
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
  
  cK->cd();

}
void QbbFit::SetRatioPrior(double ratio,double error)
{
  fIsRatioPrior=true;

  fXiRatio=ratio;
  fXiRatioError=error;

}
  
void QbbFit::GetSignalBackground(double &energy,const std::vector<double> &pars,double &S_pred,double &B_pred)
{

  double S  = pars[fParameterMap["fS"]];
  double B  = pars[fParameterMap["fB"]];
  // double xi31,xi51;
  if (fIsImp)
    {
      xi31=pars[fParameterMap["xi31"]];
      xi51=pars[fParameterMap["xi51"]];
    }
  double p1;
  if (!fFix)
    p1=pars[fParameterMap["p1"]];
  double LYp0=1,LYp1=1,PCAp0=1,PCAp1=1;
  double constant_eff=1;
  if (fFullEfficiency)
    {
      LYp0=pars[fParameterMap["LYp0"]];
      LYp1=pars[fParameterMap["LYp1"]];
      PCAp0=pars[fParameterMap["PCAp0"]];
      PCAp1=pars[fParameterMap["PCAp1"]];
      constant_eff=pars[fParameterMap["constant_eff"]];
    }
  double Q  = pars[fParameterMap["Qbb"]];
  S*=Norm(Q);
  if (fFix)
    {
      p1=0;
    }
  
  double energyBias;
  if (fFloatBias==true)
    {
      
      fBiasp0=pars[fParameterMap["Biasp0"]];
      fBiasp1=pars[fParameterMap["Biasp1"]];
      fBiasp2=pars[fParameterMap["Biasp2"]];
      
      

      fEnergyBias->SetParameters(fBiasp0,fBiasp1,fBiasp2);
      energyBias=fEnergyBias->Eval(energy);
    }
  else
    {
      energyBias=0;
    }
  // set the TF1 values                                                                                                                                                                                         
  if (!fFix)
    {
      fEfficiency->SetParameter(0,p1);
    }
  if (fFullEfficiency==true)
    {
	  
      fLYEff->SetParameters(LYp0,LYp1);
      fPCAEff->SetParameters(PCAp0,PCAp1);
    }

  double efficiency_S;
  double efficiency_B;
  if (fFullEfficiency==false)
    {
      efficiency_S=GetEfficiency(energy,constant_eff,1);
      efficiency_B=efficiency_S;
    }
  else
    {
      efficiency_S=GetEfficiency(energy,constant_eff,0);
      efficiency_B=GetEfficiency(energy,constant_eff,1);
    }
  B_pred = efficiency_B*B*fBkg->GetBinContent(fBkg->FindBin(energy));
  
  energy=-energyBias+energy;
  
  S_pred=fBinning*S*GetSpectrum(energy)*efficiency_S*pow(Q-energy,5);
 
  if (energy>Q)
    {
      
      S_pred=0;
    }
  
}

double QbbFit::GetSpectrum(double &energy)
{
  if (fIsImp)
    {
      return fShapes[0]->Eval(energy)+
	fShapes[2]->Eval(energy)*xi31+
	xi31*xi31*fShapes[22]->Eval(energy)/3.+                                                                              
	(fShapes[4]->Eval(energy)*(xi31*xi31/3.+xi51));
    }
  return fShape->Eval(energy);
}

double QbbFit::GetEfficiency(double &energy,double &constant_eff,bool &simple)
{

  double eff=1;
  if (simple &&!fFix)
    {
      // the basic model of non-constant efficency
      eff=fEfficiency->Eval(energy);

    }
  else if(!simple &&!fFix)
    {
      // we have to divide the MC by the value in MC
      eff/=fEffValMC;
      eff*=constant_eff;
      
      eff*=fLYEff->Eval(energy);
      eff*=fPCAEff->Eval(energy);
      if (fCounter%1000000==0)
	{
	  std::cout<<"eff = "<<eff<<std::endl;
	  std::cout<<"constant_eff = "<<constant_eff<<std::endl;
	  std::cout<<"LY_eff = "<<fLYEff->Eval(energy)<<std::endl;
	  std::cout<<"energy = "<<energy<<std::endl;
	  std::cout<<"PCA_eff = "<<fPCAEff->Eval(energy)<<std::endl;
	  std::cout<< " "<<std::endl;
	}
      fCounter++;
    }
  return eff;
}

TF1 * QbbFit::GetIntegralFunction(TF1 *&f,TString name)
{
  TF1 *fout = new TF1(name,"pow(x,6)*(4620*[0]+660*[1]*x+165*[2]*pow(x,2)+55*[3]*pow(x,3)+22*[4]*pow(x,4)+10*[5]*pow(x,5)+5*[6]*pow(x,6))",0,4000);

  
  fout->SetParameters(f->GetParameters());
  return fout;
}
  // ---------------------------------------------------------
QbbFit::QbbFit(const std::string& name,TString path_shape,TString path_bkg,TString path_bias,bool fix,bool floatbias,bool fullefficiency,double effvalMC,bool imp)
    : BCModel(name)
{

  fIsRatioPrior=0;
  fFix=fix;
  // get the inputs
  fFloatBias=floatbias;
  fFullEfficiency=fullefficiency;
  fEffValMC=effvalMC;
  fIsImp=imp;
  if (fFullEfficiency)
    {
      fLYEff=new TF1("fLYEff","pol1",0,4000);
      fPCAEff=new TF1("fPCAEff","pol1",0,4000);
    }

  TFile *f_shape = new TFile(path_shape);

  
  if (!imp)
    {
      fShape =(TF1*)f_shape->Get("fu");
      fIntegral=(TF1*)GetIntegralFunction(fShape,"fint");
    }
  else
    {
      f_shape->ls();


      fShapes[0]=(TF1*)f_shape->Get("fu_h0");
      fShapes[2]=(TF1*)f_shape->Get(Form("fu_h2"));
      fShapes[22]=(TF1*)f_shape->Get(Form("fu_h22"));
      fShapes[4]=(TF1*)f_shape->Get(Form("fu_h4"));

      fIntegrals[0]=(TF1*)GetIntegralFunction(fShapes[0],"f_int_0");
      fIntegrals[2]=(TF1*)GetIntegralFunction(fShapes[2],"f_int_2");
      fIntegrals[22]=(TF1*)GetIntegralFunction(fShapes[22],"f_int_22");
      fIntegrals[4]=(TF1*)GetIntegralFunction(fShapes[4],"f_int_4");
    }
  TFile *f_bkg = new TFile(path_bkg);
  fBkg = (TH1D*)f_bkg->Get("Background_model");
  fData=(TH1D*)f_bkg->Get("Experimental_data");
  std::cout<<fData->GetNbinsX()<<std::endl;
  f_bkg->ls();
  fEfficiency =new TF1("fEfficiency","1+[0]*x/3034.",0,4000);
  fEfficiencyBkg =new TF1("fEfficiencyBkg","1+[0]*x/3034.",0,4000);

  TFile *f_bias =new TFile(path_bias);
  fBias =(TH3D*)f_bias->Get("h_bias");

  fp1Sigma=0.03;
  fEnergyBias=new TF1("fEnergyBias","pol2",0,4000);
   
  // some defaults
  fUpper=3500.;
  fLower=100.;
  fCounter=0;
  fBinning=10;

  
  // Define parameters here in the constructor.
  // we have parameters
  // 1) Normalisation of signal from [0,2]
  // 2) Normalisation of bkg again [0,2]
  // 3) Efficiency fraction from [0.9,1.1]
  // 4) Qbb from [3034-100,3034+100]
  int counter=0;
  AddParameter("fS",0.5,1.5,"f_{S}"," ");
  fParameterMap["fS"]=counter;
  counter++;

  if (imp)
    {
      AddParameter("xi31",0,2,"#xi_{3,1}"," ");
      fParameterMap["xi31"]=counter;
      counter++;
      
      AddParameter("xi51",0,2,"#xi_{5,1}"," ");
      fParameterMap["xi51"]=counter;
      counter++;
    }
  AddParameter("fB",0.5,1.5,"f_{B}"," ");
  fParameterMap["fB"]=counter;
  counter++;

  if (!fFix)
    {
      AddParameter("p1",-0.2,0.2,"p1"," ");
      fParameterMap["p1"]=counter;
      counter++;
    }
  AddParameter("Qbb",3034-20,3034+20,"Q_{#beta#beta}","[keV]");
  fParameterMap["Qbb"]=counter;
  counter++;


  if (fFloatBias==true)
    {
      AddParameter("Biasp0",-2,2,"biasp0"," ");
      fParameterMap["Biasp0"]=counter;
      counter++;
      
      AddParameter("Biasp1",-2e-3,1e-3,"biasp1"," ");
      fParameterMap["Biasp1"]=counter;
      counter++;
  
      AddParameter("Biasp2",-5e-7,7e-7,"biasp2"," ");
      fParameterMap["Biasp2"]=counter;
      counter++;

    }

  
  

 
  if (fFullEfficiency==true)
    {
      AddParameter("PCAp0",.9,1.02,"PCAp0"," ");
      fParameterMap["PCAp0"]=counter;
      counter++;

      AddParameter("PCAp1",-7e-5,7e-5,"PCAp1"," ");
      fParameterMap["PCAp1"]=counter;
      counter++;

      AddParameter("LYp0",.9,1.02,"LYp0"," ");
      fParameterMap["LYp0"]=counter;
      counter++;

      AddParameter("LYp1",-7e-5,7e-5,"LYp1"," ");
      fParameterMap["LYp1"]=counter;
      counter++;

      AddParameter("constant_eff",0.9,1,"constant_eff"," ");
      fParameterMap["constant_eff"]=counter;
      counter++;

    }

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
  int counter=0;

  
  for (int i=fData->FindBin(fLower);i<fData->FindBin(fUpper);i++)
    {
      // get th
      double energy =fData->GetBinCenter(i);
      double N = fData->GetBinContent(i);
      double S_pred,B_pred;
      GetSignalBackground(energy,pars,S_pred,B_pred);
      double lambda = B_pred+S_pred;
      counter++;
      if ( fCounter%10000==0 &&counter==200 )
        {
          std::cout<<"energy = "<<energy<<std::endl;
          std::cout<<"N      = "<<N<<std::endl;
          std::cout<<"b      = "<<B_pred<<std::endl;
          std::cout<<"sig    = "<<S_pred<<std::endl;
	  std::cout<<"xi     = "<<xi31<<" , "<<xi51<<std::endl;
	  std::cout<<"lambda = "<<lambda<<std::endl;
          std::cout<<"logL   = "<<logL<<std::endl;
          std::cout<<" "<<std::endl;
        }

      if( lambda <= 0. ) continue;
      logL += -lambda + N * log( lambda ) - BCMath::LogFact( N );

      
      
     
    }
  fCounter++;
  return logL;

  
}

// ---------------------------------------------------------
 double QbbFit::LogAPrioriProbability(const std::vector<double>& pars)
 {
   // return the log of the prior probability p(pars)
   // If you use built-in priors, leave this function commented out.
   
   // We have up to 8 parameters

   double logPrior=0;

   return 1;
   // 0,1,2,3 have uniform
   double low,high;
   std::vector<std::string> par_name{"fS","fB","Qbb"};
   

   for (auto & par:par_name)
     {
       low =GetParameter(fParameterMap[par]).GetLowerLimit();
       high =GetParameter(fParameterMap[par]).GetUpperLimit();
       if (pars[fParameterMap[par]]>low &&pars[fParameterMap[par]]<high)
	 logPrior+=log(1/(high-low));
       else
	 logPrior-=1e20;
     }
 
       
   double mean=0;
   double sigma=fp1Sigma;
   // gaussian prior on par 4
   //if (fFix==0)
   //logPrior+=log(1/(sqrt(2*3.14)*sigma))-pow((pars[fParameterMap["p1"]]-mean)/sigma,2)/0.5;

   if (fIsRatioPrior==1)
     logPrior+=log(1/(sqrt(2*3.14)*fXiRatioError))-pow((pars[fParameterMap["xi51"]]/pars[fParameterMap["xi31"]]-fXiRatio)/fXiRatioError,2)/0.5;
   
   if (fFullEfficiency)
     {
       double LYp0=pars[fParameterMap["LYp0"]];
       double LYp1=pars[fParameterMap["LYp1"]];
       double PCAp0=pars[fParameterMap["PCAp0"]];
       double PCAp1=pars[fParameterMap["PCAp1"]];
       double constant_eff=pars[fParameterMap["constant_eff"]];
       
       double probLY = fLYEffHisto->GetBinContent(fLYEffHisto->FindBin(LYp1,LYp0));
       double probPCA=(fPCAEffHisto->GetBinContent(fPCAEffHisto->FindBin(PCAp1,PCAp0)));
       sigma =fConstantEffSigma;
       mean=fConstantEffMean;
       double prob_constant=log(1/(sqrt(2*3.14)*sigma))-pow((constant_eff-mean)/sigma,2)/0.5;
     
       if (probLY!=0)
	{
          logPrior+=log(probLY);
        }
      else
        logPrior-=1e20;

       if (probPCA!=0)
	 {
	   logPrior+=log(probPCA);
        }
       else
	 logPrior-=1e20;

       logPrior+=prob_constant;
     
       if (fCounter%100000==0)
	 {
	   std::cout<<"p(const) = "<<prob_constant<<std::endl;
	   std::cout<<"p(PCA) = "<<probPCA<<std::endl;
	   std::cout<<"p(lY) = "<<probLY<<std::endl;
	   std::cout<<"pars PCA = "<<PCAp0<<" , "<<PCAp1<<std::endl;
	   std::cout<<"pars LY = "<<LYp0<<" , "<<LYp1<<std::endl;
	   std::cout<<" logPrior = "<<logPrior<<std::endl;
	 }
     }
   if (fFloatBias==true)
    {
      int p2 = fParameterMap["Biasp2"];
      int p0=fParameterMap["Biasp0"];
      int p1=fParameterMap["Biasp1"];
      double prob = fBias->GetBinContent(fBias->FindBin(pars[p2],pars[p1],pars[p0]));

      if (prob!=0)
	{
	  logPrior+=log(prob);
	}
      else
	logPrior-=1e20;
    
    }
   return logPrior;
  



 }

// ---------------------------------------------------------
// void QbbFit::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }
