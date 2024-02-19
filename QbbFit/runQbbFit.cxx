// ***************************************************************
// This file was created using the bat-project script
// for project QbbFit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include "QbbFit.h"
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include "TGraphErrors.h"
#include "TF1.h"
#include <tuple>
#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include <utility>
#include "TSpectrum.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>

void Usage()
{

  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./runQbbFit"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-N name of fit"<<std::endl;
  std::cout<<"-S theory file (SSD or HSD shape) "<<std::endl;
  std::cout<<"-f bkg model (MC file) "<<std::endl;
  std::cout<<"-F bias file"<<std::endl;
  std::cout<<"-L LY file"<<std::endl;
  std::cout<<"-P PCA file"<<std::endl;
  std::cout<<"-m is improved"<<std::endl;
  std::cout<<"-e float efficiency 1 = yes 0 = no "<<std::endl;
  std::cout<<"-r prior on the xiratio, arguments are the mean and uncertainity seperated by comma"<<std::endl;
  std::cout<<"-b float bias     see above"<<std::endl;
  std::cout<<"-l full efficiency parameter is the value of efficiency in the MC"<<std::endl;
  std::cout<<"-C constant eff values  - comma seperated list mean,sigma"<<std::endl;
  std::cout<<"-B fix background parameter is the value to fix it to"<<std::endl;
  std::cout<<"-T threshold"<<std::endl;
  std::cout<<"-i binning"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
}

int main(int argc,char **argv)
{


  bool floatbias=false;
  bool floateff=false;
  bool floatbkg=true;
  double fulleff=0;
  double effvalMC;
  bool is_imp=0;
  double bkg=1;
  double threshold =100;
  TString theory_file="../inputs/SSD_AE.root";
  TString bkg_file="../inputs/M1_data_bkg_model.root";
  TString bias_file="../inputs/output_marg.root";
  TString LY_file="../inputs/LY_output_marg.root";
  TString PCA_file="../inputs/PCA_output_marg.root";
  double constant_eff_mean= 0.9409;
  double constant_eff_sigma=0.0099;
  std::string name="Fits/TestFit";
  double effMCVal=0.891;
  int binning=1;
  bool prior_ratio=0;
  double xiratio=0;
  double xiratio_error=0;
  
  {

    static struct option long_options[] = {
					   { "name", required_argument, nullptr,'N'},
                                           {"theory-file",required_argument,nullptr,'S'},
					   { "is-improved",required_argument,nullptr,'m'},
					   {"bkg-file",required_argument,nullptr,'f'},
					   {"bias-file",required_argument,nullptr,'F'},
					   {"LY-file",required_argument,nullptr,'L'},
					   {"xi-ratio-prior",required_argument,nullptr,'r'},
					   {"PCA-file",required_argument,nullptr,'P'},
					   { "float-eff",required_argument,nullptr,'e'},
                                           { "float-bias",required_argument,nullptr,'b'},
					   { "full-eff",required_argument,nullptr,'l'},
					   { "constant-eff",required_argument,nullptr,'C'},
					   { "fix-background",required_argument,nullptr,'B'},
					   { "threshold",required_argument,nullptr,'T'},
					   { "binning",required_argument,nullptr,'i'},
					   { "help",                       no_argument, nullptr, 'h' },
					   {nullptr, 0, nullptr, 0}
    };

    const char* const short_options = "N:S:f:F:l:e:b:r:B:T:L:P:C:i:m:h";
    int c;

    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
      switch (c) {
      case 'N': {
	name = optarg;
        break;
      }
      case 'S': {
        theory_file = optarg;
        break;
      }
      case 'f': {
        bkg_file = optarg;
        break;
      }
      case 'F': {
        bias_file = optarg;
        break;
      }
      case 'm':{
	is_imp=atoi(optarg);
	break;
      }
      case 'l':{
	fulleff=1;
	effvalMC=std::stod(optarg);
	break;
      }
      case 'e': {
	floateff=std::stoi(optarg);
	break;
      }
      case 'b': {
	floatbias=std::stoi(optarg);
        break;
      }
      case 'C':{
	std::string input=optarg;
	std::istringstream ss(input);
	std::string token;

	std::getline(ss, token, ',');
	constant_eff_mean=std::stod(token);
	std::getline(ss, token, ',');		
	constant_eff_sigma=std::stod(token);
	break;
      }
      case 'r':{
        std::string input=optarg;
        std::istringstream ss(input);
        std::string token;
	prior_ratio=1;
        std::getline(ss, token, ',');
        xiratio=std::stod(token);
        std::getline(ss, token, ',');
        xiratio_error=std::stod(token);
        break;
      }

      case 'B': {
	floatbkg=false;
	bkg=std::stod(optarg);
        break;
      }
      case 'T':{
	threshold=std::stod(optarg);
	break;
      }
      case 'L':{
	LY_file=optarg;
	break;
      }
      case 'P':{
	PCA_file=optarg;
	break;
      }
      case 'i':{
	binning=std::stoi(optarg);
	break;
      }
      case 'h':{
        Usage();
        return 0;
	break;
      }

	
      default: {
	exit(1);
      }
      }
    }
  }


  std::cout<<"Running the fit with: "<<std::endl;
  std::cout<<"name                : "<<name<<std::endl;
  std::cout<<"floatbias           : "<<floatbias<<std::endl;
  std::cout<<"floateff            : "<<floateff<<std::endl;
  std::cout<<"floatbkg            : "<<floatbkg<<"    "<<bkg<<std::endl;
  std::cout<<"fulleff             : "<<fulleff<<std::endl;
  std::cout<<"eff MC              : "<<effvalMC<<std::endl;
  std::cout<<"threshold           : "<<threshold<<std::endl;
  std::cout<<"binning             : "<<binning<<std::endl;
  std::cout<<"theory file         : "<<theory_file<<std::endl;
  std::cout<<"bkg file            : "<<bkg_file<<std::endl;
  std::cout<<"bias file           : "<<bias_file<<std::endl;
  std::cout<<"LY file             : "<<LY_file<<std::endl;
  std::cout<<"PCA file            : "<<PCA_file<<std::endl;
  std::cout<<"const eff           : "<<constant_eff_mean<<" , "<<constant_eff_sigma<<std::endl;
  std::cout<<"prior ratio         : "<<xiratio<<" +/- "<<xiratio_error<<std::endl;
  std::cout<<" "<<std::endl;
  
 
  

  




  
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
    
    // create new QbbFit object
    QbbFit m(std::string("fits/")+name,theory_file,bkg_file,bias_file,!floateff,floatbias,fulleff,effvalMC,is_imp);
    m.SetConstantEff(constant_eff_mean,constant_eff_sigma);
    m.SetLYPCA(LY_file,PCA_file);
    m.SetThreshold(threshold);
    m.SetBinning(binning);
    if (prior_ratio)
      m.SetRatioPrior(xiratio,xiratio_error);
    
    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    
    BCLog::OutSummary("Test model created");
    
    //////////////////////////////
    // perform your analysis here


    std::vector<double> initialPos ;
    m.FindMode(initialPos);
    std::vector<double> bestFitMinuit = m.GetBestFitParameters();                                                                                                                                       
    
    m.SetInitialPositions(bestFitMinuit);  
    int i=0;
    for (auto & par:bestFitMinuit)
      {
	std::cout<<"par("<<i<<") = "<<par<<std::endl;
	i++;
      }

    // Write Markov Chain to a ROOT file as a TTree
    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit

    m.FindMode(m.GetBestFitParameters());
    TCanvas *cK =new TCanvas();

    TCanvas *c =new TCanvas();
    m.MakePlots(c,cK,m.GetBestFitParameters());
    c->Print("fits/"+(TString)m.GetSafeName()+TString("_fit.C"));
    cK->Print("fits/"+(TString)m.GetSafeName()+TString("_kurie_fit.C"));

	
    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized("fits/"+m.GetSafeName() + "_plots.pdf");

    // print summary plots
     m.PrintParameterPlot("fits/"+m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
     m.PrintCorrelationMatrix("fits/"+m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots("fits/"+m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
