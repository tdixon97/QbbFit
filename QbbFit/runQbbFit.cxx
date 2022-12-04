// ***************************************************************
// This file was created using the bat-project script
// for project QbbFit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>

#include "QbbFit.h"

int main()
{

  // parse the options




  
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
    
    // create new QbbFit object
    bool fix=1;
    QbbFit m("TestFit","../SSD_AE.root","../M1_data_bkg_model.root",fix);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Test model created");

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    m.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit

    m.FindMode(m.GetBestFitParameters());
    TCanvas *cK =new TCanvas();

    TCanvas *c =new TCanvas();
    m.MakePlots(c,cK,m.GetBestFitParameters());
    c->Print((TString)m.GetSafeName()+TString("_fit.C"));
    cK->Print((TString)m.GetSafeName()+TString("_kurie_fit.C"));

	
    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary plots
     m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
     m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
