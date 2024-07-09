#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>


void ResetErrors(TH1D*&h,double N=1e8)
{
  h->Scale(N/h->Integral());
  for (int i=1;i<h->GetNbinsX();i++)
    {
      h->SetBinError(i,0,sqrt(h->GetBinContent(i)));
    }
}
     

void interpolate_spectrum(bool imp,std::string path,std::string plot,std::string output,std::map<std::string,double>Gs,double Qbb,int order

)
{
  // the aim of this macro is to first compare some plots of 2vbb shape shifting Qbb

  TFile * f = new TFile(path.c_str());
  TH1D *h2vbb;
  TH1D*h0;
  TH1D*h2;
  TH1D*h22;
  TH1D*h4;

  std::vector<double>G;
  std::vector<TH1D*>contributions;
  std::vector<TString> names;
  double N0;
  double N;
  if (!imp)
    {
      if (!(f->FindKey("Crystal_2n2b_SSD"))) 
      {
        throw std::runtime_error("Error you requested standard mode but do not have the neccesary objects in the file");
      }
      h2vbb= (TH1D*)f->Get("Crystal_2n2b_SSD");
      contributions.push_back(h2vbb);
      names.push_back("2vbb");
      N0=h2vbb->Integral();
    }
  else
    {
      
      if (!(f->FindKey("Crystal_2n2b_0sfs") || 
          f->FindKey("Crystal_2n2b_2sfs") || 
          f->FindKey("Crystal_2n2b_22sfs")||
          f->FindKey("Crystal_2n2b_4sfs"))
      )
      {
        f->ls();
        throw std::runtime_error("Error you requested improved model mode but do not have the neccesary objects in the file");

      }

      h0=(TH1D*)f->Get("Crystal_2n2b_0sfs");
      h2=(TH1D*)f->Get("Crystal_2n2b_2sfs");
      h22=(TH1D*)f->Get("Crystal_2n2b_22sfs");
      h4=(TH1D*)f->Get("Crystal_2n2b_4sfs");
      contributions.push_back(h0);
      contributions.push_back(h2);
      contributions.push_back(h22);
      contributions.push_back(h4);
      names.push_back("h0");
      names.push_back("h2");
      names.push_back("h22");
      names.push_back("h4");
      G.push_back(Gs["G0"]);
      G.push_back(Gs["G2"]);
      G.push_back(Gs["G22"]);
      G.push_back(Gs["G4"]);
      N = h0->Integral()+h2->Integral()+h22->Integral()+h4->Integral();
      N0=h0->Integral()/Gs["G0"];
     
      ResetErrors(h0);
      ResetErrors(h2);
      ResetErrors(h22);
      ResetErrors(h4);
    }
 
  double xi31=0.368;
  double xi51=0.4;
  TFile *fout;

  if (!imp)
    fout= new TFile(output.c_str(),"recreate");
  else
    fout=new TFile(output.c_str(),"recreate");
  
		   
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,1);
  TH1D *htmp;
  c1->Print((plot+"(").c_str(),"pdf");

  for (int i=0;i<contributions.size();i++)
    {
      htmp=contributions[i];
      
      
      
      TH1D *hA = (TH1D*)htmp->Clone(Form("hA_%s",names[i].Data()));
      
      // Create the efficiency function
      TF1 *fe = new TF1("fe",Form("1+x*[0]/%f",Qbb));
      fe->SetParameter(0,0);

      // fill a histogram of the A(E)
      for (int i=1;i<htmp->GetNbinsX()+1;i++)
      {
        hA->SetBinContent(i,fe->Eval(htmp->GetBinCenter(i))*htmp->GetBinContent(i)*pow(Qbb-htmp->GetBinCenter(i),-5));
        hA->SetBinError(i,0,fe->Eval(htmp->GetBinCenter(i))*htmp->GetBinError(i)*pow(Qbb-htmp->GetBinCenter(i),-5));
        
      }
      c1->cd(1);
      hA->SetLineColor(1);
      hA->GetXaxis()->SetRangeUser(0,2950);
      hA->SetTitle(Form("Predicted A(E) for %s",names[i].Data()));
      hA->Draw("");
      fout->cd();

      // creatae the nth order polynomial of A(E) and fit it (initial guesses)
      TF1 *fu = new TF1(Form("fu_%s",names[i].Data()),Form("pol%i",order));
      hA->Fit(fu,"REM","",0,2950);
      
      // scale the polynomial by (Q-E)^5
      TF1 *f2vbb = new TF1("f2vbb",Form("(x<[%i])*fu_%s(x)*pow(([%i]-x),5)*[%i]",order+1,names[i].Data(),order+1,order+2),0,4000);

      f2vbb->SetParameter(order+1,Qbb);
      f2vbb->FixParameter(order+1,Qbb);
      f2vbb->FixParameter(order+2,1);

  
      c1->cd(2);
      htmp->Draw();
      htmp->SetLineColor(1);

      // fit this with a likelihood fit
      htmp->Fit(f2vbb,"REML","",50,3250);
      gStyle->SetOptFit(1);
      std::cout<<f2vbb->GetChisquare()<<std::endl;
      std::cout<<f2vbb->GetNDF()<<std::endl;
      
      f2vbb->SetParameter(order+2,N0*G[i]/1e8);
      TF1 *f2= new TF1("f2",Form("f2vbb(x)*pow((%f-x),-5)",Qbb),0,4000);

      // Write just the A(E) part
      f2->Write(Form("fu_%s",names[i].Data()));
      c1->Print(plot.c_str(),"pdf");
    }
  c1->Print((plot+")").c_str(),"pdf");
  fout->Close();
}

void usage(std::string arg)
{
    std::cout<<"The usage of this function is:"<<std::endl;
    std::cout << arg << " -l [plot_path] -i [input_path] -p [polynomial order] -q [Qbb] -o [outpath] " << std::endl;
    std::cout<< "add the option -I to use the improved model and -h for this help"<<std::endl;
}
int main(int argc, char* argv[])
{
  
  std::string plot_path = "plots/beta_debug.pdf";
  bool is_improved_model = false;
  std::string input_path = "inputs/M1_data_bkg_model_2n2b_improved.root";
  int poly_order = 6;
  double Qbb = 3034.4;
  std::map<std::string,double> G;
  G["G0"]= 3.279;
  G["G2"]=1.498;
  G["G22"]=0.1972;
  G["G4"]=0.8576;
  std::string output_path= "inputs/AofE.root";

  int opt;
  while ((opt = getopt(argc, argv, "l:Ii:p:q:o:h")) != -1) 
    {
      switch (opt) 
      {
            case 'l':
            {
              plot_path  = optarg;
              break;
            }
            case 'I':
            {
              is_improved_model=true;
              break;
            }
            case 'i':
            {
              input_path=optarg;
              break;
            }
            case 'p':
            {
              poly_order=atoi(optarg);
              break;
            }
            case 'q':
            {
              Qbb = atof(optarg);
              break;
            }
            case 'o':
            {
              output_path = optarg;
              break;
            }
            case 'h':
            {
              usage(argv[0]);
              return 1;
            }
            default:
            {
              usage(argv[0]);
              return 1;
            }
  
      }
    }



  interpolate_spectrum(is_improved_model,input_path,plot_path,output_path,G,Qbb,poly_order);
 
}
