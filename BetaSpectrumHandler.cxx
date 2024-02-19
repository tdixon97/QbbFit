void ResetErrors(TH1D*&h,double N=1e8)
{
  h->Scale(N/h->Integral());
  for (int i=1;i<h->GetNbinsX();i++)
    {
      h->SetBinError(i,0,sqrt(h->GetBinContent(i)));
    }
}
     

void macro(double &Q,double &eQ,double effslope=0,bool imp=1)
{
  // the aim of this macro is to first compare some plots of 2vbb shape shifting Qbb


  double G0= 3.279;
  double G2=1.498;
  double G22=0.1972;
  double G4=0.8576;
  TFile * f;

  if (imp==0)
    f = new TFile("inputs/M1_data_bkg_model_2n2b_fix.root");
  else
    f=new TFile("inputs/M1_data_bkg_model_2n2b_improved.root");

  

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
      h2vbb= (TH1D*)f->Get("Crystal_2n2b_SSD");
      contributions.push_back(h2vbb);
      names.push_back("2vbb");
      N0=h2vbb->Integral();
    }
  else
    {
      f->Print();
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
      G.push_back(G0);
      G.push_back(G2);
      G.push_back(G22);
      G.push_back(G4);
      N = h0->Integral()+h2->Integral()+h22->Integral()+h4->Integral();
      N0=h0->Integral()/G0;
     
      ResetErrors(h0);
      ResetErrors(h2);
      ResetErrors(h22);
      ResetErrors(h4);
    }
  std::cout<<"N "<<N<<std::endl;
  std::cout<<"N0 = "<<N0<<std::endl;
  std::cout<<"N0*G0 = "<<N0*G0<<std::endl;
  double xi31=0.368;
  double xi51=0.4;
  double NT=N0*G0+N0*G2*xi31+N0*G22*xi31*xi31/3.+N0*G4*(xi31*xi31/3.+xi51);
  std::cout<<"NT = "<<NT<<std::endl;
  TFile *fout;

  if (!imp)
    fout= new TFile("inputs/SSD_AE.root","recreate");
  else
    fout=new TFile("inputs/improved_AE.root","recreate");
  
		   
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,1);
  TH1D *htmp;
  c1->Print("beta_debug.pdf(","pdf");
  for (int i=0;i<contributions.size();i++)
    {
      htmp=contributions[i];
      
      // c1->Divide(3,1);
      c1->cd(1);
      htmp->Draw("HIST");
      
      TH1D *hA = (TH1D*)htmp->Clone(Form("hA_%s",names[i].Data()));
      
      TF1 *fe = new TF1("fe","1+x*[0]/3034.");
      fe->SetParameter(0,effslope);
      for (int i=1;i<htmp->GetNbinsX()+1;i++)
	{
	  hA->SetBinContent(i,fe->Eval(htmp->GetBinCenter(i))*htmp->GetBinContent(i)*pow(3034.4-htmp->GetBinCenter(i),-5));
	  hA->SetBinError(i,0,fe->Eval(htmp->GetBinCenter(i))*htmp->GetBinError(i)*pow(3034.4-htmp->GetBinCenter(i),-5));
	  
	}
      c1->cd(2);
      hA->SetLineColor(1);
      hA->GetXaxis()->SetRangeUser(0,2950);
      hA->SetTitle(Form("Predicted A(E) for %s",names[i].Data()));
      hA->Draw("");
      fout->cd();
      TF1 *fu = new TF1(Form("fu_%s",names[i].Data()),"pol6");
      hA->Fit(fu,"REM","",0,2950);
      
      TF1 *f2vbb = new TF1("f2vbb",Form("(x<[7])*fu_%s(x)*pow(([7]-x),5)*[8]",names[i].Data()),0,4000);

      f2vbb->SetParameter(7,3034.4);
      f2vbb->FixParameter(7,3034.4);
      f2vbb->FixParameter(8,1);

      
      c1->cd(3);
      htmp->Draw();
  
      htmp->SetLineColor(1);
  
      htmp->Fit(f2vbb,"REM","",50,3250);
      gStyle->SetOptFit(1);
      std::cout<<f2vbb->GetChisquare()<<std::endl;
      std::cout<<f2vbb->GetNDF()<<std::endl;
      f2vbb->Draw("same");
      f2vbb->SetParameter(8,N0*G[i]/1e8);
      TF1 *f2= new TF1("f2","f2vbb(x)*pow((3034.4-x),-5)",0,4000);

      f2->Write(Form("fu_%s",names[i].Data()));

      c1->Print("beta_debug.pdf","pdf");
    }
  c1->Print("beta_debug.pdf)","pdf");
  fout->Close();
}


/*
      // now compare shifting the Q-value
      TCanvas *c2= new TCanvas("c2","c2");
  c2->cd();
  std::vector<double>shifts{100,20,10,0,-10,-20,-100};
  int counter=0;
  h2vbb->Draw("HIST");
  //f2vbb->Draw("Csame");
  double norm=f2vbb->Integral(0,4000);
  TLegend * l = new TLegend(0.7,0.7,0.9,0.9);
  //l->AddEntry(h2vbb,"Shift of 0 keV ");
  for (auto &shift:shifts)
    {
      TF1 *fi = new TF1(Form("f2vbb_%i",(int)shift),"(x<[7])*fu(x)*pow(([7]-x),5)*[8]",0,4000);
           fi->SetNpx(1000);
      fi->SetParameter(7,3034.4-shift);
      //      fi->SetParameter(8,1);
      fi->SetLineColor(counter+1);
      fi->SetParameter(8,1);
      l->AddEntry(fi,Form("Shift of %i keV",(int)shift));
      fi->SetParameter(8,norm/(double)fi->Integral(0,10000));
      fi->SetLineWidth(1);
      fi->Draw("Csame");

      counter++;
      
    }

  l->Draw();

  TCanvas *c3= new TCanvas("c3","c3");

  // try a simple Kurie plot

  TH1D * hdata =(TH1D*)f->Get("Experimental_data");
  TH1D * hBkg =(TH1D*)f->Get("Background_model");

  TH1D *hModel=(TH1D*)hBkg->Clone();
  hModel->Add(h2vbb);


  // draw
  TCanvas *cBkg = new TCanvas("cBkg","cBkg");
  cBkg->cd();
  // hdata->Draw();
  //hBkg->SetLineWidth(0);
  TLegend *lb= new TLegend(0.7,0.7,0.9,0.9);
  lb->AddEntry(hdata,"Data");
  lb->AddEntry(hModel,"Model");
  lb->AddEntry(hBkg,"Bkg");
  lb->AddEntry(h2vbb,"Signal");

  hModel->SetLineWidth(1);
  hdata->SetLineWidth(1);
  hBkg->SetLineWidth(1);
  h2vbb->SetLineWidth(1);
  hModel->SetLineColor(2);
  hdata->Draw("E");
  //h2vbb->SetLineWidth(0);
  h2vbb->SetLineColor(9);
  h2vbb->Draw("HISTsame");
  hBkg->Draw("HISTsame");
  hModel->Draw("HISTsame");
  lb->Draw();
  TCanvas *cBkg2=new TCanvas("cBkg2","cBkg2",700,200);
  // cBkg->Draw();
  cBkg2->cd();
  TH1D *hRatio = (TH1D*)h2vbb->Clone("hRatio");
  hRatio->Divide(hBkg);
  hRatio->Draw();
  hRatio->SetTitle(" ; Energy [keV]; Signal/Background ratio ; ");

  //hdata->Rebin((int)binning);
  // hBkg->Rebin((int)binning);
  //hdata->Draw();
  //hBkg->Draw("same");
  TH1D *hsub = (TH1D*)hdata->Clone("hsub");
  for (int i=1;i<hBkg->GetNbinsX()+1;i++)
    {
      hsub->SetBinContent(i,hsub->GetBinContent(i)-factor*hBkg->GetBinContent(i));
      hsub->SetBinError(i,0,hsub->GetBinError(i));
    }
  
  //  hsub->Rebin((int)binning);
  TCanvas* c4 = new TCanvas("c4","c4");
  hsub->SetTitle("Bkg subtracted spectrum");
  hsub->Draw();
  f2vbb->Draw("same");


  // now built the Kurie plot
  TCanvas *c5 = new TCanvas("c5","c5");


  TH1D *hKurie5 = (TH1D*)hsub->Clone("hKurie5");
  TH1D *hKurie = (TH1D*)hsub->Clone("hKurie");
  hKurie5->SetTitle("Plot of Data-Bkg/A(E) ; Energy ; Data-Bkg/A(E); ");
  for (int i=1;i<hKurie5->GetNbinsX()+1;i++)
    {
      fu->SetRange(0,4000);
      hKurie5->SetBinContent(i,hKurie5->GetBinContent(i)/(binning*fu->Eval(hKurie5->GetBinCenter(i))));
      hKurie5->SetBinError(i,0,hKurie5->GetBinError(i)/(binning*fu->Eval(hKurie5->GetBinCenter(i))));

      if (hKurie->GetBinContent(i)>0)
	hKurie->SetBinContent(i,pow(hKurie->GetBinContent(i)/(binning*fu->Eval(hKurie->GetBinCenter(i))),0.2));
      else if (hKurie->GetBinContent(i)<0)
	hKurie->SetBinContent(i,-pow(-hKurie->GetBinContent(i)/(binning*fu->Eval(hKurie->GetBinCenter(i))),0.2));
      else
	hKurie->SetBinContent(i,0);
      if (hdata->GetBinContent(i)==0)
	hKurie->SetBinContent(i,0);
      
      if (hKurie->GetBinContent(i)>0)
	hKurie->SetBinError(i,0,0.2*pow((hsub->GetBinContent(i)/(binning*fu->Eval(hKurie->GetBinCenter(i)))),-0.8)*hKurie->GetBinError(i)/(binning*fu->Eval(hKurie->GetBinCenter(i))));
      else if (hKurie->GetBinContent(i)<0)
	hKurie->SetBinError(i,0,0.2*pow(-(hsub->GetBinContent(i)/(binning*fu->Eval(hKurie->GetBinCenter(i)))),-0.8)*hKurie->GetBinError(i)/(binning*fu->Eval(hKurie->GetBinCenter(i))));
      else
	hKurie->SetBinError(i,0,0);

      if (hdata->GetBinContent(i)==0)
        hKurie->SetBinError(i,i,0);

    }

  TF1 *fpol5 = new TF1("fpol1","[0]*pow([1]-x,5)",fitlow,3200);
  fpol5->SetParameter(0,1);
  fpol5->SetParameter(1,3034.4);
  hKurie5->Fit(fpol5,"REM","",fitlow,3200);
  hKurie5->Draw();  TLatex *tlat = new TLatex();
  tlat->DrawLatex(2000,1e15,Form("Q_{#beta#beta} = %0.3f +/- %0.4f keV",fpol5->GetParameters()[1],fpol5->GetParErrors()[1]));


  gStyle->SetOptFit(1);
  
  TCanvas *c6= new TCanvas("c6","c6");
  gStyle->SetOptFit(1);
  hKurie->SetTitle("Kurie Plot for 2vbb of 100Mo; Energy [keV] ; [(Data-Background)/(A(E))]^{1/5} ; ");
  hKurie->Draw();
  TF1 *fpol1 = new TF1("fpol1","[0]*([1]-x)",fitlow,3200);
  fpol1->SetParameter(0,-1);
  fpol1->SetParameter(1,3034);
  fpol1->SetParNames("Offset","Q_{#beta#beta}");
  hKurie->Fit(fpol1,"REM","",fitlow,3200);
  gStyle->SetOptFit(1);
  tlat->DrawLatex(2000,3000,Form("Q_{#beta#beta} = %0.3f +/- %0.4f keV",fpol1->GetParameters()[1],fpol1->GetParErrors()[1]));


  // systematics


  Q=fpol1->GetParameters()[1];
  eQ =fpol1->GetParErrors()[1];

}

*/
void Loop()
{
  
  double Q;double eQ;
  macro(Q,eQ,0,1);
 
}
