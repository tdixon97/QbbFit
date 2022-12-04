void macro(double &Q,double &eQ,double binning=20,double fitlow=300,double factor=1,double effslope=0)
{
  // the aim of this macro is to first compare some plots of 2vbb shape shifting Qbb


  TFile * f= new TFile("M1_data_bkg_model.root");


  TH1D *h2vbb = (TH1D*)f->Get("Crystal_2n2b_SSD");
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,1);
  c1->cd(1);
  h2vbb->Draw("HIST");


  // now lets compute the estimate A(E)

  TH1D *hA = (TH1D*)h2vbb->Clone("hA");

  TF1 *fe = new TF1("fe","1+x*[0]/3034.");
  fe->SetParameter(0,effslope);
  for (int i=1;i<h2vbb->GetNbinsX()+1;i++)
    {
      hA->SetBinContent(i,fe->Eval(h2vbb->GetBinCenter(i))*h2vbb->GetBinContent(i)*pow(3034.4-h2vbb->GetBinCenter(i),-5));
      hA->SetBinError(i,0,fe->Eval(h2vbb->GetBinCenter(i))*h2vbb->GetBinError(i)*pow(3034.4-h2vbb->GetBinCenter(i),-5));
	    
    }
  c1->cd(2);
  hA->SetLineColor(1);
  hA->GetXaxis()->SetRangeUser(0,2950);
  hA->SetTitle("Predicted A(E)");
  hA->Draw("");
  TFile *fout = new TFile("SSD_AE.root","recreate");
  TF1 *fu = new TF1("fu","pol6");
  hA->Fit(fu,"REM","",0,2950);
  fu->Write();
  // now reconstruct the shape
  TF1 *f2vbb = new TF1("f2vbb","(x<[7])*fu(x)*pow(([7]-x),5)",0,4000);

  f2vbb->SetParameter(7,3034.4);
  f2vbb->FixParameter(7,3034.4);
  c1->cd(3);
  h2vbb->Draw();
  
  h2vbb->SetLineColor(1);
  
  h2vbb->Fit(f2vbb,"REM","",50,3250);
  gStyle->SetOptFit(1);
  std::cout<<f2vbb->GetChisquare()<<std::endl;
  std::cout<<f2vbb->GetNDF()<<std::endl;
  f2vbb->Draw("same");

  TF1 *f2= new TF1("f2","f2vbb(x)*pow((3034.4-x),-5)",0,4000);

  f2->Write("fu");

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

  // create a CUPID-output
  /*
  TFile *fout_CUPID=new TFile("M1_data_bkg_model_CUPID.root","recreate");
  TH1D *hBkg_CUPID=(TH1D*)hBkg->Clone("hBkg_CUPID");
  TH1D *h2vbb_CUPID = (TH1D*)h2vbb->Clone("Crystal_2n2b_SSD_CUPID");

  hBkg_CUPID->Scale(5000/(2.71*30.));
  h2vbb_CUPID->Scale(5000/2.71);

  int NCUPID=std::round((h2vbb_CUPID->Integral()+hBkg_CUPID->Integral()));
  TH1D *hdata_CUPID=(TH1D*)hdata->Clone("hdata_CUPID");
  TH1D *hModel_CUPID=(TH1D*)hBkg_CUPID->Clone("hModel_CUPID");
  hModel_CUPID->Add(h2vbb_CUPID);
  hdata_CUPID->Clear();
  hdata_CUPID->Reset();
  std::cout<<NCUPID<<std::endl;
  for (int i=0;i<NCUPID;i++)
    {
      double E=hModel_CUPID->GetRandom();
      hdata_CUPID->Fill(E);
    }
  hBkg_CUPID->Write("Background_model");
  h2vbb_CUPID->Write("Crystal_2n2b_SSD");
  hdata_CUPID->Write("Experimental_data");
  fout->cd();
  */

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


void Loop()
{
  
  double Q;double eQ;
  macro(Q,eQ,20,300,1);
  /*
  std::cout<<"Standard result is "<<Q<<" +/- "<<eQ<<" keV"<<std::endl;
  double Qdef=Q; 
  TGraphErrors *gBkg = new TGraphErrors();

  std::vector<double>factor{0.95,0.96,0.97,0.98,0.99,1.0,1.01,1.02,1.03,1.04,1.05};
  int counter=0;
  
  TCanvas *csyst=new TCanvas("csyst","csyst");
  csyst->cd();
  csyst->Divide(2,2);
  csyst->cd(1);
  for (auto &f:factor)
    {
      macro(Q,eQ,20,300,f);
      gBkg->SetPoint(counter,f,(Q-Qdef));
      gBkg->SetPointError(counter,0,eQ);
      counter++;
      
    }
  csyst->cd(1);
  gBkg->SetTitle("Bkg subtraction systramtic ; Bkg Subtraction factor ; Shift [keV]; ");
  gBkg->Draw("APE");

  std::vector<double>low{200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700};
  counter=0;
  csyst->cd(2);
  TGraphErrors *gRange = new TGraphErrors();
  for (auto &l:low)
    {
      macro(Q,eQ,20,l,1);
      gRange->SetPoint(counter,l,Q-Qdef);
      gRange->SetPointError(counter,0,eQ);
      counter++;
    }
  csyst->cd(2);
  gRange->SetTitle("Fit range systematic ; Low range [keV] ; Shift [keV] ; ");
  gRange->Draw("APE");
  //  csyst->Update();

  std::vector<double>binning{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35,40,45,50,100,200};
  TGraphErrors *gBinning = new TGraphErrors();
  counter=0;
  for (auto &b:binning)
    {
      macro(Q,eQ,b,300,1);
      gBinning->SetPoint(counter,b,Q-Qdef);
      gBinning->SetPointError(counter,0,eQ);
      counter++;
    }
  csyst->cd(3);
  gBinning->SetTitle("Binning range systematic ; Low range [keV] ; Shift [keV] ; ");
  gBinning->Draw("APE");


  // another possiblitiy is some slope to efficiency

  std::ifstream ifile;
  ifile.open("Eff.cfg");
  double e1,e2,p1,p2,ep1,ep2;
  ifile>>e1; ifile>>p1; ifile>>ep1;
  ifile>>e2; ifile>>p2; ifile>>ep2;

  p1/=e1; ep1/=e1;
  p2/=e2; ep2/=e2;

  TCanvas *ceff= new TCanvas("ceff","ceff");
  TF1 *feff= new TF1("feff","[0]+[1]*x",0,3034);
  feff->SetParameter(0,1);
  feff->SetParameter(1,p1+p2);
  std::cout<<"p1+p2 = "<<p1+p2<<" +/- "<<sqrt(ep1*ep1+ep2*ep2)<<std::endl;
  feff->Draw();
  
  std::vector<double>effslope{-10,-5,-2,-1,0,1,2,5,10};
  TGraphErrors *gEff = new TGraphErrors();
  counter=0;
  for (auto &e:effslope)
    {
      macro(Q,eQ,20,300,1,e/100.);
      gEff->SetPoint(counter,e,Q-Qdef);
      gEff->SetPointError(counter,0,eQ);
      counter++;
    }
  csyst->cd(4);
  gEff->SetTitle("Efficiency slope stematic ; Slope [%] ; Shift [keV] ; ");
  gEff->Draw("APE");
  */
}
