#include <stdio.h>

#include "Gaus.h"
#include "Function.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TLine.h"

int main(){
  Function prova;
  const Int_t nfunc = 2;
  prova.SetNtype(nfunc);
  
  TF1 *fAbundances = new TF1("fAbundances","0.8 - 0.3*x",0,1);
  TF1 *fTrend[2];
  fTrend[0] = new TF1("fTrend0","x",0,1);
  fTrend[1] = new TF1("fTrend1","1-x",0,1);


  TF1 *ff1 = new TF1("ff1","gaus",-20,20);
  ff1->SetParameter(0,1);
  ff1->SetParameter(1,-1);
  ff1->SetParameter(2, 1);
  ff1->SetParameter(0,1./ff1->Integral(-20,20));
  TF1 *ff2 = new TF1("ff2","gaus",-20,20);
  ff2->SetParameter(0,1);
  ff2->SetParameter(1,1);
  ff2->SetParameter(2, 1);
  ff2->SetParameter(0,1./ff2->Integral(-20,20));

  prova.SetResponseFunction(0,ff1);
  prova.SetResponseFunction(1,ff2);

  TF1 *f1 = prova.GetProbabilityDensity(0);
  TF1 *f2 = prova.GetProbabilityDensity(1);

  prova.SetMatrix();

  ff1 = prova.GetProbabilityDensityStar(0);
  ff2 = prova.GetProbabilityDensityStar(1);

  f1->SetLineColor(4);
  f2->SetLineColor(2);

  ff1->SetLineColor(4);
  ff2->SetLineColor(2);

  ff1->SetLineStyle(2);
  ff2->SetLineStyle(2);

  for(Int_t i=0;i<nfunc;i++)
    for(Int_t j=0;j<nfunc;j++)
      printf("Scalar product <%i|%i> = %f\n",i,j,prova.ScalarProduct(i,j));

  TCanvas *c = new TCanvas();
  TH1F *hback = new TH1F("hback",";S;df/dS",1000,-6,6);
  hback->Draw();
  hback->SetStats(0);
  hback->SetMaximum(1.6);
  hback->SetMinimum(-0.6);

  ff1->Draw("SAME");
  ff2->Draw("SAME");
  f1->Draw("SAME");
  f2->Draw("SAME");

  TLegend *leg = new TLegend(0.15,0.5,0.3,0.85);
  leg->SetHeader("f(S) =");
  leg->AddEntry(f1,"p_{1}(S)","l");
  leg->AddEntry(f2,"p_{2}(S)","l");
  leg->AddEntry(ff1,"#psi_{1}(S)","l");
  leg->AddEntry(ff2,"#psi_{2}(S)","l");

  leg->Draw("SAME");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);

  c->Update();

  hback->GetXaxis()->SetTitleSize(0.05);
  hback->GetYaxis()->SetTitleSize(0.05);
  hback->GetXaxis()->SetNdivisions(408);
  hback->GetYaxis()->SetNdivisions(408);
  hback->GetXaxis()->SetTitleOffset(0.9);
  hback->GetYaxis()->SetTitleOffset(0.9);

  c->Update();

  Float_t val,signal;
  Int_t isp;

  TF1 *rf[nfunc] = {f1,f2};

  TH1F *hsig = new TH1F("hsig",";signal",600,-6,6);
  TH1F *hsig1 = new TH1F("hsig1",";signal",600,-6,6);
  TH1F *hsig2 = new TH1F("hsig2",";signal",600,-6,6);

  TH1F *hsim1 = new TH1F("hsim1",";signal",600,-6,6);
  TH1F *hsim2 = new TH1F("hsim2",";signal",600,-6,6);

  hsig1->SetLineColor(4);
  hsig2->SetLineColor(2);

  hsim1->SetLineColor(4);
  hsim2->SetLineColor(2);

  hsim1->SetLineStyle(2);
  hsim2->SetLineStyle(2);

  TProfile *hampl[nfunc];
  hampl[0] = new TProfile("hampl0",";signal;amplitude",1000,-20,20);
  hampl[1] = new TProfile("hampl1",";signal;amplitude",1000,-20,20);

  TProfile *hprob[nfunc];
  hprob[0] = new TProfile("hprob0",";signal;amplitude",65,-3,3.5);
  hprob[1] = new TProfile("hprob1",";signal;amplitude",65,-3,3.5);

  TProfile *htrend[2];
  htrend[0] = new TProfile("htrend1",";x;y",100,0,1);
  htrend[1] = new TProfile("htrend2",";x;y",100,0,1);

  TH1F *hx[2];
  hx[0] = new TH1F("hx1",";x",100,0,1);
  hx[1] = new TH1F("hx2",";x",100,0,1);

  TH2F *hxy = new TH2F("hxy",";x;y",120,-0.1,1.1,120,-0.1,1.1);
  TH2F *hxyA[2];
  hxyA[0] = new TH2F("hxy1",";x;y",120,-0.1,1.1,120,-0.1,1.1);
  hxyA[1] = new TH2F("hxy2",";x;y",120,-0.1,1.1,120,-0.1,1.1);
  TH2F *hxyS[2];
  hxyS[0] = new TH2F("hxyS1",";x;y",120,-0.1,1.1,120,-0.1,1.1);
  hxyS[1] = new TH2F("hxyS2",";x;y",120,-0.1,1.1,120,-0.1,1.1);

  TH1F *hpriors[2][20];
  TProfile *htrendPriors[2][20];
  TH2F *hxyPriors[2][20];
  for(Int_t i=0;i < 20;i++){
    hpriors[0][i] = new TH1F(Form("hprior1_%i",i),";x",100,0,1);
    hpriors[1][i] = new TH1F(Form("hprior2_%i",i),";x",100,0,1);

    htrendPriors[0][i] = new TProfile(Form("htrendPriors1_%i",i),";x;y",100,0,1);
    htrendPriors[1][i] = new TProfile(Form("htrendPriors2_%i",i),";x;y",100,0,1);

    hxyPriors[0][i] = new TH2F(Form("hxyPriors1_%i",i),";x;y",120,-0.1,1.1,120,-0.1,1.1);
    hxyPriors[1][i] = new TH2F(Form("hxyPriors2_%i",i),";x;y",120,-0.1,1.1,120,-0.1,1.1);

    if(i==0)
      for(Int_t j=1;j <= 100;j++){
	hpriors[0][i]->SetBinContent(j,0.5);
	hpriors[1][i]->SetBinContent(j,0.5);
      }
  }


  Int_t counts[nfunc] = {0,0};
  Float_t recocounts[nfunc] = {0,0};

  Float_t recocountsBayes[nfunc][20];


  Double_t ampl[nfunc],norm;

  Float_t x,y;

  for(Int_t i=0;i < 1000000;i++){
    val = gRandom->Rndm();
    x = gRandom->Rndm();

    if(val < fAbundances->Eval(x)) isp=0;
    else isp=1;

    y = fTrend[isp]->Eval(x) + gRandom->Gaus(0,0.05);

    hxy->Fill(x,y);

    counts[isp]++;

    signal = rf[isp]->GetRandom();

    ampl[0] = prova.EvaluateProb(0,signal);
    ampl[1] = prova.EvaluateProb(1,signal);

    hxyA[0]->Fill(x,y,ampl[0]);
    hxyA[1]->Fill(x,y,ampl[1]);
    hxyS[isp]->Fill(x,y);

    hampl[0]->Fill(signal,ampl[0]);
    hampl[1]->Fill(signal,ampl[1]);

    recocounts[0] += ampl[0];
    recocounts[1] += ampl[1];

    hx[0]->Fill(x,ampl[0]);
    hx[1]->Fill(x,ampl[1]);

    htrend[0]->Fill(x,y,ampl[0]);
    htrend[1]->Fill(x,y,ampl[1]);
  }

  printf("counts[0] = %i -- recocounts[0] = %f -- error[%c] = %f\n",counts[0],recocounts[0],'%',(recocounts[0]/counts[0] -1)*100);
  printf("counts[1] = %i -- recocounts[1] = %f -- error[%c] = %f\n",counts[1],recocounts[1],'%',(recocounts[1]/counts[1] -1)*100);

  
  recocountsBayes[0][0]=0.5;
  recocountsBayes[1][0]=0.5;
  Float_t istep[20];
  istep[0]=1;
  for(Int_t j=1;j < 20;j++){
    istep[j] = j+1;

    recocountsBayes[0][j] =0;
    recocountsBayes[1][j] =0;

    for(Int_t i=0;i < 1000000;i++){
      val = gRandom->Rndm();
      x = gRandom->Rndm();

      if(val < fAbundances->Eval(x)) isp=0;
      else isp=1;

      y = fTrend[isp]->Eval(x) + gRandom->Gaus(0,0.05);

      signal = rf[isp]->GetRandom();
      
      if(j==19){
	hsig->Fill(signal);
	if(isp==0) hsim1->Fill(signal);
	if(isp==1) hsim2->Fill(signal);
      }
//       ampl[0] = rf[0]->Eval(signal)*recocountsBayes[0][j-1];
//       ampl[1] = rf[1]->Eval(signal)*recocountsBayes[1][j-1];

      ampl[0] = rf[0]->Eval(signal)*hpriors[0][j-1]->Interpolate(x);
      ampl[1] = rf[1]->Eval(signal)*hpriors[1][j-1]->Interpolate(x);
      
      norm = ampl[0]+ampl[1] + 1E-10;
      
      ampl[0] /= norm;
      ampl[1] /= norm;
      
      if(j==19) hsig1->Fill(signal,ampl[0]);
      if(j==19) hsig2->Fill(signal,ampl[1]);

      if(j==19) hprob[0]->Fill(signal,ampl[0]);
      if(j==19) hprob[1]->Fill(signal,ampl[1]);

      recocountsBayes[0][j] += ampl[0];
      recocountsBayes[1][j] += ampl[1];

      hpriors[0][j]->Fill(x,ampl[0]);
      hpriors[1][j]->Fill(x,ampl[1]);

      htrendPriors[0][j]->Fill(x,y,ampl[0]);
      htrendPriors[1][j]->Fill(x,y,ampl[1]);

      hxyPriors[0][j]->Fill(x,y,ampl[0]);
      hxyPriors[1][j]->Fill(x,y,ampl[1]);
    }
    Float_t tot = recocountsBayes[0][j] + recocountsBayes[1][j];
    recocountsBayes[0][j] /= tot;
    recocountsBayes[1][j] /= tot;
  }

  TGraph *g[2];
  g[0] = new TGraph(20,istep,recocountsBayes[0]);
  g[1] = new TGraph(20,istep,recocountsBayes[1]);

  g[0]->SetLineColor(4);
  g[0]->SetMarkerColor(4);
  g[0]->SetMarkerStyle(20);
  g[1]->SetLineColor(2);
  g[1]->SetMarkerColor(2);
  g[1]->SetMarkerStyle(21);

  TF1 *fProd[2][2];

  fProd[0][0] = prova.GetSPfunction(0); 
  fProd[0][1] = prova.GetSPfunction(1); 
  fProd[1][0] = prova.GetSPfunction(2); 
  fProd[1][1] = prova.GetSPfunction(3); 


  fProd[0][0]->SetRange(-6,6);
  fProd[0][1]->SetRange(-6,6);
  fProd[1][0]->SetRange(-6,6);
  fProd[1][1]->SetRange(-6,6);


  hprob[0]->Draw("SAME");
  hprob[1]->Draw("SAME");

  hprob[0]->SetLineColor(4);
  hprob[1]->SetLineColor(2);
  hprob[0]->SetMarkerColor(4);
  hprob[1]->SetMarkerColor(2);

  hprob[0]->SetMarkerStyle(25);
  hprob[1]->SetMarkerStyle(25);

  hprob[0]->SetMarkerSize(0.5);
  hprob[1]->SetMarkerSize(0.5);

  TLegend *leg4 = new TLegend(0.75,0.64,0.9,0.85);
  leg4->SetHeader("Bayesian");
  leg4->AddEntry(hprob[0],"P_{1}(S)","P");
  leg4->AddEntry(hprob[1],"P_{2}(S)","P");
  leg4->Draw("SAME");
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);

  c->Print("figGaus2.png");
  c->Print("figGaus2.eps");

  TCanvas *c2 = new TCanvas();
  c2->Divide(2,2);
  c2->cd(1);
  fProd[0][0]->Draw();
  TPaveText t1(0.5,0.4,5.5,0.6);
  t1.SetBorderSize(0);
  t1.AddText("p_{1} #psi_{1} #rightarrow #int_{} p_{2} #psi_{1} = 1");
  t1.SetFillStyle(0);
  t1.Draw("SAME");
  fProd[0][0]->SetMaximum(0.65);
  fProd[0][0]->SetMinimum(-0.2);
  c2->cd(2);
  fProd[0][1]->Draw();
  TPaveText t2(0.5,0.4,5.5,0.6);
  t2.SetBorderSize(0);
  t2.AddText("p_{2} #psi_{1} #rightarrow #int_{} p_{2} #psi_{1} = 0");
  t2.SetFillStyle(0);
  t2.Draw("SAME");
  fProd[0][1]->SetMaximum(0.65);
  fProd[0][1]->SetMinimum(-0.2);
  c2->cd(3);
  fProd[1][0]->Draw();
  TPaveText t3(-5.5,0.4,-0.5,0.6);
  t3.SetBorderSize(0);
  t3.AddText("p_{1} #psi_{2} #rightarrow #int_{} p_{1} #psi_{2} = 0");
  t3.SetFillStyle(0);
  t3.Draw("SAME");
  fProd[1][0]->SetMaximum(0.65);
  fProd[1][0]->SetMinimum(-0.2);
  c2->cd(4);
  fProd[1][1]->Draw();
  TPaveText t4(-5.5,0.4,-0.5,0.6);
  t4.SetBorderSize(0);
  t4.AddText("p_{2} #psi_{2} #rightarrow #int_{} p_{2} #psi_{2} = 1");
  t4.SetFillStyle(0);
  t4.Draw("SAME");
  fProd[1][1]->SetMaximum(0.65);
  fProd[1][1]->SetMinimum(-0.2);

  c2->Update();

  for(Int_t i=0;i<2;i++){
    for(Int_t j=0;j<2;j++){
      fProd[i][j]->GetXaxis()->SetTitleSize(0.05);
      fProd[i][j]->GetYaxis()->SetTitleSize(0.05);
      fProd[i][j]->GetXaxis()->SetNdivisions(408);
      fProd[i][j]->GetYaxis()->SetNdivisions(408);
    }
  }

  c2->Update();

  fProd[0][0]->SetTitle(";S;d(p_{1}#psi_{1})/dS");
  fProd[0][1]->SetTitle(";S;d(p_{2}#psi_{1})/dS");
  fProd[1][0]->SetTitle(";S;d(p_{1}#psi_{2})/dS");
  fProd[1][1]->SetTitle(";S;d(p_{2}#psi_{2})/dS");

  TCanvas *c3 = new TCanvas();
  c3->SetBottomMargin(0.15);
  c3->SetTopMargin(0.05);
  c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.05);
  hsig->Draw();
  hsig->SetStats(0);
  hsig->SetTitle(";S;input entries");
  hsig->GetXaxis()->SetTitleSize(0.05);
  hsig->GetXaxis()->SetNdivisions(408);
  hsig->GetYaxis()->SetTitleSize(0.05);
  hsig->GetYaxis()->SetNdivisions(408);
  hsig->GetYaxis()->SetTitleOffset(1.2);

  hsim1->Draw("SAME");
  hsim2->Draw("SAME");
  hsig->SetLineWidth(2);
  hsim1->SetLineWidth(2);
  hsim2->SetLineWidth(2);

  TLegend *leg3 = new TLegend(0.63,0.63,0.93,0.93);
  leg3->SetHeader("Simulation input");
  leg3->Draw("SAME");
  leg3->SetFillStyle(0);
  leg3->AddEntry(hsig,"Total distribution","l");
  leg3->AddEntry(hsim1,"Type 1 signals","l");
  leg3->AddEntry(hsim2,"Type 2 signals","l");

  c3->Print("figInput.png");
  c3->Print("figInput.eps");

  TCanvas *c4 = new TCanvas();
  c4->SetBottomMargin(0.15);
  c4->SetTopMargin(0.05);
  c4->SetLeftMargin(0.15);
  c4->SetRightMargin(0.05);

  g[0]->Draw("AP");
  g[1]->Draw("P");
  g[0]->SetMaximum(1);
  g[0]->SetMinimum(0);
  g[0]->SetTitle(";N_{steps} Bayesian approach;Relative abundances");
  g[0]->GetXaxis()->SetTitleSize(0.05);
  g[0]->GetXaxis()->SetNdivisions(408);
  g[0]->GetYaxis()->SetTitleSize(0.05);
  g[0]->GetYaxis()->SetNdivisions(408);

  TLine *l[2];
  l[0] = new TLine(1,recocounts[0]*1E-6,20,recocounts[0]*1E-6);
  l[1] = new TLine(1,recocounts[1]*1E-6,20,recocounts[1]*1E-6);

  l[0]->SetLineColor(4);
  l[0]->SetLineWidth(3);
  l[0]->SetLineStyle(2);
  l[1]->SetLineColor(2);
  l[1]->SetLineWidth(3);
  l[1]->SetLineStyle(2);

  l[0]->Draw("SAME");
  l[1]->Draw("SAME");

  TLegend *leg2 = new TLegend(0.2,0.75,0.9,0.93);
  leg2->SetHeader("#Delta_{12} = 2");
  leg2->Draw("SAME");
  leg2->SetFillStyle(0);
  leg2->AddEntry(g[0],"1) Bayesian iterative","lp");
  leg2->AddEntry(g[1],"2) Bayesian iterative","lp");
  leg2->AddEntry(l[0],"1) QM approach","l");
  leg2->AddEntry(l[1],"2) QM approach","l");
  leg2->SetNColumns(2);

  c4->Print("figIterative.png");
  c4->Print("figIterative.eps");

  TFile *fout = new TFile("gaus2.root","RECREATE");
  hsig->Write();
  hampl[0]->Write();
  hampl[1]->Write();
  hsig1->Write();
  hsig2->Write();
  hsim1->Write();
  hsim2->Write();
  c->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  f1->Write();
  f2->Write();
  ff1->Write();
  ff2->Write();
  fProd[0][0]->Write();
  fProd[0][1]->Write();
  fProd[1][0]->Write();
  fProd[1][1]->Write();
  hprob[0]->Write();
  hprob[1]->Write();
  hx[0]->Write();
  hx[1]->Write();
  htrend[0]->Write();
  htrend[1]->Write();
  hxy->Write();
  hxyA[0]->Write();
  hxyA[1]->Write();
  hxyS[0]->Write();
  hxyS[1]->Write();
  hpriors[0][19]->Write();
  hpriors[1][19]->Write();
  htrendPriors[0][19]->Write();
  htrendPriors[1][19]->Write();
  hxyPriors[0][19]->Write();
  hxyPriors[1][19]->Write();
  fout->Close();

  return 0;
}
