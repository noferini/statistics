#include <stdio.h>

#include "Gaus.h"
#include "Function.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TLegend.h"

int main(){
  Gaus prova;
  const Int_t nfunc = 2;
  prova.SetNtype(nfunc);
  
  Float_t abundances[nfunc] = {0.7,1};

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

  f1->SetLineColor(1);
  f2->SetLineColor(2);

  ff1->SetLineColor(1);
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
  leg->AddEntry(ff2,"#psi_{1}(S)","l");

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

  c->Print("figGaus2.png");
  c->Print("figGaus2.eps");

  Float_t val,signal;
  Int_t isp;

  TF1 *rf[nfunc] = {f1,f2};

  TH1F *hsig = new TH1F("hsig",";signal",1000,-20,20);
  TH1F *hsig1 = new TH1F("hsig1",";signal",1000,-20,20);
  TH1F *hsig2 = new TH1F("hsig2",";signal",1000,-20,20);

  TH1F *hsim1 = new TH1F("hsim1",";signal",1000,-20,20);
  TH1F *hsim2 = new TH1F("hsim2",";signal",1000,-20,20);

  hsig1->SetLineColor(6);
  hsig2->SetLineColor(2);

  hsim1->SetLineColor(6);
  hsim2->SetLineColor(2);

  hsim1->SetLineStyle(2);
  hsim2->SetLineStyle(2);

  TProfile *hampl[nfunc];
  hampl[0] = new TProfile("hampl0",";signal;amplitude",1000,-20,20);
  hampl[1] = new TProfile("hampl1",";signal;amplitude",1000,-20,20);

  Int_t counts[nfunc] = {0,0};
  Float_t recocounts[nfunc] = {0,0};


  Double_t ampl[nfunc],norm;

  for(Int_t i=0;i < 100000;i++){
    val = gRandom->Rndm();
    
    if(val < abundances[0]) isp=0;
    else isp=1;

    counts[isp]++;

    signal = rf[isp]->GetRandom();

    ampl[0] = prova.EvaluateProb(0,signal);
    ampl[1] = prova.EvaluateProb(1,signal);

    hampl[0]->Fill(signal,ampl[0]);
    hampl[1]->Fill(signal,ampl[1]);

    recocounts[0] += ampl[0];
    recocounts[1] += ampl[1];
  }

  printf("counts[0] = %i -- recocounts[0] = %f -- error[%c] = %f\n",counts[0],recocounts[0],'%',(recocounts[0]/counts[0] -1)*100);
  printf("counts[1] = %i -- recocounts[1] = %f -- error[%c] = %f\n",counts[1],recocounts[1],'%',(recocounts[1]/counts[1] -1)*100);

  for(Int_t i=0;i < 100000;i++){
    val = gRandom->Rndm();
    
    if(val < abundances[0]) isp=0;
    else isp=1;

    signal = rf[isp]->GetRandom();

    hsig->Fill(signal);
    if(isp==0) hsim1->Fill(signal);
    if(isp==1) hsim2->Fill(signal);

    ampl[0] = rf[0]->Eval(signal)*recocounts[0];
    ampl[1] = rf[1]->Eval(signal)*recocounts[1];

    norm = ampl[0]+ampl[1] + 1E-10;

    ampl[0] /= norm;
    ampl[1] /= norm;

    hsig1->Fill(signal,ampl[0]);
    hsig2->Fill(signal,ampl[1]);

  }

  TFile *fout = new TFile("gaus2.root","RECREATE");
  hsig->Write();
  hampl[0]->Write();
  hampl[1]->Write();
  hsig1->Write();
  hsig2->Write();
  hsim1->Write();
  hsim2->Write();
  c->Write();
  f1->Write();
  f2->Write();
  ff1->Write();
  ff2->Write();

  fout->Close();

  return 0;
}
