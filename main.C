#include <stdio.h>

#include "Gaus.h"
#include "Function.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TProfile.h"

int main(){
  Function prova;
  prova.SetNtype(3);
  
  Float_t abundances[3] = {0.7,0.9,1};

  TF1 *ff1 = new TF1("ff1","gaus",-20,20);
  ff1->SetParameter(0,1);
  ff1->SetParameter(1,-2);
  ff1->SetParameter(2, 1);
  ff1->SetParameter(0,1./ff1->Integral(-20,20));
  TF1 *ff2 = new TF1("ff2","gaus",-20,20);
  ff2->SetParameter(0,1);
  ff2->SetParameter(1,2);
  ff2->SetParameter(2, 1);
  ff2->SetParameter(0,1./ff2->Integral(-20,20));
  TF1 *ff3 = new TF1("ff3","gaus",-20,20);
  ff3->SetParameter(0,1);
  ff3->SetParameter(1,4);
  ff3->SetParameter(2, 2);
  ff3->SetParameter(0,1./ff3->Integral(-20,20));

  ff3 = new TF1("ff3","landau",-20,20);
  ff3->SetParameter(0,1);
  ff3->SetParameter(1,4);
  ff3->SetParameter(2, 0.5);
  ff3->SetParameter(0,1./ff3->Integral(-20,20));


  prova.SetResponseFunction(0,ff1);
  prova.SetResponseFunction(1,ff2);
  prova.SetResponseFunction(2,ff3);

  TF1 *f1 = prova.GetProbabilityDensity(0);
  TF1 *f2 = prova.GetProbabilityDensity(1);
  TF1 *f3 = prova.GetProbabilityDensity(2);

  prova.SetMatrix();

  ff1 = prova.GetProbabilityDensityStar(0);
  ff2 = prova.GetProbabilityDensityStar(1);
  ff3 = prova.GetProbabilityDensityStar(2);

  f1->SetLineColor(1);
  f2->SetLineColor(2);
  f3->SetLineColor(4);

  ff1->SetLineColor(1);
  ff2->SetLineColor(2);
  ff3->SetLineColor(4);

  ff1->SetLineStyle(2);
  ff2->SetLineStyle(2);
  ff3->SetLineStyle(2);

  for(Int_t i=0;i<3;i++)
    for(Int_t j=0;j<3;j++)
      printf("Scalar product <%i|%i> = %f\n",i,j,prova.ScalarProduct(i,j));

  TCanvas *c = new TCanvas();
  ff1->Draw();
  ff2->Draw("SAME");
  ff3->Draw("SAME");
  f1->Draw("SAME");
  f2->Draw("SAME");
  f3->Draw("SAME");

  Float_t val,signal;
  Int_t isp;

  TF1 *rf[3] = {f1,f2,f3};

  TH1F *hsig = new TH1F("hsig",";signal",1000,-20,20);
  TH1F *hsig1 = new TH1F("hsig1",";signal",1000,-20,20);
  TH1F *hsig2 = new TH1F("hsig2",";signal",1000,-20,20);
  TH1F *hsig3 = new TH1F("hsig3",";signal",1000,-20,20);

  TH1F *hsim1 = new TH1F("hsim1",";signal",1000,-20,20);
  TH1F *hsim2 = new TH1F("hsim2",";signal",1000,-20,20);
  TH1F *hsim3 = new TH1F("hsim3",";signal",1000,-20,20);

  hsig1->SetLineColor(6);
  hsig2->SetLineColor(2);
  hsig3->SetLineColor(4);

  hsim1->SetLineColor(6);
  hsim2->SetLineColor(2);
  hsim3->SetLineColor(4);

  hsim1->SetLineStyle(2);
  hsim2->SetLineStyle(2);
  hsim3->SetLineStyle(2);

  TProfile *hampl[3];
  hampl[0] = new TProfile("hampl0",";signal;amplitude",1000,-20,20);
  hampl[1] = new TProfile("hampl1",";signal;amplitude",1000,-20,20);
  hampl[2] = new TProfile("hampl2",";signal;amplitude",1000,-20,20);

  Int_t counts[3] = {0,0,0};
  Float_t recocounts[3] = {0,0,0};


  Double_t ampl[3],norm;

  for(Int_t i=0;i < 100000;i++){
    val = gRandom->Rndm();
    
    if(val < abundances[0]) isp=0;
    else if(val < abundances[1]) isp=1;
    else isp=2;

    counts[isp]++;

    signal = rf[isp]->GetRandom();

    ampl[0] = prova.EvaluateProb(0,signal);
    ampl[1] = prova.EvaluateProb(1,signal);
    ampl[2] = prova.EvaluateProb(2,signal);

    hampl[0]->Fill(signal,ampl[0]);
    hampl[1]->Fill(signal,ampl[1]);
    hampl[2]->Fill(signal,ampl[2]);

    recocounts[0] += ampl[0];
    recocounts[1] += ampl[1];
    recocounts[2] += ampl[2];
  }

  printf("counts[0] = %i -- recocounts[0] = %f -- error[%c] = %f\n",counts[0],recocounts[0],'%',(recocounts[0]/counts[0] -1)*100);
  printf("counts[1] = %i -- recocounts[1] = %f -- error[%c] = %f\n",counts[1],recocounts[1],'%',(recocounts[1]/counts[1] -1)*100);
  printf("counts[2] = %i -- recocounts[2] = %f -- error[%c] = %f\n",counts[2],recocounts[2],'%',(recocounts[2]/counts[2] -1)*100);

  for(Int_t i=0;i < 100000;i++){
    val = gRandom->Rndm();
    
    if(val < abundances[0]) isp=0;
    else if(val < abundances[1]) isp=1;
    else isp=2;

    signal = rf[isp]->GetRandom();

    hsig->Fill(signal);
    if(isp==0) hsim1->Fill(signal);
    if(isp==1) hsim2->Fill(signal);
    if(isp==2) hsim3->Fill(signal);

    ampl[0] = rf[0]->Eval(signal)*recocounts[0];
    ampl[1] = rf[1]->Eval(signal)*recocounts[1];
    ampl[2] = rf[2]->Eval(signal)*recocounts[2];

    norm = ampl[0]+ampl[1]+ampl[2] + 1E-10;

    ampl[0] /= norm;
    ampl[1] /= norm;
    ampl[2] /= norm;

    hsig1->Fill(signal,ampl[0]);
    hsig2->Fill(signal,ampl[1]);
    hsig3->Fill(signal,ampl[2]);

  }

  TFile *fout = new TFile("out.root","RECREATE");
  hsig->Write();
  hampl[0]->Write();
  hampl[1]->Write();
  hampl[2]->Write();
  hsig1->Write();
  hsig2->Write();
  hsig3->Write();
  hsim1->Write();
  hsim2->Write();
  hsim3->Write();
  c->Write();
  f1->Write();
  f2->Write();
  f3->Write();
  ff1->Write();
  ff2->Write();
  ff3->Write();
  prova.GetSPfunction(0)->Write();
  prova.GetSPfunction(4)->Write();
  prova.GetSPfunction(8)->Write();

  fout->Close();

  return 0;
}
