#include <stdio.h>

#include "Gaus.h"
#include "TFile.h"
#include "TCanvas.h"

int main(){
  Gaus prova;
  prova.SetNtype(3);
  
  TF1 *ff = new TF1("ff","gaus");
  ff->SetParameter(1,-2);
  ff->SetParameter(2, 1);
  prova.SetResponseFunction(0,ff);
  ff->SetParameter(1,0);
  ff->SetParameter(2, 1);
  prova.SetResponseFunction(1,ff);
  ff->SetParameter(1,4);
  ff->SetParameter(2, 2);
  prova.SetResponseFunction(2,ff);

  TF1 *f1 = prova.GetProbabilityDensity(0);
  TF1 *f2 = prova.GetProbabilityDensity(1);
  TF1 *f3 = prova.GetProbabilityDensity(2);

  prova.SetMatrix();

  TF1 *ff1 = prova.GetProbabilityDensityStar(0);
  TF1 *ff2 = prova.GetProbabilityDensityStar(1);
  TF1 *ff3 = prova.GetProbabilityDensityStar(2);

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

  TFile *fout = new TFile("out.root","RECREATE");
  c->Write();
  f1->Write();
  f2->Write();
  f3->Write();
  ff1->Write();
  ff2->Write();
  ff3->Write();
  fout->Close();

  return 0;
}
