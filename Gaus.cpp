#include "Gaus.h"
#include <stdio.h>
#include "TMath.h"

ClassImp(Gaus);

Gaus::Gaus(Int_t n):
  Base(n),
  fMean(NULL),
  fSigma(NULL)
{
  if(n!=0) SetNtype(n);
}

Gaus::~Gaus(){
  if(fMean) delete fMean;
  if(fSigma) delete fSigma;

}

void Gaus::SetNtype(Int_t n){
  Base::SetNtype(n);

  if(fMean) delete fMean;
  if(fSigma) delete fSigma;

  fMean = new Float_t(n);
  fSigma = new Float_t(n);
}

Double_t Gaus::EvaluateProb(Int_t itype,Float_t signal){
  return TMath::Exp(-(signal-fMean[itype])*(signal-fMean[itype])*0.5/fSigma[itype]/fSigma[itype])/fSigma[itype]/TMath::Sqrt(2*TMath::Pi());
}

TF1 *Gaus::GetProbabilityDensity(Int_t itype){
  TF1 *f = new TF1(Form("f%i",itype),"gaus",-10,10);
  f->SetParameter(0,1./fSigma[itype]/TMath::Sqrt(2*TMath::Pi()));
  f->SetParameter(1,fMean[itype]);
  f->SetParameter(2,fSigma[itype]);
  return f;
}

TF1 *Gaus::GetProbabilityDensityStar(Int_t itype){
  TString function("gaus(0)");

  for(Int_t i=1;i < GetNtype();i++)
    function.Append(Form("+gaus(%i)",i*3));
  
  TF1 *f = new TF1(Form("fstar%i",itype),function.Data(),-10,10);

  for(Int_t i=0;i < GetNtype();i++){
    f->SetParameter(i*3+0,1./fSigma[i]/TMath::Sqrt(2*TMath::Pi())*GetInvMatrix()[i][itype]);
    f->SetParameter(i*3+1,fMean[i]);
    f->SetParameter(i*3+2,fSigma[i]);
  }

  return f;
}

Double_t Gaus::ScalarProduct(Int_t itype1,Int_t itype2) {
  TString function("(gaus(0)");

  for(Int_t i=1;i < GetNtype();i++)
    function.Append(Form("+gaus(%i)",i*3));

  function.Append(Form(")*gaus(%i)",GetNtype()*3));
  
  TF1 *f = new TF1(Form("fSP%i_%i",itype1,itype2),function.Data(),-10,10);

  for(Int_t i=0;i < GetNtype();i++){
    f->SetParameter(i*3+0,1./fSigma[i]/TMath::Sqrt(2*TMath::Pi())*GetInvMatrix()[i][itype1]);
    f->SetParameter(i*3+1,fMean[i]);
    f->SetParameter(i*3+2,fSigma[i]);
  }

  f->SetParameter(GetNtype()*3+0,1./fSigma[itype2]/TMath::Sqrt(2*TMath::Pi()));
  f->SetParameter(GetNtype()*3+1,fMean[itype2]);
  f->SetParameter(GetNtype()*3+2,fSigma[itype2]);
  
  return f->Integral(-10,10);
}

void Gaus::SetResponseFunction(Int_t itype,TObject *response){
  if(itype >= GetNtype() || itype < 0) return;
  TF1 *myfunc = (TF1 *) response;

  fMean[itype] = myfunc->GetParameter(1);
  fSigma[itype] = myfunc->GetParameter(2);
}

void Gaus::SetMatrix(){
  printf("Set Matrix with Ntype = %i\n",GetNtype());


  for(Int_t i=0;i < GetNtype();i++){
    for(Int_t j=0;j < GetNtype();j++){
      GetMatrix()[i][j] = TMath::Exp(-(fMean[i]-fMean[j])*(fMean[i]-fMean[j])/(fSigma[i]*fSigma[i] + fSigma[j]*fSigma[j])*0.5)/TMath::Sqrt(2*(fSigma[i]*fSigma[i] + fSigma[j]*fSigma[j])*TMath::Pi()); // to be corrected
      GetInvMatrix()[i][j] = GetMatrix()[i][j];
    }
  }
  GetInvMatrix().Invert();

  GetMatrix().Print();
  GetInvMatrix().Print();
}

void Gaus::Print() const {
  printf("Gaus class\n");
  for(Int_t i=0;i < GetNtype();i++)
    printf("type %i) mean = %f -- sigma = %f\n",i,fMean[i],fSigma[i]);

}
