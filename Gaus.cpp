#include "Gaus.h"

ClassImp(Gaus);

void Gaus::SetNtype(Int_t n){
  Base::SetNtype(n);

}

Double_t Gaus::EvaluateProb(Int_t itype,Float_t signal){}

TF1 *Gaus::GetProbabilityDensity(Int_t itype){}

void Gaus::SetResponseFunction(Int_t itype,TObject *response){
  if(itype >= GetNtype() || itype < 0) return;
  TF1 *myfunc = (TF1 *) response;


}

void Gaus::SetMatrix(){}
