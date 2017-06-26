#ifndef GAUS
#define GAUS
#include "Base.h"

/*
  Assuming Gaussian function with mean and sigma varying with itype
 */

class Gaus: public Base{
 public:

  void SetNtype(Int_t n);

  Double_t EvaluateProb(Int_t itype,Float_t signal);
  TF1 *GetProbabilityDensity(Int_t itype);
  void SetResponseFunction(Int_t itype,TObject *response); // pass a TF1 gaus function (mean, sigma already set)
  void SetMatrix();

 private:
  ClassDef(Gaus,1);
};
#endif
