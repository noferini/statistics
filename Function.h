#ifndef FUNCTION
#define FUNCTION
#include "Base.h"

/*
  Assuming general functions
 */

class Function: public Base{
 public:

  Function(Int_t n=0);
  ~Function();

  void SetNtype(Int_t n);

  Double_t EvaluateProb(Int_t itype,Float_t signal);
  TF1 *GetProbabilityDensity(Int_t itype);
  TF1 *GetProbabilityDensityStar(Int_t itype);
  void SetResponseFunction(Int_t itype,TObject *response); // pass a TF1 gaus function (mean, sigma already set)
  void SetMatrix();
  void Print() const;
  Double_t ScalarProduct(Int_t itype1,Int_t itype2);
  Double_t ScalarProductOrFunc(Int_t itype1,Int_t itype2);

  Double_t GetIntegral(TF1 *f);

  TF1 *GetSPfunction(Int_t ival) {return fFuncSP[ival];}

 private:
  TF1 **fFunc; // functions Original
  TF1 **fFuncStar; // functions Transponed
  TF1 **fFuncSPor; // functions scalar products
  TF1 **fFuncSP; // functions scalar products Origin

  ClassDef(Function,1);
};
#endif
