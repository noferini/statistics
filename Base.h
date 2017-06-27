#ifndef BASE
#define BASE

#include "TObject.h"
#include "TF1.h"
#include "TMatrix.h"

class Base:public TObject{
 public:
  Base(Int_t n=0);
  ~Base();

  Int_t GetNtype() const {return fNtype;}
  virtual void SetNtype(Int_t n);

  virtual Double_t EvaluateProb(Int_t itype,Float_t signal) = 0;
  virtual TF1 *GetProbabilityDensity(Int_t itype) = 0;
  virtual TF1 *GetProbabilityDensityStar(Int_t itype) = 0;
  virtual void SetResponseFunction(Int_t itype,TObject *response) = 0;
  virtual void SetMatrix() = 0;
  virtual void Print() const = 0;
  virtual Double_t ScalarProduct(Int_t itype1,Int_t itype2) = 0;

  TMatrix& GetMatrix() {return *fMatrix;}
  TMatrix& GetInvMatrix() {return *fInvMatrix;}

  Float_t *GetAbundances() {return fAbundances;}
  void SetAbundances(Int_t itype,Float_t val) {fAbundances[itype]=val;}

 private:
  Int_t fNtype; // number of types 
  TMatrix *fMatrix; // matrix with scalar products 
  TMatrix *fInvMatrix; // inverse matrix
  Float_t *fAbundances; // simulated abundances

  ClassDef(Base,1);
};

#endif
