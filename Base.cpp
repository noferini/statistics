#include "Base.h"

ClassImp(Base);

Base::Base(Int_t n):
  fNtype(n),
  fMatrix(NULL),
  fInvMatrix(NULL),
  fAbundances(NULL)
{
  if(n!=0) SetNtype(n);

}

Base::~Base(){
  if(fMatrix) delete fMatrix;
  if(fInvMatrix) delete fInvMatrix;
  if(fAbundances) delete fAbundances;
}

void Base::SetNtype(Int_t n){
  fNtype=n;

  if(fMatrix) delete fMatrix;
  if(fInvMatrix) delete fInvMatrix;
  if(fAbundances) delete fAbundances;

  fMatrix = new TMatrix(n,n);
  fInvMatrix = new TMatrix(n,n);
  fAbundances = new Float_t(n);
}
