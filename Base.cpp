#include "Base.h"

ClassImp(Base);

Base::Base(Int_t n):
  fNtype(n),
  fMatrix(NULL),
  fInvMatrix(NULL)
{
  if(n!=0) SetNtype(n);

}

Base::~Base(){
  if(fMatrix) delete fMatrix;
  if(fInvMatrix) delete fInvMatrix;
}

void Base::SetNtype(Int_t n){
  fNtype=n;

  if(fMatrix) delete fMatrix;
  if(fInvMatrix) delete fInvMatrix;

  fMatrix = new TMatrix(n,n);
  fInvMatrix = new TMatrix(n,n);
}
