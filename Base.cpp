#include "Base.h"

ClassImp(Base);

Base::Base(Int_t n):
  fNtype(n),
  fMatrix(NULL),
  fInvMatrix(NULL)
{

}

Base::~Base(){
  if(fMatrix) delete fMatrix;
  if(fInvMatrix) delete fInvMatrix;
}

void Base::SetNtype(Int_t n){
  fNtype=n;
}
