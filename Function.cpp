#include "Function.h"
#include <stdio.h>
#include "TMath.h"
// #include "Math/GSLIntegrator.h"
// #include "Math/WrappedTF1.h"

ClassImp(Function);

Function::Function(Int_t n):
  Base(n),
  fFunc(NULL),
  fFuncStar(NULL),
  fFuncSP(NULL),
  fFuncSPor(NULL)
{
  if(n!=0) SetNtype(n);
}

Function::~Function(){
  if(fFunc) delete fFunc;
  if(fFuncStar) delete fFuncStar;
  if(fFuncSP) delete fFuncSP;
  if(fFuncSPor) delete fFuncSPor;
}

void Function::SetNtype(Int_t n){
  Base::SetNtype(n);

  if(fFunc) delete fFunc;
  if(fFuncStar) delete fFuncStar;
  if(fFuncSP) delete fFuncSP;
  if(fFuncSPor) delete fFuncSPor;

  fFunc = new TF1*[n];
  fFuncStar = new TF1*[n];
  fFuncSP = new TF1*[n*n];
  fFuncSPor = new TF1*[n*n];

  for(Int_t i=0; i < n;i++){
    fFunc[i] = NULL;
    fFuncStar[i] = NULL;
    for(Int_t j=0; j < n;j++){
      fFuncSP[i*3+j] = NULL;
      fFuncSPor[i*3+j] = NULL;
    }
  }

}

Double_t Function::EvaluateProb(Int_t itype,Float_t signal){
  return GetProbabilityDensityStar(itype)->Eval(signal);
}

TF1 *Function::GetProbabilityDensity(Int_t itype){
  return fFunc[itype];
}

TF1 *Function::GetProbabilityDensityStar(Int_t itype){
  return fFuncStar[itype];
}

Double_t Function::ScalarProductOrFunc(Int_t itype1,Int_t itype2) {
  if(!fFuncSPor[itype1*GetNtype() + itype2]){
    TString function("(");
    function.Append(GetProbabilityDensity(itype1)->GetExpFormula("par"));
    function.Append(") * (");
    function.Append(GetProbabilityDensity(itype2)->GetExpFormula("par"));
    function.Append(")");

    TString namefunc(GetProbabilityDensity(itype1)->GetName());
    namefunc.Append("_");
    namefunc.Append(GetProbabilityDensity(itype2)->GetName());

    fFuncSPor[itype1*GetNtype() + itype2] = new TF1(namefunc.Data(),function.Data(),-20,20);


    for(Int_t i=0;i < GetProbabilityDensity(itype1)->GetNpar();i++)
      fFuncSPor[itype1*GetNtype() + itype2]->SetParameter(i,GetProbabilityDensity(itype1)->GetParameter(i));

    for(Int_t i=0;i < GetProbabilityDensity(itype2)->GetNpar();i++)
      fFuncSPor[itype1*GetNtype() + itype2]->SetParameter(i+GetProbabilityDensity(itype1)->GetNpar(),GetProbabilityDensity(itype2)->GetParameter(i));
    
    fFuncSPor[itype1*GetNtype() + itype2]->SetNpx(10000);

  }
  
  return GetIntegral(fFuncSPor[itype1*GetNtype() + itype2]);
}

Double_t Function::ScalarProduct(Int_t itype1,Int_t itype2) {
  if(!fFuncSP[itype1*GetNtype() + itype2]){

    return 0;
  }
  return GetIntegral(fFuncSP[itype1*GetNtype() + itype2]);
}

void Function::SetResponseFunction(Int_t itype,TObject *response){
  if(itype >= GetNtype() || itype < 0) return;
  fFunc[itype] = (TF1 *) response;
}

void Function::SetMatrix(){
  printf("Set Matrix with Ntype = %i\n",GetNtype());


  for(Int_t i=0;i < GetNtype();i++){
    for(Int_t j=0;j < GetNtype();j++){
      GetMatrix()[i][j] = ScalarProductOrFunc(i,j);
      GetInvMatrix()[i][j] = GetMatrix()[i][j];
    }
  }
  GetMatrix().Print();

  GetInvMatrix().Invert();

  GetInvMatrix().Print();

  for(Int_t i=0;i < GetNtype();i++){
    TString namefunc("T_");
    namefunc.Append(GetProbabilityDensity(i)->GetName());

    TString function("");
    for(Int_t j=0;j < GetNtype();j++){
      if(j >0 && GetInvMatrix()[i][j] >= 0)  function.Append("+");
      function.Append(Form("%f*(",GetInvMatrix()[i][j]));
      function.Append(GetProbabilityDensity(j)->GetExpFormula("par"));
      function.Append(")");
    }

    fFuncStar[i] = new TF1(namefunc.Data(),function.Data(),-20,20);

    for(Int_t j=0;j < GetNtype();j++){
      TString namefunc2("SP_");
      TString function2("(");
      function2.Append(function.Data());
      function2.Append(")");

      namefunc2.Append(GetProbabilityDensityStar(i)->GetName());
      namefunc2.Append("_");
      namefunc2.Append(GetProbabilityDensity(j)->GetName());

      function2.Append("*(");

      //function2.Append(GetProbabilityDensity(j)->GetName());
      function2.Append(GetProbabilityDensity(j)->GetExpFormula("par"));
      
      function2.Append(")");

      fFuncSP[i*GetNtype() + j] = new TF1(namefunc2.Data(),function2.Data(),-20,20);
    }
  }
}

void Function::Print() const {
  printf("Function class\n");
  for(Int_t i=0;i < GetNtype();i++)
    fFunc[i]->Print();

}

Double_t Function::GetIntegral(TF1 *f){
  Double_t *param = 0;

//   ROOT::Math::GSLIntegrator ig(1.E-9,1.E-9,10000);
//   ROOT::Math::WrappedTF1 wf(*f);
//   ig.SetFunction(wf);
//   double val2 = ig.Integral(-20,20);
  return f->Integral(-20,20);//,param,1e-10);
}
