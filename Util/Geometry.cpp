#include "Geometry.h"

Int_t GetVaneID( Double_t f_phi ){
  ///*
  Int_t vaneID = (Int_t)(f_phi/(2.0*TMath::Pi()/n_vane));
  if( vaneID==n_vane ) vaneID = 0;
  return vaneID;
  //*/
  /* uk codes
  Double_t tmpVane=(f_phi/(2*TMath::Pi()/n_vane));
  Int_t VNum = (Int_t)tmpVane;
  if(tmpVane-VNum >=0.5){
    VNum++;
    if(VNum==n_vane) VNum=0;
  }
  return VNum;
  */
}
