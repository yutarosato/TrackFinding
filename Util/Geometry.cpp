#include "Geometry.h"

Int_t GetVaneID( Double_t f_phi ){
  Int_t vaneID;
  Double_t tmpVane= f_phi/(2.0*TMath::Pi()/n_vane);
  vaneID = (Int_t)(tmpVane+0.5);
  if( vaneID==n_vane ) vaneID = 0;
  return vaneID;
}
