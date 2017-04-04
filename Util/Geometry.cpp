#include "Geometry.h"

Int_t GetVaneID( Double_t f_phi ){
  Int_t vaneID = (Int_t)(f_phi/(2.0*TMath::Pi()/n_vane));
  if( vaneID==n_vane ) vaneID = 0;

  return vaneID;
}
