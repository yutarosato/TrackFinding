#include "root.hh"

void hesscol(){
   UInt_t Number = 5;
   Double_t Red[5]   = { 0.00, 0.00, 0.60, 1.00, 1.00 };
   Double_t Green[5] = { 0.00, 0.00, 0.00, 0.00, 1.00 };
   Double_t Blue[5]  = { 0.00, 0.60, 0.60, 0.00, 0.00 };
   Double_t Stops[5] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
   //TStyle::CreateGradientColorTable(Number,Stops, Red, Green, Blue, 256);
   TColor::CreateGradientColorTable(Number,Stops, Red, Green, Blue, 256);

}

void cpalette(){
  const Int_t   colNum = 32;
  Int_t         palette[colNum+1];
  Float_t       tmp;
  Int_t         n = colNum/4;
 
  Float_t       high = 1.0;
 
  for(Int_t i=0; i<n; i++){
    tmp = ((float) i/n)*high;

    if(! gROOT->GetColor(230+i)){
      //      TColor *color = new TColor(230+i, 0, tmp, high, "");
    }
    else{
      TColor *color = gROOT->GetColor(230+i);
      color->SetRGB(0, tmp, high);
    }

    palette[i] = 230+i;
  }

  for(Int_t i=n; i<2*n; i++){
    tmp = (1 - (float)(i-n)/n)*high;

    if(! gROOT->GetColor(230+i)){
      //TColor *color = new TColor(230+i, 0, high, tmp, "");
    }
    else{
      TColor *color = gROOT->GetColor(230+i);
      color->SetRGB(0, high, tmp);
    }

    palette[i] = 230+i;
  }

  for(Int_t i=2*n; i<3*n; i++){
    tmp = ((float)(i-2*n)/n)*high;

    if(! gROOT->GetColor(230+i)){
      //TColor *color = new TColor(230+i, tmp, high, 0, "");
    }
    else{
      TColor *color = gROOT->GetColor(230+i);
      color->SetRGB(tmp, high, 0);
    }

    palette[i] = 230+i;
  }

  for(Int_t i=3*n; i<=4*n; i++){
    tmp = (1 - (float)(i-3*n)/n)*high;

    if(! gROOT->GetColor(230+i)){
      //TColor *color = new TColor(230+i, high, tmp, 0, "");
    }
    else{
      TColor *color = gROOT->GetColor(230+i);
      color->SetRGB(high, tmp, 0);
    }

    palette[i] = 230+i;
  }

  gStyle->SetPalette(colNum, palette);
}
