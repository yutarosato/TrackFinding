 void SetHESSColor(){
   UInt_t Number = 5;
   Double_t Red[5]   = { 0.00, 0.00, 0.60, 1.00, 1.00 };
   Double_t Green[5] = { 0.00, 0.00, 0.00, 0.00, 1.00 };
   Double_t Blue[5]  = { 0.00, 0.60, 0.60, 0.00, 0.00 };
   Double_t Stops[5] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
   //  TStyle::CreateGradientColorTable(Number,Stops, Red, Green, Blue, 256);
   TColor::CreateGradientColorTable(Number,Stops, Red, Green, Blue, 256);
   return;
 }
