#include "Hough.h"

void HoughTransform( std::vector<Double_t> &f_X,   std::vector<Double_t> &f_Y,
		     std::vector<Double_t> &f_afX, std::vector<Double_t> &f_afY ){
  Int_t nstep = 180;

  for( Int_t ihit=0; ihit<f_X.size(); ihit++ ){
    for( Int_t istep=0; istep<nstep; istep++ ){
      if( f_Y[ihit]<-100 ) continue;
      if( f_Y[ihit]> 100 ) continue;
      f_afY.push_back(f_X[ihit]*180/TMath::Pi()*TMath::Cos((double)istep*TMath::Pi()/180)
		      +f_Y[ihit]*TMath::Sin((double)istep*TMath::Pi()/180));
      f_afX.push_back(istep);
    }
  }
  return;
}

void HoughFit_One( TH2D *f_hist, Double_t &f_par0, Double_t &f_par1){
  Int_t max_xbin, max_ybin, max_zbin;
  f_hist->GetMaximumBin(max_xbin,max_ybin, max_zbin);

  Double_t Rtheta = f_hist->GetXaxis()->GetBinLowEdge(max_xbin) + f_hist->GetXaxis()->GetBinWidth(max_xbin)/2.0;
  Double_t Rr     = f_hist->GetYaxis()->GetBinLowEdge(max_ybin) + f_hist->GetXaxis()->GetBinWidth(max_ybin)/2.0;
  f_hist->SetTitle( Form("%s:Max(%.1f,%.1f)",f_hist->GetTitle(),Rtheta,Rr) );
  
  f_par0 = Rr/TMath::Sin(Rtheta*TMath::Pi()/180);
  f_par1 = -1/TMath::Tan(Rtheta*TMath::Pi()/180);

  return;
}

void GetPol1Residual(TF1 *f_func, std::vector<Double_t> &f_X, std::vector<Double_t> &f_Y, std::vector<Double_t>& f_resi ){
  Double_t par0 = f_func->GetParameter(0);
  Double_t par1 = f_func->GetParameter(1);
  for( Int_t ivec=0; ivec<f_X.size(); ivec++ ) f_resi.push_back( f_Y[ivec]-(par0+par1*f_X[ivec]) );

  return;   
}

void GetPol1ResidualY(TF1 *f_func, std::vector<Double_t> &f_X, std::vector<Double_t> &f_Y, std::vector<Double_t>& f_resi ){
  Double_t par0 = f_func->GetParameter(0);
  Double_t par1 = f_func->GetParameter(1);
  
  for( Int_t ivec=0; ivec<f_X.size(); ivec++ ) f_resi.push_back( (f_X[ivec]-((f_Y[ivec]-par0)/par1)) *180/TMath::Pi() );

  return; 
}
