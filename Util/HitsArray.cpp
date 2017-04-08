#include "HitsArray.h"

HitsArray::HitsArray(){};
HitsArray::~HitsArray(){};

void HitsArray::InputHits( Int_t index, Double_t x, Double_t y, Double_t z, Double_t pt, Double_t gt, Int_t pID, Double_t EachDepE ){
  m_Index.push_back ( index                   );
  m_X.push_back     ( x                       );
  m_Y.push_back     ( y                       );
  m_Z.push_back     ( z                       );
  m_R.push_back     ( sqrt(pow(x,2)+pow(y,2)) );
  m_Phi.push_back   ( Phi_uk(y,x)             );
  m_VaneID.push_back( GetVaneID(Phi_uk(y,x))  );
  m_gT.push_back    ( gt                      );
  m_pT.push_back    ( pt                      );
  m_pID.push_back   ( pID                     );
  m_EachDepE.push_back(EachDepE);
}

void HitsArray::Print( Int_t fl_message ){
  if( fl_message < 2 ) return;
  for( Int_t ihit=0; ihit<m_X.size(); ihit++ ){
    std::cout << "              "
	      << std::setw(3) << std::right << m_Index.at(ihit) << " : pID = "
	      << m_pID.at(ihit)   << ", (x,y,z) = ("
	      << std::setw(9) << std::right << Form("%.3f",m_X.at       (ihit)) << ", "
	      << std::setw(9) << std::right << Form("%.3f",m_Y.at       (ihit)) << ", "
	      << std::setw(9) << std::right << Form("%.3f",m_Z.at       (ihit)) << ", phi = "
	      << std::setw(6) << std::right << Form("%.3f",m_Phi.at     (ihit)) << ", vane-ID = "
	      << std::setw(2) << std::right << Form("%d",  m_VaneID.at  (ihit)) << ", t(proper) = "
	      << std::setw(7) << std::right << Form("%.5f",m_pT.at      (ihit)) << ", t(global) = "
	      << std::setw(8) << std::right << Form("%.5f",m_gT.at      (ihit)) << ", Edep = "
	      << std::setw(7) << std::right << Form("%.5f",m_EachDepE.at(ihit))
	      << std::endl;
  }
  return;
}

void HitsArray::CalcOrder(){
  std::multimap<Int_t,Double_t> tMap_VaneID;
  std::multimap<Int_t,Double_t> tMap_gT;
  for( Int_t ivec=0; ivec< m_Index.size(); ivec++ ){
    tMap_VaneID.insert( std::make_pair(m_VaneID.at(ivec),m_Index.at(ivec)) );
    tMap_gT.insert    ( std::make_pair(m_VaneID.at(ivec),m_Index.at(ivec)) );
  }
  std::multimap<Int_t,Double_t>::iterator it_VaneID = tMap_VaneID.begin();
  std::multimap<Int_t,Double_t>::iterator it_gT     = tMap_gT.begin    ();
  while( it_VaneID != tMap_VaneID.end() ){ m_order_VaneID.push_back( (*it_VaneID).second ); it_VaneID++; }
  while( it_gT     != tMap_gT.end    () ){ m_order_gT.push_back    ( (*it_gT    ).second ); it_gT++;     }
}

Double_t HitsArray::Phi_uk(Double_t y, Double_t x){ // 0~phi~2pi : this definition is used in Ueno-san's codes
  Double_t phi = TMath::ATan2(y,x);
  if( phi<0 ) phi += 2*TMath::Pi();
  return phi;
}

Int_t HitsArray::GetVaneID( Double_t f_phi ){
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
