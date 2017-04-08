#ifndef HITSARRAY_HH
#define HITSARRAY_HH

#include "const.h"

class HitsArray{

public:
  HitsArray();
  ~HitsArray();
  void InputHits( Int_t index, Double_t x, Double_t y, Double_t z, Double_t pt, Double_t gt, Int_t pID, Double_t EachDepE );
  void Print( Int_t fl_message );
  void CalcOrder();
  
  private:
  std::vector<Int_t>    m_Index;
  std::vector<Double_t> m_X;
  std::vector<Double_t> m_Y;
  std::vector<Double_t> m_Z;
  std::vector<Double_t> m_R;
  std::vector<Double_t> m_Phi;
  std::vector<Int_t>    m_VaneID;
  std::vector<Double_t> m_gT;
  std::vector<Double_t> m_pT;
  std::vector<Int_t>    m_pID;
  std::vector<Double_t> m_EachDepE;
  //std::vector<Double_t> m_residual_PhiZ;

  std::vector<Int_t>    m_order_VaneID;
  std::vector<Int_t>    m_order_gT;

  std::vector<Double_t> m_zphi_par0;
  std::vector<Double_t> m_zphi_par1;

  Double_t Phi_uk   ( Double_t y, Double_t x );
  Int_t    GetVaneID( Double_t f_phi );
};

#endif
