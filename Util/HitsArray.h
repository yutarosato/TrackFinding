#ifndef HITSARRAY_HH
#define HITSARRAY_HH

#include "const.h"

class HitsArray{

public:
  HitsArray();
  ~HitsArray();

  void InputHits( Int_t index, Double_t x, Double_t y, Double_t z, Double_t pt, Double_t gt, Int_t pID, Double_t EachDepE );
  void ClearEvent();
  void Print              ( Int_t fl_message=2 );
  void Print_VaneID_Order ( Int_t fl_message=2 );
  void Print_gT_Order     ( Int_t fl_message=2 );
  void CalcOrder();
  void HoughTransform_phiz();
  void HoughFit_phiz();
  void CalcHoughResidual_phiz();

  void Clustering();

  Double_t Phi_uk   ( Double_t y, Double_t x );
  Int_t    GetVaneID( Double_t f_phi );
  Int_t    GetNVane      (){ return m_geom_nvane;             }
  Int_t    GetNhits      (){ return m_X.size();               }
  Int_t    GetNHoughLines(){ return m_hough_phiz_par0.size(); }

  Double_t GetX       ( Int_t index ){ if( index>=m_X.size       () ){ std::cerr << "[ABORT] Wrong index for X"        << std::endl, abort(); } return m_X.at       (index); }
  Double_t GetY       ( Int_t index ){ if( index>=m_Y.size       () ){ std::cerr << "[ABORT] Wrong index for Y"        << std::endl, abort(); } return m_Y.at       (index); }
  Double_t GetZ       ( Int_t index ){ if( index>=m_Z.size       () ){ std::cerr << "[ABORT] Wrong index for Z"        << std::endl, abort(); } return m_Z.at       (index); }
  Double_t GetR       ( Int_t index ){ if( index>=m_R.size       () ){ std::cerr << "[ABORT] Wrong index for R"        << std::endl, abort(); } return m_R.at       (index); }
  Double_t GetPhi     ( Int_t index ){ if( index>=m_Phi.size     () ){ std::cerr << "[ABORT] Wrong index for Phi"      << std::endl, abort(); } return m_Phi.at     (index); }
  Int_t    GetVaneID  ( Int_t index ){ if( index>=m_VaneID.size  () ){ std::cerr << "[ABORT] Wrong index for VaneID"   << std::endl, abort(); } return m_VaneID.at  (index); }
  Double_t GetgT      ( Int_t index ){ if( index>=m_gT.size      () ){ std::cerr << "[ABORT] Wrong index for gT"       << std::endl, abort(); } return m_gT.at      (index); }
  Double_t GetpT      ( Int_t index ){ if( index>=m_pT.size      () ){ std::cerr << "[ABORT] Wrong index for pT"       << std::endl, abort(); } return m_pT.at      (index); }
  Int_t    GetpID     ( Int_t index ){ if( index>=m_pID.size     () ){ std::cerr << "[ABORT] Wrong index for pID"      << std::endl, abort(); } return m_pID.at     (index); }
  Double_t GetEachDepE( Int_t index ){ if( index>=m_EachDepE.size() ){ std::cerr << "[ABORT] Wrong index for EachDepE" << std::endl, abort(); } return m_EachDepE.at(index); }
  Int_t GetClose_Hough_phiz( Int_t index ){ if( index>=m_close_hough_phiz.size() ){ std::cerr << "[ABORT] Wrong index for close-hough-line(phiz)" << std::endl, abort(); } return m_close_hough_phiz.at(index); }
  Int_t GetClusterNo  ( Int_t index ){ if( index>=m_clusterNo.size()){ std::cerr << "[ABORT] Wrong index for clusterNo" << std::endl, abort(); } return m_clusterNo.at(index); }


  
  TH1D*    GetHist_Hough_phiz_Residual( Int_t index ){ if( index>=m_hist_hough_phiz_resi.size() ){ std::cerr << "[ABORT] Wrong index for slope of hough-line(phiz)" << std::endl, abort(); } return m_hist_hough_phiz_resi.at(index); }
  TH1D*    GetHist_Hough_phiz_Slope( Int_t index ){ if( index>=m_hist_hough_phiz_slope.size() ){ std::cerr << "[ABORT] Wrong index for offset of hough-line(phiz)" << std::endl, abort(); } return m_hist_hough_phiz_slope.at(index); }
  TH1D*    GetHist_Hough_phiz_Offset( Int_t index ){ if( index>=m_hist_hough_phiz_offset.size() ){ std::cerr << "[ABORT] Wrong index for slope-offset of hough-line(phiz)" << std::endl, abort(); } return m_hist_hough_phiz_offset.at(index); }
  TH2D*    GetHist_Hough_phiz_Slope_Offset( Int_t index ){ if( index>=m_hist_hough_phiz_slope_offset.size() ){ std::cerr << "[ABORT] Wrong index for residual of hough-line(phiz)" << std::endl, abort(); } return m_hist_hough_phiz_slope_offset.at(index); }
  TH2D*    GetHist_Hough_phiz    ( Int_t index ){ if( index>=m_hist_hough_phiz.size     () ){ std::cerr << "[ABORT] Wrong index for hough-plane(phiz)"            << std::endl, abort(); } return m_hist_hough_phiz.at     (index); }
  TF1*     GetFunc_Hough_phiz    ( Int_t index ){ if( index>=m_func_hough_phiz.size     () ){ std::cerr << "[ABORT] Wrong index for hough-fit line(phiz)"         << std::endl, abort(); } return m_func_hough_phiz.at     (index); }
  Double_t GetPar0_Hough_phiz( Int_t index ){ if( index>=m_hough_phiz_par0.size() ){ std::cerr << "[ABORT] Wrong index for hough-line parameter0(phiz)" << std::endl, abort(); } return m_hough_phiz_par0.at(index); }
  Double_t GetPar1_Hough_phiz( Int_t index ){ if( index>=m_hough_phiz_par1.size() ){ std::cerr << "[ABORT] Wrong index for hough-line parameter1(phiz)" << std::endl, abort(); } return m_hough_phiz_par1.at(index); }

  private:

  const Int_t    m_geom_nvane;
  const Double_t m_geom_r_inner;
  const Double_t m_geom_r_outer;
  const Double_t m_geom_z_min;
  const Double_t m_geom_z_max;
  Double_t*      m_geom_phi;
  const Int_t    m_hough_nstep_theta;
  // Hits information
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
  std::vector<Int_t>    m_fl_hough_phiz;    // flag if this point is Hough-transformded or not.
  std::vector<Int_t>    m_close_hough_phiz; // indicate which hough-line is clonsed (initial value:-999) 
  // Sort
  std::vector<Int_t>    m_order_VaneID;
  std::vector<Int_t>    m_order_gT;

  // Hough Fit
  std::vector<std::vector<Double_t> > m_hough_phiz_rho;
  std::vector<Double_t>               m_hough_phiz_theta;
  
  std::vector<Double_t> m_hough_phiz_par0;
  std::vector<Double_t> m_hough_phiz_par1;
  std::vector<Double_t> m_hough_phiz_rho_max;
  std::vector<Double_t> m_hough_phiz_theta_max;

 
  std::vector<TH2D*> m_hist_hough_phiz;
  std::vector<TH2D*> m_hist_hough_phiz_slope_offset;
  std::vector<TH1D*> m_hist_hough_phiz_slope;
  std::vector<TH1D*> m_hist_hough_phiz_offset;
  std::vector<TF1*>  m_func_hough_phiz;
  std::vector<TH1D*> m_hist_hough_phiz_resi;

  // Clustering
  std::vector<Int_t> m_clusterNo;
  

 public:
  // Extrapolation
  void  CalcPerpLineSeg( Double_t x1,  Double_t y1,  Double_t x2, Double_t y2, Double_t& slope, Double_t& offset );
  void  CircleBy3Point ( Double_t x1,  Double_t y1,  Double_t x2, Double_t y2, Double_t x3, Double_t y3,
			 Double_t& x0, Double_t& y0, Double_t& r );
  Int_t IntersectionCircleLine( Double_t x0, Double_t y0, Double_t r, Double_t slope, Double_t offset, Double_t& x1, Double_t& y1, Double_t& x2, Double_t& y2 );
  Int_t Extrapolation              ( Bool_t fl_clockwise, Int_t index1, Int_t index2, Int_t index3, Int_t VaneID, Double_t& extrap_r, Double_t& extrap_z );
  Int_t Extrapolation_Clockwise    ( Int_t index1, Int_t index2, Int_t index3, Int_t VaneID, Double_t& extrap_r, Double_t& extrap_z );
  Int_t Extrapolation_AntiClockwise( Int_t index1, Int_t index2, Int_t index3, Int_t VaneID, Double_t& extrap_r, Double_t& extrap_z );

  void Test();
};

#endif
