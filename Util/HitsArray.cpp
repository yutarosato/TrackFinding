#include "HitsArray.h"

HitsArray::HitsArray():
  m_geom_nvane       (  48),
  m_geom_r_inner     (  20),
  m_geom_r_outer     (  80),
  m_geom_z_min       (-800),
  m_geom_z_max       ( 800),
  m_hough_nstep_theta( 180)
{
  m_geom_phi = new Double_t[m_geom_nvane];
  for( Int_t ivane=0; ivane<m_geom_nvane; ivane++ ){
    m_geom_phi[ivane] = 2.0*TMath::Pi()/(m_geom_nvane/2.0)*((Int_t)((ivane+1)/2));
    if( ivane%2 ) m_geom_phi[ivane] -= 0.01; // tmppppp
    else          m_geom_phi[ivane] += 0.01; // tmppppp
    //std::cout << std::setw(3) << std::right << ivane << " : "
    //<< std::setw(10) << std::right << m_geom_phi[ivane] << std::endl;
  }
}

  HitsArray::~HitsArray(){
  delete m_geom_phi;
};

void HitsArray::InputHits( Int_t index, Double_t x, Double_t y, Double_t z, Double_t pt, Double_t gt, Int_t pID, Double_t EachDepE ){
  m_Index.push_back           ( index                   );
  m_X.push_back               ( x                       );
  m_Y.push_back               ( y                       );
  m_Z.push_back               ( z                       );
  m_R.push_back               ( sqrt(pow(x,2)+pow(y,2)) );
  m_Phi.push_back             ( Phi_uk(y,x)             );
  m_VaneID.push_back          ( GetVaneID(Phi_uk(y,x))  );
  m_gT.push_back              ( gt                      );
  m_pT.push_back              ( pt                      );
  m_pID.push_back             ( pID                     );
  m_EachDepE.push_back        ( EachDepE                );
  m_fl_hough_phiz.push_back   ( 1                       );
  m_close_hough_phiz.push_back( -999                    );
}

void HitsArray::ClearEvent(){
  m_Index.clear();
  m_X.clear();
  m_Y.clear();
  m_Z.clear();
  m_R.clear();
  m_Phi.clear();
  m_VaneID.clear();
  m_gT.clear();
  m_pT.clear();
  m_pID.clear();
  m_EachDepE.clear();
  m_fl_hough_phiz.clear();
  m_close_hough_phiz.clear();

  m_order_VaneID.clear();
  m_order_gT.clear();
  
  m_hough_phiz_rho.clear();
  m_hough_phiz_theta.clear();

  m_hough_phiz_par0.clear();
  m_hough_phiz_par1.clear();

  for( Int_t ivec=0; ivec<m_hist_hough_phiz.size     (); ivec++ ) delete m_hist_hough_phiz.at     (ivec);
  for( Int_t ivec=0; ivec<m_hist_hough_phiz_resi.size(); ivec++ ) delete m_hist_hough_phiz_resi.at(ivec);
  for( Int_t ivec=0; ivec<m_func_hough_phiz.size     (); ivec++ ) delete m_func_hough_phiz.at     (ivec);
  m_hist_hough_phiz.clear     ();
  m_hist_hough_phiz_resi.clear();
  m_func_hough_phiz.clear     ();  



  return;
}

void HitsArray::Print( Int_t fl_message ){
  if( fl_message < 2 ) return;
  for( Int_t ivec=0; ivec<m_X.size(); ivec++ ){
    std::cout << "              "
	      << std::setw(3) << std::right << m_Index.at(ivec) << " : pID = "
	      << m_pID.at(ivec)   << ", (x,y,z) = ("
	      << std::setw(7) << std::right << Form("%.2f",m_X.at       (ivec)) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Y.at       (ivec)) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Z.at       (ivec)) << "), phi = "
	      << std::setw(5) << std::right << Form("%.2f",m_Phi.at     (ivec)) << ", vane-ID = "
	      << std::setw(2) << std::right << Form("%d",  m_VaneID.at  (ivec)) << ", t(proper) = "
	      << std::setw(8) << std::right << Form("%.4f",m_pT.at      (ivec)) << ", t(global) = "
	      << std::setw(8) << std::right << Form("%.3f",m_gT.at      (ivec)) << ", Edep = "
	      << std::setw(5) << std::right << Form("%.5f",m_EachDepE.at(ivec)) << " : "
      	      << std::setw(2) << std::right << m_fl_hough_phiz.at   (ivec)    << ", "
      	      << std::setw(4) << std::right << m_close_hough_phiz.at(ivec)
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

void HitsArray::HoughTransform_phiz(){
  TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",m_hist_hough_phiz.size()), "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 500, -250, 250 );

  for( Int_t ivec=0; ivec<m_X.size(); ivec++ ){
    std::vector<Double_t> tmp_vector;
    if( m_Z[ivec]<-100 || m_Z[ivec]> 100 ){ // tmppppppp
      m_fl_hough_phiz[ivec] = 0;
      continue;
    }

    for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ){
      Double_t rho = (m_Phi[ivec]*180/TMath::Pi()*TMath::Cos((double)istep*TMath::Pi()/180)
		      +m_Z[ivec]*TMath::Sin((double)istep*TMath::Pi()/180));
      Double_t theta = istep;
      tmp_vector.push_back( rho );
      hist_hough_phiz->Fill( theta, rho );
    }
    m_hough_phiz_rho.push_back( tmp_vector );
  }

  m_hist_hough_phiz.push_back( hist_hough_phiz );
  for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ) m_hough_phiz_theta.push_back(istep);

  return;
}
void HitsArray::HoughFit_phiz(){
  // Search Line
  Int_t max_xbin, max_ybin, max_zbin;
  while(1){
    m_hist_hough_phiz.at(0)->GetMaximumBin(max_xbin,max_ybin, max_zbin);
    Double_t rho   = m_hist_hough_phiz.at(0)->GetYaxis()->GetBinLowEdge(max_ybin) + m_hist_hough_phiz.at(0)->GetXaxis()->GetBinWidth(max_ybin)/2.0;
    Double_t theta = m_hist_hough_phiz.at(0)->GetXaxis()->GetBinLowEdge(max_xbin) + m_hist_hough_phiz.at(0)->GetXaxis()->GetBinWidth(max_xbin)/2.0;
    
    Double_t par0 =  rho/TMath::Sin(theta*TMath::Pi()/180.0);
    Double_t par1 = -1.0/TMath::Tan(theta*TMath::Pi()/180.0)*180.0/TMath::Pi();
    m_hough_phiz_par0.push_back( par0 );
    m_hough_phiz_par1.push_back( par1 );
    TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",m_hist_hough_phiz.size()), "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 500, -250, 250 );
    for( Int_t ivec=0; ivec<m_hough_phiz_rho.size(); ivec++ ){
      if( TMath::Abs(m_hough_phiz_rho[ivec][max_xbin-1]==m_hist_hough_phiz.at(0)->GetYaxis()->GetBinLowEdge(max_ybin)) ) continue;
      for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ) hist_hough_phiz->Fill( m_hough_phiz_theta.at(ivec), m_hough_phiz_rho.at(ivec).at(istep) );
    }
    m_hist_hough_phiz.push_back( hist_hough_phiz );
    break;
  }


  // Make Objects(TF1)
  for( Int_t iline=0; iline<m_hough_phiz_par0.size(); iline++ ){
    TF1* func_hough_phiz = new TF1( Form("func_hough_phiz_%d",iline), "[0]+[1]*x", 0.0, 2*TMath::Pi() );
    func_hough_phiz->SetLineColor( iline+1 );
    func_hough_phiz->SetParameter( 0, m_hough_phiz_par0.at(iline) );
    func_hough_phiz->SetParameter( 1, m_hough_phiz_par1.at(iline) );
    m_func_hough_phiz.push_back( func_hough_phiz );
  }

  return;
}

void HitsArray::CalcHoughResidual_phiz(){
  for( Int_t iline=m_hough_phiz_par0.size()-1; iline>=0; iline-- ){
    TH1D* hist_resi_phiz = new TH1D( Form("hist_resi_phiz_%d",iline),  Form("Residual_%d(#phi-Z);",iline),  100, -20, 20 );
    Double_t par0 = m_hough_phiz_par0.at(iline);
    Double_t par1 = m_hough_phiz_par1.at(iline);
    hist_resi_phiz->SetLineColor( m_hist_hough_phiz_resi.size()+1 );
    Double_t fl_angle = TMath::ATan(par1*180.0/TMath::Pi())*180.0/TMath::Pi();
    for( Int_t ivec=0; ivec<m_X.size(); ivec++ ){
      Double_t residual;
      if( fl_angle>89.9 ) residual = m_Z.at(ivec) - ( par0+par1*m_Phi.at(ivec) );
      else                residual = m_Phi.at(ivec) - ((m_Z.at(ivec) - par0)/par1) *180.0/TMath::Pi();
      hist_resi_phiz->Fill( residual );
      if( TMath::Abs(residual) < 10 ) m_close_hough_phiz[ivec] = iline;
    }
    m_hist_hough_phiz_resi.push_back( hist_resi_phiz );
  }
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t HitsArray::Phi_uk(Double_t y, Double_t x){ // 0~phi~2pi : this definition is used in Ueno-san's codes
  Double_t phi = TMath::ATan2(y,x);
  if( phi<0 ) phi += 2*TMath::Pi();
  return phi;
}

Int_t HitsArray::GetVaneID( Double_t f_phi ){
  ///*
  Int_t vaneID = (Int_t)(f_phi/(2.0*TMath::Pi()/m_geom_nvane));
  if( vaneID==m_geom_nvane ) vaneID = 0;
  return vaneID;
  //*/
  /* uk codes
  Double_t tmpVane=(f_phi/(2*TMath::Pi()/m_geom_nvane));
  Int_t VNum = (Int_t)tmpVane;
  if(tmpVane-VNum >=0.5){
    VNum++;
    if(VNum==m_geom_nvane) VNum=0;
  }
  return VNum;
  */
}
