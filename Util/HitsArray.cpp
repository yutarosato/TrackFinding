#include "HitsArray.h"

HitsArray::HitsArray():
  m_geom_nvane       (  48),
  m_geom_r_inner     (  64),
  m_geom_r_outer     ( 264),
  m_geom_z_min       (-200),
  m_geom_z_max       ( 200),
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
  ClearEvent();
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
  m_clusterNo.push_back       ( -999                    );
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
  m_hough_phiz_rho_max.clear();
  m_hough_phiz_theta_max.clear();

  for( Int_t ivec=0; ivec<m_hist_hough_phiz.size       (); ivec++ ) delete m_hist_hough_phiz.at       (ivec);
  for( Int_t ivec=0; ivec<m_hist_hough_phiz_slope_offset.size (); ivec++ ) delete m_hist_hough_phiz_slope_offset.at (ivec);
  for( Int_t ivec=0; ivec<m_hist_hough_phiz_slope.size (); ivec++ ) delete m_hist_hough_phiz_slope.at (ivec);
  for( Int_t ivec=0; ivec<m_hist_hough_phiz_offset.size(); ivec++ ) delete m_hist_hough_phiz_offset.at(ivec);
  for( Int_t ivec=0; ivec<m_func_hough_phiz.size       (); ivec++ ) delete m_func_hough_phiz.at       (ivec);
  for( Int_t ivec=0; ivec<m_hist_hough_phiz_resi.size  (); ivec++ ) delete m_hist_hough_phiz_resi.at  (ivec);
  m_hist_hough_phiz_slope_offset.clear ();
  m_hist_hough_phiz_slope.clear ();
  m_hist_hough_phiz_offset.clear();
  m_hist_hough_phiz.clear       ();
  m_func_hough_phiz.clear       ();  
  m_hist_hough_phiz_resi.clear  ();

  m_clusterNo.clear();
  
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
	      << std::setw(7) << std::right << Form("%.2f",m_Z.at       (ivec)) << "), (R,phi) = ("
	      << std::setw(5) << std::right << Form("%.2f",m_R.at       (ivec)) << ", "
	      << std::setw(5) << std::right << Form("%.2f",m_Phi.at     (ivec)) << "), vane-ID = "
	      << std::setw(2) << std::right << Form("%d",  m_VaneID.at  (ivec)) << ", t(proper) = "
	      << std::setw(8) << std::right << Form("%.4f",m_pT.at      (ivec)) << ", t(global) = "
	      << std::setw(8) << std::right << Form("%.3f",m_gT.at      (ivec)) << ", Edep = "
	      << std::setw(5) << std::right << Form("%.5f",m_EachDepE.at(ivec)) << " : "
      	      << std::setw(2) << std::right << m_fl_hough_phiz.at   (ivec)    << ", "
      	      << std::setw(4) << std::right << m_close_hough_phiz.at(ivec)    << ", "
      	      << std::setw(4) << std::right << m_clusterNo.at       (ivec)    << ", "
	      << std::endl;
  }
  return;
}

void HitsArray::Print_VaneID_Order( Int_t fl_message ){
  if( fl_message < 2 ) return;
  for( Int_t ivane=0; ivane<m_order_VaneID.size(); ivane++ ){
    std::cout << "              "
	      << std::setw(3) << std::right << m_Index.at(m_order_VaneID.at(ivane)) << " : pID = "
	      << m_pID.at(m_order_VaneID.at(ivane))   << ", (x,y,z) = ("
	      << std::setw(7) << std::right << Form("%.2f",m_X.at       (m_order_VaneID.at(ivane))) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Y.at       (m_order_VaneID.at(ivane))) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Z.at       (m_order_VaneID.at(ivane))) << "), (R,phi) = ("
	      << std::setw(5) << std::right << Form("%.2f",m_R.at       (m_order_VaneID.at(ivane))) << ", "
	      << std::setw(5) << std::right << Form("%.2f",m_Phi.at     (m_order_VaneID.at(ivane))) << "), vane-ID = "
	      << std::setw(2) << std::right << Form("%d",  m_VaneID.at  (m_order_VaneID.at(ivane))) << ", t(proper) = "
	      << std::setw(8) << std::right << Form("%.4f",m_pT.at      (m_order_VaneID.at(ivane))) << ", t(global) = "
	      << std::setw(8) << std::right << Form("%.3f",m_gT.at      (m_order_VaneID.at(ivane))) << ", Edep = "
	      << std::setw(5) << std::right << Form("%.5f",m_EachDepE.at(m_order_VaneID.at(ivane))) << " : "
      	      << std::setw(2) << std::right << m_fl_hough_phiz.at   (m_order_VaneID.at(ivane))    << ", "
      	      << std::setw(4) << std::right << m_close_hough_phiz.at(m_order_VaneID.at(ivane))    << ", "
      	      << std::setw(4) << std::right << m_clusterNo.at       (m_order_VaneID.at(ivane))    << ", "
	      << std::endl;
  }
  return;
}

void HitsArray::Print_gT_Order( Int_t fl_message ){
  if( fl_message < 2 ) return;
  for( Int_t ivane=0; ivane<m_order_gT.size(); ivane++ ){
    std::cout << "              "
	      << std::setw(3) << std::right << m_Index.at(m_order_gT.at(ivane)) << " : pID = "
	      << m_pID.at(m_order_gT.at(ivane))   << ", (x,y,z) = ("
	      << std::setw(7) << std::right << Form("%.2f",m_X.at       (m_order_gT.at(ivane))) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Y.at       (m_order_gT.at(ivane))) << ", "
	      << std::setw(7) << std::right << Form("%.2f",m_Z.at       (m_order_gT.at(ivane))) << "), (R,phi) = ("
	      << std::setw(5) << std::right << Form("%.2f",m_R.at       (m_order_gT.at(ivane))) << ", "
	      << std::setw(5) << std::right << Form("%.2f",m_Phi.at     (m_order_gT.at(ivane))) << "), vane-ID = "
	      << std::setw(2) << std::right << Form("%d",  m_VaneID.at  (m_order_gT.at(ivane))) << ", t(proper) = "
	      << std::setw(8) << std::right << Form("%.4f",m_pT.at      (m_order_gT.at(ivane))) << ", t(global) = "
	      << std::setw(8) << std::right << Form("%.3f",m_gT.at      (m_order_gT.at(ivane))) << ", Edep = "
	      << std::setw(5) << std::right << Form("%.5f",m_EachDepE.at(m_order_gT.at(ivane))) << " : "
      	      << std::setw(2) << std::right << m_fl_hough_phiz.at   (m_order_gT.at(ivane))    << ", "
      	      << std::setw(4) << std::right << m_close_hough_phiz.at(m_order_gT.at(ivane))    << ", "
      	      << std::setw(4) << std::right << m_clusterNo.at       (m_order_gT.at(ivane))    << ", "
	      << std::endl;
  }
  return;
}

void HitsArray::CalcOrder(){
  if( m_order_VaneID.size() ) m_order_VaneID.clear();
  if( m_order_gT.size()     ) m_order_gT.clear();
  
  std::multimap<Int_t,Double_t> tMap_VaneID;
  std::multimap<Int_t,Double_t> tMap_gT;
  for( Int_t ivec=0; ivec< m_Index.size(); ivec++ ){
    tMap_VaneID.insert( std::make_pair(m_VaneID.at(ivec),m_Index.at(ivec)) );
    tMap_gT.insert    ( std::make_pair(m_gT.at    (ivec),m_Index.at(ivec)) );
  }
  std::multimap<Int_t,Double_t>::iterator it_VaneID = tMap_VaneID.begin();
  std::multimap<Int_t,Double_t>::iterator it_gT     = tMap_gT.begin    ();
  while( it_VaneID != tMap_VaneID.end() ){ m_order_VaneID.push_back( (*it_VaneID).second ); it_VaneID++; }
  while( it_gT     != tMap_gT.end    () ){ m_order_gT.push_back    ( (*it_gT    ).second ); it_gT++;     }
}

void HitsArray::HoughTransform_phiz(){
  TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",m_hist_hough_phiz.size()), "Hough(#phi-Z #rightarrow #theta-#rho);#theta [#circ];#rho [mm]", 180, 0, 180, 1000, -500, 500 );

  for( Int_t ivec=0; ivec<m_X.size(); ivec++ ){
    std::vector<Double_t> tmp_vector;
    //if( m_Z.at(ivec)<-100 || m_Z.at(ivec)> 100 ){ // tmppppppp
    //m_fl_hough_phiz.at(ivec) = 0;
    //continue;
    //}

    for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ){
      Double_t rho = (m_Phi.at(ivec)*180/TMath::Pi()*TMath::Cos((double)istep*TMath::Pi()/180)
		      +m_Z.at(ivec)*TMath::Sin((double)istep*TMath::Pi()/180));
      Double_t theta = istep;
      tmp_vector.push_back( rho );
      if( theta < 30.0 || theta > 150 ) continue; // tmpppppp
      hist_hough_phiz->Fill( theta+1.0e-5, rho );
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
    TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",m_hist_hough_phiz.size()), "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 1000, -500, 500 );
    for( Int_t ivec=0; ivec<m_hough_phiz_rho.size(); ivec++ ){
      if( TMath::Abs(m_hough_phiz_rho.at(ivec).at(max_xbin-1)==m_hist_hough_phiz.at(0)->GetYaxis()->GetBinLowEdge(max_ybin)) ) continue;
      for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ){
	hist_hough_phiz->Fill( m_hough_phiz_theta.at(istep), m_hough_phiz_rho.at(ivec).at(istep) );
      }
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
  
  /*
void HitsArray::HoughFit_phiz(){
  // Search Line
  // setting parameter
  //const Int_t range_max_theta = 150;
  //const Int_t range_min_theta =  30;
  
  Int_t max_xbin, max_ybin, max_zbin;
  Int_t cnt = 0;
  while(1){
    TH1D* hist_hough_phiz_slope  = new TH1D( Form("hist_hough_phiz_slope_%d", m_hist_hough_phiz_slope.size ()), "Slope at hough-point(#phi-Z);Slope",   100, -10, 10 );
    TH1D* hist_hough_phiz_offset = new TH1D( Form("hist_hough_phiz_offset_%d",m_hist_hough_phiz_offset.size()), "Offset at hough-point(#phi-Z);Offset", 100, -100, 100 );
    TH2D* hist_hough_phiz_slope_offset = new TH2D( Form("hist_hough_phiz_slope_offset_%d",m_hist_hough_phiz_slope_offset.size()), "Slope-Offset at hough-point(#phi-Z);Slope;Offset", 100, -10, 10, 100, -100, 100 );
    m_hist_hough_phiz.at(cnt)->GetMaximumBin(max_xbin,max_ybin,max_zbin);
    //std::cout << "max_xbin = " << max_xbin << ", "
    //<< "max_ybin = " << max_ybin << ", "
    //<< "max_zbin = " << max_zbin << std::endl;
    //std::cout << "Maximum valu = " << m_hist_hough_phiz.at(cnt)->GetBinContent(max_xbin,max_ybin) << std::endl;
    Double_t rho   = m_hist_hough_phiz.at(cnt)->GetYaxis()->GetBinLowEdge(max_ybin) + m_hist_hough_phiz.at(cnt)->GetYaxis()->GetBinWidth(max_ybin)/2.0;
    Double_t theta = m_hist_hough_phiz.at(cnt)->GetXaxis()->GetBinLowEdge(max_xbin) + m_hist_hough_phiz.at(cnt)->GetXaxis()->GetBinWidth(max_xbin)/2.0;
    //std::cout << "rho = " << rho << ", theta = " << theta << std::endl;    
    Double_t par0 =  rho/TMath::Sin(theta*TMath::Pi()/180.0);
    Double_t par1 = -1.0/TMath::Tan(theta*TMath::Pi()/180.0)*180.0/TMath::Pi();
    //std::cout << "par0 = " << par0 << ", par1 = " << par1 << std::endl;    

    for( Int_t ivec=0; ivec<m_hough_phiz_rho.size(); ivec++ ){
      Double_t offset = m_hough_phiz_rho.at(ivec).at(max_xbin-1) - rho;
      if( TMath::Abs(offset) >100 ) continue;
      Double_t del_rho   = m_hough_phiz_rho.at(ivec).at(max_xbin) - m_hough_phiz_rho.at(ivec).at(max_xbin-1);
      Double_t del_theta = m_hough_phiz_theta.at(max_xbin)        - m_hough_phiz_theta.at(max_xbin-1);
      Double_t slope = del_rho/del_theta;
      std::cout << "    ivec = " << ivec << "  slope = " << slope << ", offset = " << offset << std::endl;
      hist_hough_phiz_slope ->Fill(slope);
      hist_hough_phiz_offset->Fill(offset);
      hist_hough_phiz_slope_offset->Fill(slope,offset);
    }

    m_hist_hough_phiz_slope.push_back( hist_hough_phiz_slope );
    m_hist_hough_phiz_offset.push_back( hist_hough_phiz_offset );
    m_hist_hough_phiz_slope_offset.push_back( hist_hough_phiz_slope_offset );
    
    m_hough_phiz_par0.push_back( par0 );
    m_hough_phiz_par1.push_back( par1 );
    TH2D* hist_hough_phiz = new TH2D( Form("hist_hough_phiz_%d",m_hist_hough_phiz.size()), "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 1000, -500, 500 );
    for( Int_t ivec=0; ivec<m_hough_phiz_rho.size(); ivec++ ){
      if( TMath::Abs(m_hough_phiz_rho.at(ivec).at(max_xbin-1)==m_hist_hough_phiz.at(0)->GetYaxis()->GetBinLowEdge(max_ybin)) ) continue; // removed points which is already used.
      for( Int_t istep=0; istep<m_hough_nstep_theta; istep++ ) hist_hough_phiz->Fill( m_hough_phiz_theta.at(istep), m_hough_phiz_rho.at(ivec).at(istep) );
    }
    m_hist_hough_phiz.push_back( hist_hough_phiz );
    break;
    cnt++;
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
  */
void HitsArray::CalcHoughResidual_phiz(){
  for( Int_t iline=m_hough_phiz_par0.size()-1; iline>=0; iline-- ){
    TH1D* hist_resi_phiz = new TH1D( Form("hist_resi_phiz_%d",iline),  Form("Residual_%d(#phi-Z);",iline),  100, -20, 20 );
    Double_t par0 = m_hough_phiz_par0.at(iline);
    Double_t par1 = m_hough_phiz_par1.at(iline);
    //std::cout << "   par0 = " << par0 << ", par1 = " << par1 << std::endl;
    hist_resi_phiz->SetLineColor( m_hist_hough_phiz_resi.size()+1 );
    for( Int_t ivec=0; ivec<m_X.size(); ivec++ ){
      Double_t residual = m_Z.at(ivec) - ( par0+par1*m_Phi.at(ivec) );
      hist_resi_phiz->Fill( residual );
      if( TMath::Abs(residual) < 10 ) m_close_hough_phiz.at(ivec) = iline;
      //std::cout << "               ivec = " << ivec << ", residual = " << residual << " : Phi = " << m_Phi.at(ivec) << ", Z = " << m_Z.at(ivec) << std::endl;
    }
    m_hist_hough_phiz_resi.push_back( hist_resi_phiz );
  }
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HitsArray::Clustering( Int_t fl_message ){
  if( fl_message > 1 ) std::cout << "Clustering Start" << std::endl
				 << "Nline = " << m_hough_phiz_par0.size() << std::endl;

  for( Int_t iline=0; iline<m_hough_phiz_par0.size(); iline++ ){ // START LINE-LOOP
    if( fl_message > 1 ) std::cout << "   iline = " << iline << ", Vane-ID : ";

    std::vector<Int_t> cluster_index; // index for Vane-ID order
    for( Int_t ivane=0; ivane<m_order_VaneID.size(); ivane++ ){
      if( m_close_hough_phiz.at(m_order_VaneID.at(ivane))==iline ){
	m_clusterNo.at( m_order_VaneID.at(ivane) ) = iline;
	cluster_index.push_back(ivane);
	if( fl_message > 1 ) std::cout << m_VaneID.at( m_order_VaneID.at(ivane) ) << ", ";
      }
    }
    if( fl_message > 1 ) std::cout << "seed sluster : " << cluster_index.size() << std::endl;
    if( cluster_index.size()<4 ) return;

    // forward ++++++++++++++++++++++++++;
    Int_t cnt_miss_forward = 0;
    Double_t pre_target_VaneID = -999;
    Double_t x0;
    Double_t y0;
    Double_t r;
    Double_t dphi;
    while( cnt_miss_forward < 5 ){
      Double_t extrap_r;
      Double_t extrap_z;
      Int_t index1 = m_order_VaneID.at(cluster_index.at(cluster_index.size()-3));
      Int_t index2 = m_order_VaneID.at(cluster_index.at(cluster_index.size()-2));
      Int_t index3 = m_order_VaneID.at(cluster_index.at(cluster_index.size()-1));
      Int_t current_VaneID = m_VaneID.at(index3);
      Int_t target_VaneID  = (pre_target_VaneID < 0 ? m_VaneID.at(index3)+1 : pre_target_VaneID+1 );
      if( fl_message > 1 ) std::cout << std::endl
				     << "Extrapolate(forward) from "
				     << m_VaneID.at( index1 ) << "&"
				     << m_VaneID.at( index2 ) << "&"
				     << m_VaneID.at( index3 ) << std::endl;
      
      // Find next hit point by extrapolation
      while(1){
	if( target_VaneID==m_geom_nvane ) target_VaneID = 0;
	if( current_VaneID==target_VaneID ){
	  target_VaneID = -999;
	  break;
	}else if( Extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi )>0 ){
	  break;
	}else{
	  target_VaneID++;
	}
      }
      if( target_VaneID < 0 ) break;
      pre_target_VaneID = target_VaneID;
      if( fl_message > 1 ) std::cout << "     target VaneID = "
				     << target_VaneID << " : R = "
				     << extrap_r      << ", Z = "
				     << extrap_z      << " : x0 = "
				     << x0            << ", y0 = "
				     << y0            << ", r = "
				     << r             << ", dphi = "
				     << dphi          << std::endl;
      
      // compare actual hits with extrapolated point.
      Double_t dev_min   = 10000;
      Double_t index_min = -999;
      for( Int_t ivane=0; ivane<m_order_VaneID.size(); ivane++ ){
	if( m_VaneID.at(m_order_VaneID.at(ivane))!=target_VaneID ) continue;
	if( m_clusterNo.at(m_order_VaneID.at(ivane))>=0          ) continue;
	Double_t dev_r = TMath::Abs( extrap_r - m_R.at(m_order_VaneID.at(ivane)) );
	Double_t dev_z = TMath::Abs( extrap_z - m_Z.at(m_order_VaneID.at(ivane)) );
	if( fl_message > 1 ) std::cout << "       dev(r) = " << dev_r << ", dev(z) = " << dev_z << ", dev_min = " << sqrt(pow(dev_r,2)+pow(dev_z,2)) << std::endl;
	if( dev_min > sqrt(pow(dev_r,2)+pow(dev_z,2)) ){
	  dev_min   = sqrt(pow(dev_r,2)+pow(dev_z,2));
	  index_min = ivane;
	}
      }
      if( dev_min < 5.0*dphi || dev_min < 10 ){ // tmpppp
	m_clusterNo.at( m_order_VaneID.at(index_min) ) = iline;
	cluster_index.push_back(index_min);
	cnt_miss_forward = 0;
	if( fl_message > 1 ) std::cout << "       => added the hit into the cluster" << std::endl;
      }else{
	if( !(extrap_r < m_geom_r_inner || extrap_r > m_geom_r_outer) ) cnt_miss_forward++;
	if( fl_message > 1 ) std::cout << "       => can not find hit points by extrapolation : dev_min = " << dev_min << ", cnt_miss_foward = " << cnt_miss_forward << std::endl;
      }
    }
    // backward ++++++++++++++++++++++++++;
    Int_t cnt_miss_backward = 0;
    pre_target_VaneID = -999;
    while( cnt_miss_backward < 5 ){
      Double_t extrap_r;
      Double_t extrap_z;
      Int_t index1 = m_order_VaneID.at(cluster_index.at(2));
      Int_t index2 = m_order_VaneID.at(cluster_index.at(1));
      Int_t index3 = m_order_VaneID.at(cluster_index.at(0));
      Int_t current_VaneID = m_VaneID.at(index3);
      Int_t target_VaneID  = (pre_target_VaneID < 0 ? m_VaneID.at(index3)-1 : pre_target_VaneID-1 );
      if( fl_message > 1 ) std::cout << std::endl
				     << "Extrapolate(backward) from "
				     << m_VaneID.at( index1 ) << "&"
				     << m_VaneID.at( index2 ) << "&"
				     << m_VaneID.at( index3 ) << std::endl;
      // Find next hit point by extrapolation
      while(1){
	if( target_VaneID==-1 ) target_VaneID = m_geom_nvane-1;
	if( current_VaneID==target_VaneID                                                        ){
	  target_VaneID = -999;
	  break;
	}else if( Extrapolation( index1, index2, index3, target_VaneID, extrap_r, extrap_z, x0, y0, r, dphi )>0 ){
	  break;
	}else{
	  target_VaneID--;
	}
      }
      if( target_VaneID < 0 ) break;
      pre_target_VaneID = target_VaneID;
      if( fl_message > 1 ) std::cout << "     target VaneID = "
				     << target_VaneID << " : R = "
				     << extrap_r      << ", Z = "
				     << extrap_z      << " : x0 = "
				     << x0            << ", y0 = "
				     << y0            << ", r = "
				     << r             << ", dphi = "
				     << dphi          << std::endl;

      // compare actual hits with extrapolated point.
      Double_t dev_min   = 10000;
      Double_t index_min = -999;
      for( Int_t ivane=0; ivane<m_order_VaneID.size(); ivane++ ){
	if( m_VaneID.at(m_order_VaneID.at(ivane))!=target_VaneID ) continue;
	if( m_clusterNo.at(m_order_VaneID.at(ivane))>=0          ) continue;
	Double_t dev_r = TMath::Abs( extrap_r - m_R.at(m_order_VaneID.at(ivane)) );
	Double_t dev_z = TMath::Abs( extrap_z - m_Z.at(m_order_VaneID.at(ivane)) );
	if( fl_message > 1 ) std::cout << "       dev(r) = " << dev_r << ", dev(z) = " << dev_z << ", dev_min = " << sqrt(pow(dev_r,2)+pow(dev_z,2)) << std::endl;
	if( dev_min > sqrt(pow(dev_r,2)+pow(dev_z,2)) ){
	  dev_min   = sqrt(pow(dev_r,2)+pow(dev_z,2));
	  index_min = ivane;
	}
      }
      if( dev_min < 5.0*dphi || dev_min < 10 ){ // tmpppp
	m_clusterNo.at( m_order_VaneID.at(index_min) ) = iline;
	cluster_index.insert(cluster_index.begin(), index_min);
	cnt_miss_backward = 0;
	if( fl_message > 1 ) std::cout << "       => added the hit into the cluster" << std::endl;
      }else{
	if( !(extrap_r < m_geom_r_inner || extrap_r > m_geom_r_outer) ) cnt_miss_backward++;
	if( fl_message > 1 ) std::cout << "can not find hit points by extrapolation : dev_min = " << dev_min << std::endl;
      }
    }
    
  } // END LINE-LOOP
  if( fl_message > 1 ) std::cout << "Clustering finish" << std::endl;
  
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

Int_t HitsArray::CalcPerpLineSeg( Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t& slope, Double_t& offset ){
  if( TMath::Abs(x1-x2)<1.0e-5 && TMath::Abs(y1-y2)<1.0e-5 ){
    //std::cerr << "[ABORT] Same points are input" << std::endl;
    return -1;
  }
  Double_t slope1 = ( x2!=x1 ? (y2-y1)/(x2-x1) : 0.0 );
  if( TMath::Abs(slope1)<1.0e-5 ){
    //std::cerr << "[ABORT] Invalid slope : " << slope1 << std::endl;
    return -1;
  }
  slope = -1.0/slope1;
  offset = (y1+y2)/2.0 - slope*(x1+x2)/2.0;
  return 1;
}

Int_t HitsArray::CircleBy3Point( Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3,
				Double_t& x0, Double_t& y0, Double_t& r ){
  if( TMath::Abs(x1-x2)<1.0e-5 && TMath::Abs(y1-y2)<1.0e-5 ||
      TMath::Abs(x2-x3)<1.0e-5 && TMath::Abs(y2-y3)<1.0e-5 ||
      TMath::Abs(x3-x1)<1.0e-5 && TMath::Abs(y3-y1)<1.0e-5 ){
    //std::cerr << "[ABORT] Same points are input" << std::endl;
    return -1;
  }

  Double_t slope1;
  Double_t slope2;
  Double_t offset1;
  Double_t offset2;

  CalcPerpLineSeg( x1, y1, x2, y2, slope1, offset1 );
  CalcPerpLineSeg( x2, y2, x3, y3, slope2, offset2 );
  if( TMath::Abs(slope2-slope1)<1.0e-5 ){
    //std::cerr << "[ABORT] Three points on straight line" << std::endl;
    return -1;
  }
  x0 = (offset2-offset1)/(slope1-slope2);
  y0 = slope1 * x0 + offset1;
  r = sqrt( pow(x1-x0,2) + pow(y1-y0,2) );
  return 1;
}

Int_t HitsArray::IntersectionCircleLine( Double_t x0, Double_t y0, Double_t r, Double_t slope, Double_t offset, Double_t& x1, Double_t& y1, Double_t& x2, Double_t& y2 ){
  // (x-x0)^2 + (y-y0)^2 = r^2
  // y = slope*x + offset
  // c2*x^2 + c1*x^1 + c0 = 0;
  Double_t c2 = pow(slope,2) + 1.0;
  Double_t c1 = 2.0*(offset-y0)*slope - 2.0* x0;
  Double_t c0 = pow(x0,2) + pow(offset-y0,2) - pow(r,2);
  if( pow(c1,2) - 4*c2*c0 < 0 ) return -1;
  
  x1 = (-c1 + sqrt(pow(c1,2)-4*c2*c0))/(2*c2);
  x2 = (-c1 - sqrt(pow(c1,2)-4*c2*c0))/(2*c2);
  y1 = slope*x1 + offset;
  y2 = slope*x2 + offset;
  
  return 1;

}


Int_t HitsArray::Extrapolation( Int_t index1, Int_t index2, Int_t index3, Int_t VaneID, Double_t& extrap_r, Double_t& extrap_z, Double_t& x0, Double_t& y0, Double_t& r, Double_t& dphi ){
  //if     ( VaneID==m_geom_nvane ) VaneID = 0;
  //else if( VaneID==-1           ) VaneID = m_geom_nvane-1;

  if( CircleBy3Point( m_X.at(index1), m_Y.at(index1), m_X.at(index2), m_Y.at(index2), m_X.at(index3), m_Y.at(index3), x0, y0, r )<0 ) return -1;
  Double_t slope_vane = TMath::Tan(m_geom_phi[VaneID]);

  Double_t phi_org1 = Phi_uk( m_Y.at(index1)-y0, m_X.at(index1)-x0 );
  Double_t phi_org2 = Phi_uk( m_Y.at(index2)-y0, m_X.at(index2)-x0 );
  Double_t phi_org3 = Phi_uk( m_Y.at(index3)-y0, m_X.at(index3)-x0 );


  // judgement of clockwise(phi1 > phi2 > phi3) or anti-clockwise(phi1 < phi2 < phi3)
  Double_t tmp_phi_org1_clockwise = phi_org1;
  Double_t tmp_phi_org2_clockwise = phi_org2;
  Double_t tmp_phi_org3_clockwise = phi_org3;
  while( tmp_phi_org2_clockwise > tmp_phi_org1_clockwise ) tmp_phi_org2_clockwise -= 2.0*TMath::Pi();
  while( tmp_phi_org3_clockwise > tmp_phi_org2_clockwise ) tmp_phi_org3_clockwise -= 2.0*TMath::Pi();
  Double_t dphi_clockwise = TMath::Abs( tmp_phi_org3_clockwise - tmp_phi_org1_clockwise );

  Double_t tmp_phi_org1_anticlockwise = phi_org1;
  Double_t tmp_phi_org2_anticlockwise = phi_org2;
  Double_t tmp_phi_org3_anticlockwise = phi_org3;
  while( tmp_phi_org2_anticlockwise < tmp_phi_org1_anticlockwise ) tmp_phi_org2_anticlockwise += 2.0*TMath::Pi(); 
  while( tmp_phi_org3_anticlockwise < tmp_phi_org2_anticlockwise ) tmp_phi_org3_anticlockwise += 2.0*TMath::Pi();
  Double_t dphi_anticlockwise = TMath::Abs( tmp_phi_org3_anticlockwise - tmp_phi_org1_anticlockwise );

  Bool_t fl_clockwise;
  if( dphi_clockwise < dphi_anticlockwise ){
    fl_clockwise = true;
    phi_org1 = tmp_phi_org1_clockwise;
    phi_org2 = tmp_phi_org2_clockwise;
    phi_org3 = tmp_phi_org3_clockwise;
  }else{
    fl_clockwise = false;
    phi_org1 = tmp_phi_org1_anticlockwise;
    phi_org2 = tmp_phi_org2_anticlockwise;
    phi_org3 = tmp_phi_org3_anticlockwise;
  }

  Double_t x1;
  Double_t y1;
  Double_t x2;
  Double_t y2;
  Int_t    fl_answer = IntersectionCircleLine( x0, y0, r, slope_vane, 0.0, x1, y1, x2, y2 );
  if( fl_answer<0 ) return -1;
  Double_t phi1 = Phi_uk( y1-y0, x1-x0 );
  Double_t phi2 = Phi_uk( y2-y0, x2-x0 );

  //std::cout << "     Extrapolated Circle : x0 = " << x0 << ", y0 = " << y0 << ", r  = " << r << " : fl_clockwise = " << fl_clockwise << std::endl;
  //std::cout << "     x1 = " << m_X.at(index1) << ", y1 = " << m_Y.at(index1) << ", z1 = " << m_Z.at(index1) << ", phi(org1) = " << phi_org1 << " : Vane-ID = " << m_VaneID.at(index1) << std::endl
  //<< "     x2 = " << m_X.at(index2) << ", y2 = " << m_Y.at(index2) << ", z2 = " << m_Z.at(index2) << ", phi(org2) = " << phi_org2 << " : Vane-ID = " << m_VaneID.at(index2) << std::endl
  //<< "     x3 = " << m_X.at(index3) << ", y3 = " << m_Y.at(index3) << ", z3 = " << m_Z.at(index3) << ", phi(org3) = " << phi_org3 << " : Vane-ID = " << m_VaneID.at(index3) << std::endl;
  Double_t extrap_x;
  Double_t extrap_y;
  Double_t extrap_phi;
  if( VaneID/(m_geom_nvane/2)==0 ){
    if     ( y1>0 ){ extrap_x = x1; extrap_y = y1; extrap_phi = phi1; }
    else if( y2>0 ){ extrap_x = x2; extrap_y = y2; extrap_phi = phi2; }
    else           return -1;
  }else if( VaneID/(m_geom_nvane/2)==1 ){
    if     ( y1<0 ){ extrap_x = x1; extrap_y = y1; extrap_phi = phi1; }
    else if( y2<0 ){ extrap_x = x2; extrap_y = y2; extrap_phi = phi2; }
    else           return -1;
  }else{
    return -1;
  }
  if( fl_clockwise ){ // phi3 > extrap_phi
    while( extrap_phi > phi_org3 ) extrap_phi -= 2.0*TMath::Pi();
  }else{ // phi3 < extrap_phi
    while( extrap_phi < phi_org3 ) extrap_phi += 2.0*TMath::Pi();
  }
  dphi = TMath::Abs( phi_org3 - extrap_phi );

  extrap_z = m_Z.at(index3) + (m_Z.at(index3)-m_Z.at(index1))/(phi_org3-phi_org1)*(extrap_phi-phi_org3);  
  extrap_r = sqrt(pow(extrap_x,2) + pow(extrap_y,2));

  //std::cout << "     Extrapx = " << extrap_x << ", extrapy = " << extrap_y << ", extrapz = " << extrap_z << ", extrap_phi = " << extrap_phi << std::endl;
  if( extrap_z < m_geom_z_min   || extrap_z > m_geom_z_max   ) return -1;
  if( extrap_r < m_geom_r_inner || extrap_r > m_geom_r_outer ) return 2; 

  return 1; // 
}

void HitsArray::Test(){
  //Double_t extrap_r;
  //Double_t extrap_z;
  //Extrapolation( 18, 17, 16, 3, extrap_r, extrap_z );
  
  return;
}
