#include "setting.h"
const Int_t    fl_message        = 1;
const Int_t    fl_show           = 200;
const Double_t th_show_energy    = 100.0;
const Int_t    threshold_success = 3; // Hit definition : >= threshold_success/range_success
const Int_t    range_success     = 6;
const Int_t    fl_batch          = 2; // 0(show), 1(batch), 2(batch&save)

// Objects
std::vector<Double_t> v_X;
std::vector<Double_t> v_Y;
std::vector<Double_t> v_Phi;
std::vector<Double_t> v_Z;
std::vector<Int_t>    v_VaneID;
std::vector<Double_t> v_residual_PhiZ;

std::vector<Int_t>    v_closeIndex;
std::vector<Double_t> v_closeX;
std::vector<Double_t> v_closeY;
std::vector<Double_t> v_closePhi;
std::vector<Double_t> v_closeZ;
std::vector<Int_t>    v_closeVaneID;

std::vector<Int_t>    v_clusterIndex;
std::vector<Double_t> v_clusterX;
std::vector<Double_t> v_clusterY;
std::vector<Double_t> v_clusterPhi;
std::vector<Double_t> v_clusterZ;
std::vector<Int_t>    v_clusterVaneID;

std::vector<Double_t> h_X;
std::vector<Double_t> h_Y;
std::vector<Double_t> h_Phi;
std::vector<Double_t> h_Z;

std::vector<Double_t> v_gT; // global time
std::vector<Double_t> v_pT; // proper time

void ResetObject(){
  v_X.clear();
  v_Y.clear();
  v_Phi.clear();
  v_Z.clear();
  v_VaneID.clear();
  v_residual_PhiZ.clear();

  v_closeIndex.clear();
  v_closeX.clear();
  v_closeY.clear();
  v_closePhi.clear();
  v_closeZ.clear();
  v_closeVaneID.clear();

  v_clusterIndex.clear();
  v_clusterX.clear();
  v_clusterY.clear();
  v_clusterPhi.clear();
  v_clusterZ.clear();
  v_clusterVaneID.clear();

  v_gT.clear(); // global time
  v_pT.clear(); // proper time

  h_X.clear();
  h_Y.clear();
  h_Phi.clear();
  h_Z.clear();

  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t main( Int_t argc, Char_t** argv ){
  gROOT->SetBatch(fl_batch);
  TStyle* sty = Style(1);
  TApplication app( "app", &argc, argv );
  if( !(app.Argc()==2) )
    std::cerr << "Wrong input" << std::endl
	      << "Usage : " << app.Argv(0)
	      << " (char*)infilename" << std::endl
	      << "[e.g]" << std::endl
	      << app.Argv(0) << " test.root" << std::endl
	      << std::endl, abort();
  
  Char_t* infilename = app.Argv(1);
  
  TChain* tree_body  = new TChain( "ntupleBody"  );
  TChain* tree_decay = new TChain( "ntupleDecay" );
  tree_body ->Add( infilename );
  tree_decay->Add( infilename );
  std::cout << "[infile] " << infilename << " : "
	    << tree_body ->GetEntries() << " entries(body), "
	    << tree_decay->GetEntries() << " entries(decay)" << std::endl;
  
  set_readbranch_body ( tree_body  );
  set_readbranch_decay( tree_decay );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Make Canvas
  TCanvas* can_1evt = new TCanvas( "can_1evt", "can_1evt", 1200, 750 );
  can_1evt->Divide(3,2);
  can_1evt->Draw();
  Int_t cnt_show = 0;

   // Objects
  TH2D*   hist_Epos_Nhit       = new TH2D( "hist_Epos_Nhit", "E_{e^{+}}v.s.N_{hit};E_{e^{+}} [MeV];N_{hit}", 70, 0, 350, 100, 0, 200 );
  TH1D*   hist_Epos            = new TH1D( "hist_Epos",      "E_{e^{+}};E_{e^{+}} [MeV]",                    50, 0, 350 );
  TH1D*   hist_Nhit            = new TH1D( "hist_Nhit",      "N_{hit};N_{hit}",                              20, 0, 200 );
  TH1D*   hist_Nrec            = new TH1D( "hist_Nrec",      "Rec. Evt;E_{e^{+}} [MeV]; Rec. Events",        50, 0, 350 );
  TH1D*   hist_eff             = new TH1D( "hist_eff",       "Rec. Eff.;E_{e^{+}} [MeV]; Rec. Eff.",         50, 0, 350 );
  TGraph* g_hitpoint_xy_int    = new TGraph();
  TGraph* g_hitpoint_phiz_int  = new TGraph();
  TGraph* g_hitpoint_vanez_int = new TGraph();

  hist_Epos->SetLineColor  (3);
  hist_Nhit->SetLineColor  (3);
  hist_Nrec->SetLineColor  (2);
  hist_eff ->SetLineColor  (2);
  hist_eff ->SetMarkerColor(2);

  // Muon Orbit
  TArc* g_orbit = new TArc( 0, 0, 330 );
  g_orbit->SetFillStyle(0);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Int_t nevt       = tree_body->GetEntries();
  Int_t cnt_signal = 0;
  if( fl_batch==2 ) can_1evt->Print("pic/tracking.pdf[");
  for( Int_t ievt=0; ievt<nevt; ievt++ ){ // START EVENT-LOOP
    
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "+++++++++++++++ ievt = " << ievt << " ++++++++++++++++++++" << std::endl;
    // read event
    ResetObject();
    tree_body ->GetEntry(ievt);
    tree_decay->GetEntry(ievt);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Muon Decay Information
    TGraph* g_decpoint_xy   = new TGraph(); g_decpoint_xy  ->SetMarkerColor(2);
    TGraph* g_decpoint_phiz = new TGraph(); g_decpoint_phiz->SetMarkerColor(2);
    g_decpoint_xy  ->SetPoint( g_decpoint_xy  ->GetN(), td_Dpos_x[0],                      td_Dpos_y[0] );
    g_decpoint_phiz->SetPoint( g_decpoint_phiz->GetN(), phi_uk(td_Dpos_y[0],td_Dpos_x[0]), td_Dpos_z[0] );

    TArrow* g_decvec_xy = new TArrow( td_Dpos_x[0],                   td_Dpos_y[0],
				      td_Dpos_x[0]+100*td_Dmomv_x[1], td_Dpos_y[0]+100*td_Dmomv_y[1],
				      0.01,">"
				      );
    TArrow* g_decvec_phiz = new TArrow( phi_uk(td_Dpos_y[0],td_Dpos_x[0]), td_Dpos_z[0],
					phi_uk(td_Dpos_y[0],td_Dpos_x[0]) + 
					( phi_uk(td_Dpos_y[0]+td_Dmomv_y[1],td_Dpos_x[0]+td_Dmomv_x[1]) - phi_uk(td_Dpos_y[0],td_Dpos_x[0]) )*300,
					td_Dpos_z[0] + 300*td_Dmomv_z[1],
					0.01,">"
					);
    g_decvec_xy  ->SetLineColor(2);
    g_decvec_phiz->SetLineColor(2);
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Hits
    TGraph* g_hitpoint_xy          = new TGraph(); // all (e+,e-,gamma)
    TGraph* g_hitpoint_phiz        = new TGraph(); // all (e+,e-,gamma)
    TGraph* g_hitpoint_vanez       = new TGraph(); // all (e+,e-,gamma)
    TGraph* g_hitpoint_xy_other    = new TGraph(); // except e+
    TGraph* g_hitpoint_phiz_other  = new TGraph(); // except e+
    TGraph* g_hitpoint_vanez_other = new TGraph(); // except e+
    g_hitpoint_xy   ->SetMarkerColor(3);
    g_hitpoint_phiz ->SetMarkerColor(3);
    g_hitpoint_vanez->SetMarkerColor(3);
    g_hitpoint_xy   ->SetMarkerStyle(24);
    g_hitpoint_phiz ->SetMarkerStyle(24);
    g_hitpoint_vanez->SetMarkerStyle(24);
    g_hitpoint_xy_other   ->SetMarkerColor(4);
    g_hitpoint_phiz_other ->SetMarkerColor(4);
    g_hitpoint_vanez_other->SetMarkerColor(4);
    g_hitpoint_xy_other   ->SetMarkerStyle(24);
    g_hitpoint_phiz_other ->SetMarkerStyle(24);
    g_hitpoint_vanez_other->SetMarkerStyle(24);
    g_hitpoint_xy_other   ->SetLineWidth(2);
    g_hitpoint_phiz_other ->SetLineWidth(2);
    g_hitpoint_vanez_other->SetLineWidth(2);

    Int_t cnt_hit=0;
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "Nvec = " << tb_pos_x->size() << std::endl;
    for( Int_t ihit=0; ihit<tb_pos_x->size(); ihit++ ){ // START HIT-LOOP
      if( tb_EachDepE  ->at(ihit)<=0    ) continue; // not zero energy-deposit
      if( tb_bodyStatus->at(ihit)!=0    ) continue; // injection hit-point (veto outgoing hit-point)
      if( tb_bodyTyp   ->at(ihit)<= 100 ) continue; // hit on vane
      if( tb_bodyTyp   ->at(ihit)>=1000 ) continue; // hit on vane

      // input to vector-objects
      v_X.push_back     ( tb_pos_x->at(ihit)                                       );
      v_Y.push_back     ( tb_pos_y->at(ihit)                                       );
      v_Z.push_back     ( tb_pos_z->at(ihit)                                       );
      v_Phi.push_back   ( phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit))            );
      v_VaneID.push_back( GetVaneID(phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit))) );
      v_gT.push_back    ( tb_gtime->at(ihit)                                       );
      v_pT.push_back    ( tb_ptime->at(ihit)                                       );

      g_hitpoint_xy       ->SetPoint( g_hitpoint_xy       ->GetN(), v_X.at     (v_X.size     ()-1), v_Y.at(v_Y.size()-1) );
      g_hitpoint_phiz     ->SetPoint( g_hitpoint_phiz     ->GetN(), v_Phi.at   (v_Phi.size   ()-1), v_Z.at(v_Z.size()-1) );
      g_hitpoint_vanez    ->SetPoint( g_hitpoint_vanez    ->GetN(), v_VaneID.at(v_VaneID.size()-1), v_Z.at(v_Z.size()-1) );
      g_hitpoint_xy_int   ->SetPoint( g_hitpoint_xy_int   ->GetN(), v_X.at     (v_X.size     ()-1), v_Y.at(v_Y.size()-1) );
      g_hitpoint_phiz_int ->SetPoint( g_hitpoint_phiz_int ->GetN(), v_Phi.at   (v_Phi.size   ()-1), v_Z.at(v_Z.size()-1) );
      g_hitpoint_vanez_int->SetPoint( g_hitpoint_vanez_int->GetN(), v_VaneID.at(v_VaneID.size()-1), v_Z.at(v_Z.size()-1) );

      if( tb_pID->at(ihit)!=2 ){
	g_hitpoint_xy_other   ->SetPoint( g_hitpoint_xy_other   ->GetN(), v_X.at     (v_X.size     ()-1), v_Y.at(v_Y.size()-1) );
	g_hitpoint_phiz_other ->SetPoint( g_hitpoint_phiz_other ->GetN(), v_Phi.at   (v_Phi.size   ()-1), v_Z.at(v_Z.size()-1) );
	g_hitpoint_vanez_other->SetPoint( g_hitpoint_vanez_other->GetN(), v_VaneID.at(v_VaneID.size()-1), v_Z.at(v_Z.size()-1) );
      }
      
      if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "              "
									     << std::setw(3) << std::right << cnt_hit << " : pID = "
									     << tb_pID->at(ihit)   << ", (x,y,z) = ("
									     << std::setw(9) << std::right << Form("%.3f",v_X.at     (v_X.size     ()-1)) << ", "
									     << std::setw(9) << std::right << Form("%.3f",v_Y.at     (v_Y.size     ()-1)) << ", "
									     << std::setw(9) << std::right << Form("%.3f",v_Z.at     (v_Z.size     ()-1)) << ", phi = "
									     << std::setw(6) << std::right << Form("%.3f",v_Phi.at   (v_Phi.size   ()-1)) << ", vane-ID = "
									     << std::setw(2) << std::right << Form("%d",  v_VaneID.at(v_VaneID.size()-1)) << ", t(proper) = "
									     << std::setw(7) << std::right << Form("%.5f",v_pT.at    (v_pT.size    ()-1)) << ", t(global) = "
									     << std::setw(8) << std::right << Form("%.5f",v_gT.at    (v_gT.size    ()-1)) << ", Edep = "
									     << std::setw(7) << std::right << Form("%.5f",tb_EachDepE->at(ihit))
									     << std::endl;
      cnt_hit++;
    } // END HIT-LOOP
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ){
      std::cout << "Nhit  = " << cnt_hit        << std::endl;
      std::cout << "E(e+) = " << td_DtEnergy[1] << " MeV" << std::endl;
    }
    
    hist_Epos     ->Fill( td_DtEnergy[1] );
    hist_Nhit     ->Fill( cnt_hit        );
    hist_Epos_Nhit->Fill( td_DtEnergy[1], cnt_hit );
    if( td_DtEnergy[1] > 150 && cnt_hit > 0 ) cnt_signal++;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Hough Transformation (phi-Z)
    HoughTransform( v_Phi,v_Z, h_Phi, h_Z );
    TH2D* hist_hough_phiz = new TH2D("hist_hough_phiz", "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 500, -250, 250 );
    for( Int_t ivec=0; ivec<h_Phi.size(); ivec++ ) hist_hough_phiz->Fill( h_Phi.at(ivec), h_Z.at(ivec) );
    
    Double_t par0;
    Double_t par1;
    HoughFit_One( hist_hough_phiz, par0, par1 );
    TF1* func_hough_phiz = new TF1("func_hough_phiz","[0]+[1]*x", 0.0, 2*TMath::Pi() );
    func_hough_phiz->SetParameter( 0, par0 );
    func_hough_phiz->SetParameter( 1, par1*180/TMath::Pi() );
    
    TH1D* hist_resi_phiz  = new TH1D("hist_resi_phiz",  "Residual(#phi-Z);",  100, -20, 20 );
    Double_t fl_angle = TMath::ATan(par1*180/TMath::Pi())*180/TMath::Pi();
    // calculate residual from the result of Hough-Transformation
    if( TMath::Abs(fl_angle)>89.9 ) GetPol1ResidualY( func_hough_phiz, v_Phi, v_Z, v_residual_PhiZ );
    else                            GetPol1Residual ( func_hough_phiz, v_Phi, v_Z, v_residual_PhiZ );
    // make residual-histogram
    for( Int_t ivec=0; ivec<v_residual_PhiZ.size(); ivec++ ) hist_resi_phiz->Fill( v_residual_PhiZ.at(ivec) );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // select the hit-points close to the Hough-Fit-line
    TGraph* g_hitpoint_xy_close    = new TGraph();
    TGraph* g_hitpoint_phiz_close  = new TGraph();
    TGraph* g_hitpoint_vanez_close = new TGraph();
    g_hitpoint_xy_close   ->SetMarkerColor(3);
    g_hitpoint_phiz_close ->SetMarkerColor(3);
    g_hitpoint_vanez_close->SetMarkerColor(3);
    g_hitpoint_xy_close   ->SetMarkerStyle(20);
    g_hitpoint_phiz_close ->SetMarkerStyle(20);
    g_hitpoint_vanez_close->SetMarkerStyle(20);
    

    //if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Hits Close to Hough-Line>***************" << std::endl;
    for( Int_t ivec=0; ivec<v_residual_PhiZ.size(); ivec++ ){
      if( TMath::Abs(v_residual_PhiZ.at(ivec) - hist_resi_phiz->GetMean()) <= 3*(hist_resi_phiz->GetRMS()) ){
	v_closeIndex.push_back ( ivec              );
	v_closeX.push_back     ( v_X.at     (ivec) );
	v_closeY.push_back     ( v_Y.at     (ivec) );
	v_closeZ.push_back     ( v_Z.at     (ivec) );
	v_closePhi.push_back   ( v_Phi.at   (ivec) );
	v_closeVaneID.push_back( v_VaneID.at(ivec) );
	g_hitpoint_xy_close   ->SetPoint( g_hitpoint_xy_close   ->GetN(), v_closeX.at     (v_closeX.size     ()-1), v_closeY.at(v_closeY.size()-1) );
	g_hitpoint_phiz_close ->SetPoint( g_hitpoint_phiz_close ->GetN(), v_closePhi.at   (v_closePhi.size   ()-1), v_closeZ.at(v_closeZ.size()-1) );
	g_hitpoint_vanez_close->SetPoint( g_hitpoint_vanez_close->GetN(), v_closeVaneID.at(v_closeVaneID.size()-1), v_closeZ.at(v_closeZ.size()-1) );
	/*
	if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                  "
									       << std::setw(10) << std::right << ivec << " : (X,Y,Z,phi,Vane-ID) = ("
									       << std::setw(9)  << std::right << Form("%.3f",v_closeX.at     (v_closeX.size     ()-1)) << ","
									       << std::setw(9)  << std::right << Form("%.3f",v_closeY.at     (v_closeY.size     ()-1)) << ","
									       << std::setw(9)  << std::right << Form("%.3f",v_closeZ.at     (v_closeZ.size     ()-1)) << ","
								       	       << std::setw(6)  << std::right << Form("%.3f",v_closePhi.at   (v_closePhi.size   ()-1)) << ","
								       	       << std::setw(2)  << std::right << Form("%d",  v_closeVaneID.at(v_closeVaneID.size()-1)) << ")"
									       << std::endl;
	*/
      }
    }

    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   =====> #Hits close to Hough-Line : " << v_closeZ.size() << std::endl;
    //if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "fl_angle = " << fl_angle << std::endl;

    // Sorting
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Sorted Hits Close to Hough-Line>***************" << std::endl;
    std::multimap<Int_t,Double_t> tMap_Index;
    std::multimap<Int_t,Double_t> tMap_X;
    std::multimap<Int_t,Double_t> tMap_Y;
    std::multimap<Int_t,Double_t> tMap_Z;
    std::multimap<Int_t,Double_t> tMap_Phi;
    for( Int_t ivec=0; ivec<v_closeVaneID.size(); ivec++ ){
      tMap_Index.insert( std::make_pair(v_closeVaneID.at(ivec),v_closeIndex.at(ivec)) );
      tMap_X.insert    ( std::make_pair(v_closeVaneID.at(ivec),v_closeX.at    (ivec)) );
      tMap_Y.insert    ( std::make_pair(v_closeVaneID.at(ivec),v_closeY.at    (ivec)) );
      tMap_Z.insert    ( std::make_pair(v_closeVaneID.at(ivec),v_closeZ.at    (ivec)) );
      tMap_Phi.insert  ( std::make_pair(v_closeVaneID.at(ivec),v_closePhi.at  (ivec)) );
    }
    v_closeIndex.clear ();
    v_closeX.clear     ();
    v_closeY.clear     ();
    v_closeZ.clear     ();
    v_closePhi.clear   ();
    v_closeVaneID.clear();

    std::multimap<Int_t,Double_t>::iterator it_Index = tMap_Index.begin();
    std::multimap<Int_t,Double_t>::iterator it_X     = tMap_X.begin    ();
    std::multimap<Int_t,Double_t>::iterator it_Y     = tMap_Y.begin    ();
    std::multimap<Int_t,Double_t>::iterator it_Z     = tMap_Z.begin    ();
    std::multimap<Int_t,Double_t>::iterator it_Phi   = tMap_Phi.begin  ();
    while( it_Index != tMap_Index.end() ){ v_closeVaneID.push_back( (*it_Index).first  ); v_closeIndex.push_back( (*it_Index).second ); it_Index++; }
    while( it_X     != tMap_X.end    () ){                                                v_closeX.push_back    ( (*it_X    ).second ); it_X++;     }
    while( it_Y     != tMap_Y.end    () ){                                                v_closeY.push_back    ( (*it_Y    ).second ); it_Y++;     }
    while( it_Z     != tMap_Z.end    () ){                                                v_closeZ.push_back    ( (*it_Z    ).second ); it_Z++;     }
    while( it_Phi   != tMap_Phi.end  () ){                                                v_closePhi.push_back  ( (*it_Phi  ).second ); it_Phi++;   }
    
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ){
      for( Int_t ivec=0; ivec<v_closePhi.size(); ivec++ ){
	std::cout << "                  "
		  << std::setw(3)  << std::right << ivec << " : "
		  << std::setw(3)  << std::right << v_closeIndex.at(ivec) << " : (X,Y,Z,phi,Vane-ID) = ("
		  << std::setw(9)  << std::right << Form("%.3f",v_closeX.at     (ivec)) << ","
		  << std::setw(9)  << std::right << Form("%.3f",v_closeY.at     (ivec)) << ","
		  << std::setw(9)  << std::right << Form("%.3f",v_closeZ.at     (ivec)) << ","
		  << std::setw(6)  << std::right << Form("%.3f",v_closePhi.at   (ivec)) << ","
		  << std::setw(3)  << std::right << Form("%d",  v_closeVaneID.at(ivec)) << ")"
		  << std::endl;
      }
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // clustering
    Int_t group [100] = {0};
    Int_t groupS[100] = {0};
    Int_t groupE[100] = {0};
    Int_t groupnum  = 0;
    Int_t groupflag = 0;
    
    for( Int_t ii=v_closeVaneID[0]; ii<=v_closeVaneID[v_closeVaneID.size()-1];ii++ ){
      //std::cout << "        ii = " << ii << ", groupnum = " << groupnum << ", group[] = " << group[groupnum] << ", groupflag = " << groupflag << std::endl;
      if( tMap_Phi.count(ii)!=1 ){ // no hit or multiple hits in one vane
	if( groupflag==0 && ii>v_closeVaneID[0] ) groupE[groupnum]=ii-1;
	groupflag++;
      }else if( tMap_Phi.count(ii)==1 ){ // move to next group
	if( groupflag>0 ){
	  if( group[groupnum]>0 ) groupnum++;
	  //std::cout << "New Group : " << groupnum << std::endl;
	  if( group[groupnum]==0 ){
	    groupS[groupnum]=ii;
	    group[groupnum]++;
	  }
	}else{ // added to current group
	  if(group[groupnum]==0) groupS[groupnum]=ii;
	  group[groupnum]++;
	  if( ii==v_closeVaneID[v_closeVaneID.size()-1] ) groupE[groupnum]=ii;
	}
	groupflag=0;
      }else{
	group[groupnum]++;
      }
    }
    
    //-get group with maximum number
    Int_t tmpmax=0;
    Int_t maxgroup=0;
    if( group[0]>0 ){
      for( Int_t ii=0; ii<groupnum+1;ii++ ){
	if( tmpmax<group[ii] ){
	  tmpmax=group[ii];
	  maxgroup=ii;
	}
      }
      //std::cout << "group-ID : " << maxgroup << ", #vane : " << tmpmax << ", S = " << groupS[maxgroup] << ", E = " << groupE[maxgroup] << std::endl;
    }

    
    //-cut group with small number & clustering
    //std::vector<Double_t> clsVaneNum;
    //std::vector<Double_t> clsZ;
    if( tmpmax>1 ){
      if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Clusterd seed Hits>***************" << std::endl;
      for( Int_t ii=0; ii<v_closeVaneID.size(); ii++ ){
	if( v_closeVaneID[ii]>=groupS[maxgroup] && v_closeVaneID[ii]<=groupE[maxgroup] ){
	  v_clusterIndex.push_back ( v_closeIndex [ii] );
	  v_clusterX.push_back     ( v_closeX     [ii] );
	  v_clusterY.push_back     ( v_closeY     [ii] );
	  v_clusterPhi.push_back   ( v_closePhi   [ii] );
	  v_clusterZ.push_back     ( v_closeZ     [ii] );
	  v_clusterVaneID.push_back( v_closeVaneID[ii] );
	  if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                 "
										 << std::setw(10) << std::right << v_clusterIndex.at(v_clusterIndex.size()-1) << " : (X,Y,Z,phi,Vane-ID) = ("
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterX.at     (v_clusterX.size     ()-1)) << ","
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterY.at     (v_clusterY.size     ()-1)) << ","
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterZ.at     (v_clusterZ.size     ()-1)) << ","
										 << std::setw(6)  << std::right << Form("%.3f",v_clusterPhi.at   (v_clusterPhi.size   ()-1)) << ","
										 << std::setw(3)  << std::right << Form("%d",  v_clusterVaneID.at(v_clusterVaneID.size()-1)) << ")"
										 << std::endl;
	}
      }
      
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-steering
      std::set<Int_t> extrap_index;
      for( Int_t ivec=0; ivec<v_clusterIndex.size(); ivec++ ) extrap_index.insert( v_clusterIndex.at(ivec) );
      while(1){
	Int_t    edge_vane_L     = v_clusterVaneID.at(0);
	Int_t    pre_edge_vane_L = v_clusterVaneID.at(1);
	Int_t    edge_vane_H     = v_clusterVaneID.at(v_clusterVaneID.size()-1);
	Int_t    pre_edge_vane_H = v_clusterVaneID.at(v_clusterVaneID.size()-2);
	Int_t    edge_Z_L        = v_clusterZ.at(0);
	Int_t    pre_edge_Z_L    = v_clusterZ.at(1);
	Int_t    edge_Z_H        = v_clusterZ.at(v_clusterZ.size()-1);
	Int_t    pre_edge_Z_H    = v_clusterZ.at(v_clusterZ.size()-2);
	Double_t extrap_H = 2*edge_Z_H - pre_edge_Z_H;
	Double_t extrap_L = 2*edge_Z_L - pre_edge_Z_L;
	if( edge_vane_H == n_vane-1 ) edge_vane_H = edge_vane_H - n_vane;
	if( edge_vane_L == 0        ) edge_vane_L = edge_vane_L + n_vane;
	Double_t min_dev_H = 1000;
	Double_t min_dev_L = 1000;
	Double_t tmp_X_H;
	Double_t tmp_X_L;
	Double_t tmp_Y_H;
	Double_t tmp_Y_L;
	Double_t tmp_Z_H;
	Double_t tmp_Z_L;
	Double_t tmp_Phi_H;
	Double_t tmp_Phi_L;
	Double_t tmp_VaneID_H;
	Double_t tmp_VaneID_L;
	Double_t tmp_index_H;
	Double_t tmp_index_L;
	/*
	std::cout << "***********************Edges : " << edge_vane_L << " & " << edge_vane_H << std::endl;
	std::cout << " +++size =  " << v_clusterVaneID.size() << std::endl;
	std::cout << " +++size =  " << extrap_index.size()    << std::endl;
	std::cout << " +++extrap_L = " << extrap_L << std::endl;
	std::cout << " +++extrap_H = " << extrap_H << std::endl;
	*/

	for( Int_t ivane=0; ivane<v_VaneID.size(); ivane++ ){ // START HIT-LOOP for steering
	  if( v_VaneID.at(ivane)!=edge_vane_H + 1 && 
	      v_VaneID.at(ivane)!=edge_vane_L - 1 ) continue;
	  if( extrap_index.count(ivane) ) continue;
	  
	  if( v_VaneID.at(ivane)==edge_vane_H+1 && TMath::Abs( extrap_H - v_Z.at(ivane) ) < min_dev_H ){
	    min_dev_H  = TMath::Abs( extrap_H - v_Z.at(ivane) );
	    tmp_index_H  = ivane;
	    tmp_X_H      = v_X.at     (ivane);
	    tmp_Y_H      = v_Y.at     (ivane);
	    tmp_Phi_H    = v_Phi.at   (ivane);
	    tmp_Z_H      = v_Z.at     (ivane);
	    tmp_VaneID_H = v_VaneID.at(ivane);
	    //std::cout << " hhhhhhhhhhhhhhhhhh  " << ivane << " : vane-ID = " << v_VaneID.at(ivane) << ", Z = " << tmp_Z_H << ", min_dev = " << min_dev_H << std::endl;
	  }

	  if( v_VaneID.at(ivane)==edge_vane_L-1 && TMath::Abs( extrap_L - v_Z.at(ivane) ) < min_dev_L ){
	    min_dev_L    = TMath::Abs( extrap_L - v_Z.at(ivane) );
	    tmp_index_L  = ivane;
	    tmp_X_L      = v_X.at     (ivane);
	    tmp_Y_L      = v_Y.at     (ivane);
	    tmp_Phi_L    = v_Phi.at   (ivane);
	    tmp_Z_L      = v_Z.at     (ivane);
	    tmp_VaneID_L = v_VaneID.at(ivane);
	    //std::cout << " lllllllllllllllllll  " << ivane << " : vane-ID = " << v_VaneID.at(ivane) << ", Z = " << tmp_Z_L << ", min_dev = " << min_dev_L << std::endl;
	  }
	} // END HIT-LOOP for steering
	//std::cout << "min_dev(H&L) = " << min_dev_H << " & " << min_dev_L << std::endl;
	if( min_dev_H < 10 ){
	  //std::cout << "HHHHHHHHHHHHHHH : " << tmp_VaneID_H << std::endl;
	  extrap_index.insert      ( tmp_index_H  );
	  v_clusterIndex.insert    ( v_clusterIndex.begin(), tmp_index_H  );
	  v_clusterX.push_back     ( tmp_X_H      );
	  v_clusterY.push_back     ( tmp_Y_H      );
	  v_clusterPhi.push_back   ( tmp_Phi_H    );
	  v_clusterZ.push_back     ( tmp_Z_H      );
	  v_clusterVaneID.push_back( tmp_VaneID_H );
	  if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                 Added Hits to Cluster(High-side) : "
										 << std::setw(3)  << std::right << v_clusterIndex.at(v_clusterIndex.size()-1) << " : (X,Y,Z,phi,Vane-ID) = ("
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterX.at     (v_clusterX.size     ()-1)) << ","
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterY.at     (v_clusterY.size     ()-1)) << ","
										 << std::setw(9)  << std::right << Form("%.3f",v_clusterZ.at     (v_clusterZ.size     ()-1)) << ","
										 << std::setw(6)  << std::right << Form("%.3f",v_clusterPhi.at   (v_clusterPhi.size   ()-1)) << ","
										 << std::setw(3)  << std::right << Form("%d",  v_clusterVaneID.at(v_clusterVaneID.size()-1)) << "), min_dev = "
										 << min_dev_H
										 << std::endl;
	}
	if( min_dev_L < 10 ){
	  //std::cout << "LLLLLLLLLLLLLLLL : " << tmp_VaneID_L << std::endl;
	  extrap_index.insert   ( tmp_index_L  );
	  v_clusterIndex.insert ( v_clusterIndex.begin(),  tmp_index_L  );
	  v_clusterX.insert     ( v_clusterX.begin(),      tmp_X_L      );
	  v_clusterY.insert     ( v_clusterY.begin(),      tmp_Y_L      );
	  v_clusterPhi.insert   ( v_clusterPhi.begin(),    tmp_Phi_L    );
	  v_clusterZ.insert     ( v_clusterZ.begin(),      tmp_Z_L      );
	  v_clusterVaneID.insert( v_clusterVaneID.begin(), tmp_VaneID_L );
	  if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                           Added Hits to Cluster(Low-side) : "
										 << std::setw(3) << std::right << v_clusterIndex.at(0) << " : (X,Y,Z,phi,Vane-ID) = ("
										 << std::setw(9) << std::right << Form("%.3f",v_clusterX.at     (0)) << ","
										 << std::setw(9) << std::right << Form("%.3f",v_clusterY.at     (0)) << ","
										 << std::setw(9) << std::right << Form("%.3f",v_clusterZ.at     (0)) << ","
										 << std::setw(6) << std::right << Form("%.3f",v_clusterPhi.at   (0)) << ","
										 << std::setw(3) << std::right << Form("%d",  v_clusterVaneID.at(0)) << "), min_dev = "
										 << min_dev_L
										 << std::endl;
	}
	if( !(min_dev_H < 10 || min_dev_L < 10) ){
	  if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                            ====> steering is finished : "
										 << "min_dev(L)=" << min_dev_L << " at Vane-ID=" << edge_vane_L -1 << " & "
										 << "min_dev(H)=" << min_dev_H << " at Vane-ID=" << edge_vane_H +1
										 << std::endl;
										 
	  break;
	}
      }
    }
    
    TGraph* g_hitcluster_vanez = new TGraph();
    g_hitcluster_vanez->SetMarkerColor(3);
    g_hitcluster_vanez->SetMarkerStyle(20);

    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Clustered Hits>***************" << std::endl;
    for( Int_t ivec=0; ivec<v_clusterVaneID.size(); ivec++ ){
      g_hitcluster_vanez->SetPoint( g_hitcluster_vanez->GetN(), v_clusterVaneID.at(ivec), v_clusterZ.at(ivec) );
      if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                 "
									     << std::setw(3)  << std::right << ivec << " : "
									     << std::setw(3)  << std::right << v_clusterIndex.at(ivec) << " : (X,Y,Z,phi,Vane-ID) = ("
									     << std::setw(9)  << std::right << Form("%.3f",v_clusterX.at     (ivec)) << ","
									     << std::setw(9)  << std::right << Form("%.3f",v_clusterY.at     (ivec)) << ","
									     << std::setw(9)  << std::right << Form("%.3f",v_clusterZ.at     (ivec)) << ","
									     << std::setw(6)  << std::right << Form("%.3f",v_clusterPhi.at   (ivec)) << ","
									     << std::setw(3)  << std::right << Form("%d",  v_clusterVaneID.at(ivec)) << ")"
									     << std::endl;
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Judgement of success or false
    std::map<Double_t,Int_t> map_time_sort;
    for( Int_t ihit=0; ihit<(int)v_Z.size(); ihit++ ) map_time_sort.insert( std::pair<Double_t,Int_t>(v_gT.at(ihit), ihit) );

    const Int_t n_primary_hit = 10;
    Int_t fl_success[n_primary_hit] = {0};

    Int_t tmp_cnt = 0;
    for( std::map<Double_t,Int_t>::iterator itime = map_time_sort.begin(); itime != map_time_sort.end(); itime++ ){
      for( Int_t icls=0;icls<v_clusterZ.size();icls++ ){
	if( v_VaneID[itime->second]==v_clusterVaneID[icls] && TMath::Abs(v_Z[itime->second]-v_clusterZ[icls])<1e-5 ) fl_success[tmp_cnt]++; // check if how many hits are detected from hits(1st~3rd)
      }
      tmp_cnt++;
      if( tmp_cnt==n_primary_hit ) break;
    }
    
    Int_t fl_success_integral[n_primary_hit] = {0};
    for( Int_t ip=0; ip<n_primary_hit; ip++ ){
      if( ip ) fl_success_integral[ip] = fl_success_integral[ip-1] + fl_success[ip];
      else     fl_success_integral[ip] =                             fl_success[ip];
    }

    Int_t fl_fin_success = 0; // final judgement of success or false
    if( fl_success_integral[range_success-1]>=threshold_success ) fl_fin_success = 1;

    if( fl_fin_success ){
      hist_Nrec->Fill( td_DtEnergy[1] );
      
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TText* tex = new TText();
    tex->SetTextColor(3);
    tex->SetTextSize(0.05);

    // Draw
    if( ((cnt_show < fl_show || ievt==nevt-1) || fl_batch==2) && td_DtEnergy[1] > th_show_energy ){
      can_1evt->cd(1);
      gPad->DrawFrame(-350,-350,350,350, Form("EvtNo:%d, E(e+)=%.1f MeV, P(e+) = (%.1f, %.1f, %.1f);X [mm];Y [mm]",td_eventNum,td_DtEnergy[1],td_Dmom_x[1],td_Dmom_y[1],td_Dmom_z[1]));
      g_orbit      ->Draw("Lsame");
      g_decpoint_xy->Draw("Psame");
      g_decvec_xy  ->Draw();
      if( g_hitpoint_xy      ->GetN() ) g_hitpoint_xy      ->Draw("Psame");
      if( g_hitpoint_xy_close->GetN() ) g_hitpoint_xy_close->Draw("Psame");
      if( g_hitpoint_xy_other->GetN() ) g_hitpoint_xy_other->Draw("Psame");

      
      can_1evt->cd(2);
      gPad           ->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, Form("EvtNo:%d, E(e+)=%.1f MeV;#phi [rad];Z [mm]",td_eventNum,td_DtEnergy[1]));
      g_decpoint_phiz->Draw("Psame");
      g_decvec_phiz  ->Draw();
      if( g_hitpoint_phiz      ->GetN() ) g_hitpoint_phiz      ->Draw("Psame");
      if( g_hitpoint_phiz_close->GetN() ) g_hitpoint_phiz_close->Draw("Psame");
      if( g_hitpoint_phiz_other->GetN() ) g_hitpoint_phiz_other->Draw("Psame");
      func_hough_phiz      ->Draw("same");
      
      can_1evt->cd(3);
      hist_hough_phiz->Draw("COLZ");
      
      can_1evt->cd(4);
      hist_resi_phiz ->Draw();
      
      can_1evt->cd(5);
      gPad->DrawFrame(0.0,-250,n_vane,250, Form("EvtNo:%d, E(e+)=%.1f MeV;Vane-ID;Z [mm]",td_eventNum,td_DtEnergy[1]));
      if( g_hitcluster_vanez->GetN() ) g_hitcluster_vanez->Draw("Psame");

      can_1evt->cd(6);
      can_1evt->cd(6)->Clear();
      tex->DrawTextNDC( 0.2,0.80, Form("EvtNo = %d", td_eventNum   ) );
      tex->DrawTextNDC( 0.2,0.75, Form("E(e+) = %.1f MeV", td_DtEnergy[1]) );
      tex->DrawTextNDC( 0.2,0.70, Form("P(e+) = (%.1f, %.1f, %.1f) [MeV]",td_Dmom_x[1],td_Dmom_y[1],td_Dmom_z[1]) );
      tex->DrawTextNDC( 0.2,0.65, Form("Nhit(  all  ) = %d",  v_X.size       ()) );
      tex->DrawTextNDC( 0.2,0.60, Form("Nhit( close ) = %d",  v_closeX.size  ()) );
      tex->DrawTextNDC( 0.2,0.55, Form("Nhit(cluster) = %d",  v_clusterX.size()) );
      tex->DrawTextNDC( 0.2,0.50, Form("%s : %d%d%d%d%d %d%d%d%d%d",(fl_fin_success ? "Success" : "False"),
				       fl_success[0],fl_success[1],fl_success[2],fl_success[3],fl_success[4],
				       fl_success[5],fl_success[6],fl_success[7],fl_success[8],fl_success[9])
			);
      if( v_X.size()>0 ) tex->DrawTextNDC( 0.2,0.40, Form("        (X,  Y,  Z,  phi,  Vane-ID)") );
      if( v_X.size()>0 ) tex->DrawTextNDC( 0.2,0.35, Form("1st : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",v_X.at(0),v_Y.at(0),v_Z.at(0),v_Phi.at(0),v_VaneID.at(0)) );
      if( v_X.size()>1 ) tex->DrawTextNDC( 0.2,0.30, Form("2nd : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",v_X.at(1),v_Y.at(1),v_Z.at(1),v_Phi.at(1),v_VaneID.at(1)) );
      if( v_X.size()>2 ) tex->DrawTextNDC( 0.2,0.25, Form("3rd : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",v_X.at(2),v_Y.at(2),v_Z.at(2),v_Phi.at(2),v_VaneID.at(2)) );
      if( v_X.size()>3 ) tex->DrawTextNDC( 0.2,0.20, Form("4th : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",v_X.at(3),v_Y.at(3),v_Z.at(3),v_Phi.at(3),v_VaneID.at(3)) );
      

      if( ievt!=nevt-1 && !gROOT->IsBatch() ){
	can_1evt->Update();
	can_1evt->WaitPrimitive();
      }
      if( fl_batch==2 ) can_1evt->Print("pic/tracking.pdf");

      cnt_show++;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Delete
    if( ievt!=nevt-1 ){
      delete g_decpoint_xy;
      delete g_decpoint_phiz;
      delete g_decvec_xy;
      delete g_decvec_phiz;
      delete g_hitpoint_xy;
      delete g_hitpoint_phiz;
      delete g_hitpoint_vanez;
      delete g_hitpoint_xy_other;
      delete g_hitpoint_phiz_other;
      delete g_hitpoint_vanez_other;
      delete g_hitpoint_xy_close;
      delete g_hitpoint_phiz_close;
      delete g_hitpoint_vanez_close;

      delete g_hitcluster_vanez;
      
      delete hist_hough_phiz;
      delete func_hough_phiz;
      delete hist_resi_phiz;

      delete tex;
    }
    
  } // END EVENT-LOOP
  if( fl_batch==2 ) can_1evt->Print("pic/tracking.pdf]");

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Calculation of Rec. Eff.
  for( Int_t ibin=0; ibin<hist_eff->GetNbinsX(); ibin++ ){
    Double_t rec_eff  = ( hist_Epos->GetBinContent(ibin+1) ? hist_Nrec->GetBinContent(ibin+1)/hist_Epos->GetBinContent(ibin+1) : 0.0 );
    Double_t rec_effE = ( hist_Nrec->GetBinContent(ibin+1) ? rec_eff*sqrt((1-rec_eff)/hist_Nrec->GetBinContent(ibin+1))        : 0.0 );
    hist_eff->SetBinContent( ibin+1, rec_eff  );
    hist_eff->SetBinError  ( ibin+1, rec_effE );
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Draw
  TCanvas* can = new TCanvas( "can", "can", 1600, 750 );
  can->Divide(4,2);
  can->Draw();
  can->cd(1);
  hist_Epos->Draw();
  hist_Nrec->Draw("same");

  can->cd(2);
  hist_Nhit->Draw();

  can->cd(3);
  hist_Epos_Nhit->Draw("COLZ");

  can->cd(4);
  gPad->DrawFrame( 0.0, 0.0, 350, 1.0, Form("Rec. Eff.(Hit definition >= %d/%d);E_{e^{+}} [MeV]; Rec. Eff.",threshold_success,range_success) );
  hist_eff->Draw("same");

  can->cd(5);
  gPad->DrawFrame( -350,-350,350,350, "Hits(x-y) [integral];X [mm];Y [mm]" );
  g_hitpoint_xy_int->Draw("Psame");

  can->cd(6);
  gPad->DrawFrame( 0.0,-250,2.0*TMath::Pi(),250, "Hits(#phi-Z) [integral];#phi;Z [mm]" );
  g_hitpoint_phiz_int->Draw("Psame");

  can->cd(7);
  gPad->DrawFrame( 0.0,-250,n_vane,250, "Hits(#phi-VaneID) [integral];Vane-ID;Z [mm]" );
  g_hitpoint_vanez_int->Draw("Psame");

  can->Update();

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();

  delete tree_body;
  delete tree_decay;
  delete g_orbit;
  
  return 0;
  
}

