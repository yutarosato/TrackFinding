#include "setting.h"
const Int_t fl_message = 2;
const Int_t fl_show    = 0;

// Objects
std::vector<Double_t> v_X;
std::vector<Double_t> v_Y;
std::vector<Double_t> v_Phi;
std::vector<Double_t> v_Z;
std::vector<Double_t> v_gT;
std::vector<Double_t> v_pT;
std::vector<Double_t> v_residual_PhiZ;
std::vector<Double_t> v_closePhi;
std::vector<Double_t> v_closeZ;
std::vector<Int_t>    v_VaneID;
std::vector<Int_t>    v_closeVaneID;
std::vector<Double_t> v_clusterZ;
std::vector<Int_t>    v_clusterVaneID;
std::vector<Double_t> h_X;
std::vector<Double_t> h_Y;
std::vector<Double_t> h_Phi;
std::vector<Double_t> h_Z;

void ResetObject(){
  v_X.clear();
  v_Y.clear();
  v_Phi.clear();
  v_Z.clear();
  v_gT.clear(); // global time
  v_pT.clear(); // proper time
  v_residual_PhiZ.clear();
  v_closePhi.clear();
  v_closeZ.clear();
  v_VaneID.clear();
  v_closeVaneID.clear();
  v_clusterZ.clear();
  v_clusterVaneID.clear();
  h_X.clear();
  h_Y.clear();
  h_Phi.clear();
  h_Z.clear();
  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t main( Int_t argc, Char_t** argv ){
  //gROOT->SetBatch(true);
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
  TH1D*   hist_Epos           = new TH1D( "hist_Epos",      "E_{e^{+}};E_{e^{+}} [MeV]", 50, 0, 350 );
  TH1D*   hist_Nhit           = new TH1D( "hist_Nhit",      "N_{hit};N_{hit}",           20, 0, 200 );
  TH2D*   hist_Epos_Nhit      = new TH2D( "hist_Epos_Nhit", "E_{e^{+}}v.s.N_{hit};E_{e^{+}} [MeV];N_{hit}", 50, 0, 350, 20, 0, 200 );
  TH1D*   hist_Nrec           = new TH1D( "hist_Nrec",      "Rec. Evt;E_{e^{+}} [MeV]; Rec. Events", 50, 0, 350 );
  TH1D*   hist_eff            = new TH1D( "hist_eff",       "Rec. Eff.;E_{e^{+}} [MeV]; Rec. Eff.",  50, 0, 350 );
  TGraph* g_hitpoint_xy_int   = new TGraph();
  TGraph* g_hitpoint_phiz_int = new TGraph();

  // Muon Orbit
  TArc* g_orbit = new TArc( 0, 0, 330 );
  g_orbit->SetFillStyle(0);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Int_t nevt       = tree_body->GetEntries();
  Int_t cnt_signal = 0;
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
    TGraph* g_hitpoint_xy         = new TGraph();
    TGraph* g_hitpoint_phiz       = new TGraph();
    TGraph* g_hitpoint_xy_other   = new TGraph();
    TGraph* g_hitpoint_phiz_other = new TGraph();
    g_hitpoint_xy  ->SetMarkerColor(3);
    g_hitpoint_phiz->SetMarkerColor(3);
    g_hitpoint_xy  ->SetMarkerStyle(24);
    g_hitpoint_phiz->SetMarkerStyle(24);
    g_hitpoint_xy_other  ->SetMarkerColor(4);
    g_hitpoint_phiz_other->SetMarkerColor(4);
    g_hitpoint_xy_other  ->SetMarkerStyle(24);
    g_hitpoint_phiz_other->SetMarkerStyle(24);
    g_hitpoint_xy_other  ->SetLineWidth(2);
    g_hitpoint_phiz_other->SetLineWidth(2);

    Int_t cnt_hit=0;
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "Nvec = " << tb_pos_x->size() << std::endl;
    for( Int_t ihit=0; ihit<tb_pos_x->size(); ihit++ ){ // START HIT-LOOP
      if( tb_EachDepE  ->at(ihit)<=0    ) continue; // not zero energy-deposit
      if( tb_bodyStatus->at(ihit)!=0    ) continue; // injection hit-point (veto outgoing hit-point)
      if( tb_bodyTyp   ->at(ihit)<= 100 ) continue; // hit on vane
      if( tb_bodyTyp   ->at(ihit)>=1000 ) continue; // hit on vane
      g_hitpoint_xy      ->SetPoint( g_hitpoint_xy      ->GetN(), tb_pos_x->at(ihit),                            tb_pos_y->at(ihit) );
      g_hitpoint_phiz    ->SetPoint( g_hitpoint_phiz    ->GetN(), phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)), tb_pos_z->at(ihit) );
      g_hitpoint_xy_int  ->SetPoint( g_hitpoint_xy_int  ->GetN(), tb_pos_x->at(ihit),                            tb_pos_y->at(ihit) );
      g_hitpoint_phiz_int->SetPoint( g_hitpoint_phiz_int->GetN(), phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)), tb_pos_z->at(ihit) );
      if( tb_pID->at(ihit)!=2 ){
	g_hitpoint_xy_other  ->SetPoint( g_hitpoint_xy  ->GetN(), tb_pos_x->at(ihit),                            tb_pos_y->at(ihit) );
	g_hitpoint_phiz_other->SetPoint( g_hitpoint_phiz->GetN(), phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)), tb_pos_z->at(ihit) );
      }
      
      // input to vector-objects
      v_X.push_back  ( tb_pos_x->at(ihit)                            );
      v_Y.push_back  ( tb_pos_y->at(ihit)                            );
      v_Z.push_back  ( tb_pos_z->at(ihit)                            );
      v_gT.push_back ( tb_gtime->at(ihit)                            );
      v_pT.push_back ( tb_ptime->at(ihit)                            );
      v_Phi.push_back( phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)) );
      
      cnt_hit++;
      if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "              "
									     << std::setw(3) << std::right << cnt_hit << " : pID = "
									     << tb_pID->at(ihit)   << ", (x,y,z) = ("
									     << std::setw(10) << std::right << Form("%.3f",tb_pos_x->at(ihit)) << ", "
									     << std::setw(10) << std::right << Form("%.3f",tb_pos_y->at(ihit)) << ", "
									     << std::setw(10) << std::right << Form("%.3f",tb_pos_z->at(ihit)) << ", phi = "
									     << std::setw(10) << std::right << Form("%.5f",phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit))) << ", Edep = "
									     << std::setw(10) << std::right << Form("%.7f",tb_EachDepE->at(ihit))                         << ", t(proper) = "
									     << std::setw(10) << std::right << Form("%.7f",tb_ptime->at(ihit))    << ", t(global) = "
									     << std::setw(10) << std::right << Form("%.7f",tb_gtime->at(ihit))
									     << std::endl;
    } // END HIT-LOOP
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ){
      std::cout << "Nhit  = " << cnt_hit        << std::endl;
      std::cout << "E(e+) = " << td_DtEnergy[1] << " MeV" << std::endl;
    }
    hist_Epos     ->Fill( td_DtEnergy[1] );
    hist_Nhit     ->Fill( cnt_hit        );
    hist_Epos_Nhit->Fill( td_DtEnergy[1], cnt_hit );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if( td_DtEnergy[1] > 150 && cnt_hit > 0 ) cnt_signal++;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Hough Transformation (phi-Z)
    HoughTransform( v_Phi,v_Z, h_Phi, h_Z );
    TH2D* hist_hough_phiz = new TH2D("hist_hough_phiz", "Hough(#phi-Z);Hough(#phi) [#circ];Hough(Z) [mm]", 180, 0, 180, 500, -250, 250 );
    for( Int_t ivec=0; ivec<h_Phi.size(); ivec++ ) hist_hough_phiz->Fill( h_Phi.at(ivec), h_Z.at(ivec) );
    
    Double_t par0, par1;
    HoughFit_One( hist_hough_phiz, par0, par1 );
    TF1* func_hough_phiz = new TF1("func_hough_phiz","[0]+[1]*x", 0.0, 2*TMath::Pi() );
    func_hough_phiz->SetParameter( 0, par0 );
    func_hough_phiz->SetParameter( 1, par1*180/TMath::Pi() );
    
    TH1D* hist_resi_phiz  = new TH1D("hist_resi_phiz",  "Residual(#phi-Z);",  100, -20, 20 );
    Double_t fl_angle = TMath::ATan(par1*180/TMath::Pi())*180/TMath::Pi();
    // calculate residual from the result of Hough-Transformation
    if( abs(fl_angle)>89.9 ) GetPol1ResidualY( func_hough_phiz, v_Phi, v_Z, v_residual_PhiZ );
    else                     GetPol1Residual ( func_hough_phiz, v_Phi, v_Z, v_residual_PhiZ );
    // make residual-histogram
    for( Int_t ivec=0; ivec<v_residual_PhiZ.size(); ivec++ ) hist_resi_phiz->Fill( v_residual_PhiZ.at(ivec) );

    // select the hit-points close to the Hough-Fit-line
    TGraph* g_hitpoint_phiz_close  = new TGraph();
    g_hitpoint_phiz_close->SetMarkerColor(3);
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Hits Close to Hough-Line>***************" << std::endl;
    for( Int_t ivec=0; ivec<v_residual_PhiZ.size(); ivec++ ){
      if( fabs(v_residual_PhiZ.at(ivec) - hist_resi_phiz->GetMean()) <= 3*(hist_resi_phiz->GetRMS()) ){
	v_closeZ.push_back  ( v_Z.at  (ivec) );
	v_closePhi.push_back( v_Phi.at(ivec) );
	g_hitpoint_phiz_close->SetPoint( g_hitpoint_phiz_close->GetN(), v_Phi.at(ivec), v_Z.at(ivec) );
	/*
	if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                  "
									       << std::setw(10) << std::right << ivec << " : (Z,phi) = ("
									       << std::setw(10) << std::right << Form("%.5f",v_closeZ.at  (v_closeZ.size  ()-1)) << ","
								       	       << std::setw(10) << std::right << Form("%.5f",v_closePhi.at(v_closePhi.size()-1)) << ")" << std::endl;
	*/
      }
    }

    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   =====> #Hits close to Hough-Line : " << v_closeZ.size() << std::endl;
    TGraph* g_digihit_phiz       = new TGraph();
    TGraph* g_digihit_phiz_close = new TGraph();
    g_digihit_phiz      ->SetMarkerStyle(24);
    g_digihit_phiz      ->SetMarkerColor(3);
    g_digihit_phiz_close->SetMarkerColor(3);

    // convert Z to VaneID
    for( Int_t ivec=0; ivec<v_Phi.size(); ivec++ ){
      v_VaneID.push_back( GetVaneID(v_Phi.at(ivec)) );
      g_digihit_phiz->SetPoint( g_digihit_phiz->GetN(), GetVaneID(v_Phi.at(ivec)), v_Z.at(ivec) );
    }
    // convert Z(close) to VaneID(close)
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Hits Close to Hough-Line with VaneID>***************" << std::endl;
    for( Int_t ivec=0; ivec<v_closePhi.size(); ivec++ ){
      v_closeVaneID.push_back( GetVaneID(v_closePhi.at(ivec)) );
      g_digihit_phiz_close->SetPoint( g_digihit_phiz_close->GetN(), GetVaneID(v_closePhi.at(ivec)), v_closeZ.at(ivec) );
      if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "                  "
									     << std::setw(3)  << std::right << ivec << " : (Z,phi,ID) = ("
									     << std::setw(10) << std::right << Form("%.5f",v_closeZ.at     (ivec)) << ","
									     << std::setw(10) << std::right << Form("%.5f",v_closePhi.at   (ivec)) << ","
									     << std::setw(4)  << std::right << Form("%d",  v_closeVaneID.at(ivec)) << ")"
									     << std::endl;
    }

    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "fl_angle = " << fl_angle << std::endl;

    // Sorting
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "   ************<Sorted Hits Close to Hough-Line with VaneID>***************" << std::endl;
    std::multimap<Double_t,Double_t> tMap;
    for( Int_t ivec=0; ivec<(Int_t)v_closeVaneID.size(); ivec++ ) tMap.insert( std::make_pair(v_closeVaneID.at(ivec),v_closeZ.at(ivec)) );
    v_closeVaneID.clear();
    v_closeZ.clear     ();
    std::multimap<Double_t,Double_t>::iterator it = tMap.begin();
    while( it != tMap.end() ){
      v_closeVaneID.push_back( (*it).first  );
      v_closeZ.push_back     ( (*it).second );
      it++;
    }
    if( fl_message > 1 && (cnt_show < fl_show || ievt==nevt-1) ){
      for( Int_t ivec=0; ivec<v_closePhi.size(); ivec++ ){
	std::cout << "                  "
		  << std::setw(3) << std::right << ivec << " : (Z,phi,ID) = ("
		  << std::setw(10) << std::right << Form("%.7f",v_closeZ.at     (ivec)) << ","
		  << std::setw(10) << std::right << Form("%.7f",v_closePhi.at   (ivec)) << ","
		  << std::setw(4)  << std::right << Form("%d",  v_closeVaneID.at(ivec)) << ")"
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
    
    for( Int_t ii=(Int_t)v_closeVaneID[0]; ii<=(Int_t)v_closeVaneID[v_closeVaneID.size()-1];ii++ ){
      //std::cout << "        ii = " << ii << ", groupnum = " << groupnum << ", group[] = " << group[groupnum] << ", groupflag = " << groupflag << std::endl;
      if( tMap.count(ii)!=1 ){ // no hit or multiple hits in one vane
	if( groupflag==0 && ii>(Int_t)v_closeVaneID[0] ) groupE[groupnum]=ii-1;
	groupflag++;
      }else if( tMap.count(ii)==1 ){ // move to next group
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
	  if( ii==(Int_t)v_closeVaneID[v_closeVaneID.size()-1]) groupE[groupnum]=ii;
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
      //std::cout << "group-ID : " << maxgroup << ", #vane : " << tmpmax << std::endl;
    }
    
    //-cut group with small number & clustering
    //std::vector<Double_t> clsVaneNum;
    //std::vector<Double_t> clsZ;
    if( tmpmax>1 ){
      for( Int_t ii=0; ii<(Int_t)v_closeVaneID.size(); ii++ ){
	if( v_closeVaneID[ii]>=groupS[maxgroup]&&
	    v_closeVaneID[ii]<=groupE[maxgroup]){
	  v_clusterVaneID.push_back(v_closeVaneID[ii]);
	  v_clusterZ.push_back(v_closeZ[ii]);
	}
      }

      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-steering
      //-back
      Double_t hashi,mae,hashiY,maeY,exval;
      Double_t subt,minsub,minsubY;
      Int_t overii[1000*v_VaneID.size()];
      Int_t overiicnt=0;

      for(Int_t ii=0;ii<(Int_t)v_VaneID.size()*1000;ii++) overii[ii]=v_VaneID.size();
      Int_t mplflag=0;
      while(1){
	minsub = 1000;
	hashi  = v_clusterVaneID[v_clusterVaneID.size()-1];
	hashiY = v_clusterZ[v_clusterVaneID.size()-1];
	maeY   = v_clusterZ[v_clusterVaneID.size()-2];
	exval  = hashiY + hashiY - maeY;
	//std::cout << "hashi = " << hashi << std::endl;
	if( hashi==n_vane-1 ) hashi=hashi-n_vane;
	for( Int_t ii=0; ii<(Int_t)v_VaneID.size(); ii++ ){
	  //std::cout << "ii = " << ii << " : " << v_VaneID[ii] << std::endl;
	  if( v_VaneID[ii] == hashi+1 ){
	    //std::cout << "   hashi+1 = " << v_VaneID[ii] << std::endl;
	    subt = TMath::Abs(exval-v_Z[ii]);
	    if(subt<minsub){
	      //std::cout << "         " << subt << " < " << minsub << ", overiicnt = " << overiicnt << std::endl;
	      minsub=subt;
	      minsubY=v_Z[ii];
	      for(Int_t jj=0; jj<overiicnt;jj++){
		if(overii[jj]==ii)minsub=1000;
	      }
	      overii[overiicnt]=ii;
	      overiicnt++;
	    }
	  }
	}
	//std::cout << "minsub = " << minsub << std::endl;
	if( minsub>10 ) break;
	v_clusterVaneID.push_back(hashi+1);
	v_clusterZ.push_back(minsubY);
      }
      for( Int_t ii=0;ii<(Int_t)v_VaneID.size()*1000;ii++ ) overii[ii]=v_VaneID.size();
      overiicnt=0;
      //std::cout << "END BACK" << std::endl;
      //////

      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-front
      while(1){
	minsub = 1000;
	hashi  = v_clusterVaneID[0];
	hashiY = v_clusterZ[0];
	maeY   = v_clusterZ[1];
	exval  = hashiY + hashiY - maeY;
	if( hashi==0 ) hashi=hashi+n_vane;
	for( Int_t ii=0;ii<(Int_t)v_VaneID.size();ii++ ){
	  if( v_VaneID[ii]==hashi-1 ){
	    subt=TMath::Abs(exval-v_Z[ii]);
	    if( subt<minsub ){
	      minsub=subt;
	      minsubY=v_Z[ii];
	      for(Int_t jj=0; jj<overiicnt;jj++){
		if( overii[jj]==ii ) minsub=1000;
	      }
	      overii[overiicnt]=ii;
	      overiicnt++;
	    }
	  }
	}
	if( minsub>10 ) break;
	v_clusterVaneID.insert(v_clusterVaneID.begin(),hashi-1);
	v_clusterZ.insert(v_clusterZ.begin(),minsubY);
      }
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // judgement of success or false
    std::map<Double_t,Int_t> map_time_sort;
    for( Int_t ihit=0; ihit<(int)v_Z.size(); ihit++ ) map_time_sort.insert( std::pair<Double_t,Int_t>(v_gT.at(ihit), ihit) );

    const Int_t n_primary_hit = 10;
    Int_t fl_success[n_primary_hit] = {0};
    

    Int_t tmp_cnt = 0;
    for( std::map<Double_t,Int_t>::iterator itime = map_time_sort.begin(); itime != map_time_sort.end(); itime++ ){
      for( Int_t icls=0;icls<(Int_t)v_clusterZ.size();icls++){ // to be checked . fabs is need ?? tmppppp
	if( v_VaneID[itime->second]==v_clusterVaneID[icls] && v_Z[itime->second]==v_clusterZ[icls] ) fl_success[tmp_cnt]++; // check if how many hits are detected from hits(1st~3rd)
      }
      tmp_cnt++;
      if( tmp_cnt==n_primary_hit ) break;
    }

    Int_t fl_success_integral[n_primary_hit] = {0};
    for( Int_t ip=0; ip<n_primary_hit; ip++ ){
      if( ip ) fl_success_integral[ip] = fl_success_integral[ip-1] + fl_success[ip];
      else     fl_success_integral[ip] =                             fl_success[ip];
    }

    const Int_t index_success = 3;
    if( fl_success_integral[index_success-1]==index_success ){// success
      hist_Nrec->Fill( td_DtEnergy[1] );
    }
    //if( fl_success_integral[0] ){
    //
    //}

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw
    if( cnt_show < fl_show || ievt==nevt-1 ){
      can_1evt->cd(1);
      gPad->DrawFrame(-350,-350,350,350, Form("EvtNo:%d, E(e+)=%.1f MeV, P(e+) = (%.1f, %.1f, %.1f);X [mm];Y [mm]",td_eventNum,td_DtEnergy[1],td_Dmom_x[1],td_Dmom_y[1],td_Dmom_z[1]));
      g_orbit      ->Draw("Lsame");
      g_decpoint_xy->Draw("Psame");
      g_decvec_xy  ->Draw();
      g_hitpoint_xy->Draw("Psame");
      g_hitpoint_xy_other->Draw("Psame");
      
      can_1evt->cd(2);
      gPad                 ->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, Form("EvtNo:%d, E(e+)=%.1f MeV;#phi [rad];Z [mm]",td_eventNum,td_DtEnergy[1]));
      g_decpoint_phiz      ->Draw("Psame");
      g_decvec_phiz        ->Draw();
      g_hitpoint_phiz      ->Draw("Psame");
      func_hough_phiz      ->Draw("same");
      g_hitpoint_phiz_close->Draw("Psame");
      g_hitpoint_phiz_other->Draw("Psame");
      
      can_1evt->cd(3);
      hist_hough_phiz->Draw("COLZ");
      
      can_1evt->cd(4);
      hist_resi_phiz ->Draw();
      
      can_1evt->cd(5);
      gPad                 ->DrawFrame(0.0,-250,n_vane,250, Form("EvtNo:%d, E(e+)=%.1f MeV;Vane-ID;Z [mm]",td_eventNum,td_DtEnergy[1]));
      g_digihit_phiz      ->Draw("Psame");
      g_digihit_phiz_close->Draw("Psame");
      
      can_1evt->Update();
      if( ievt!=nevt-1 ) can_1evt->WaitPrimitive();
      cnt_show++;
    }

    // Delete
    if( ievt!=nevt-1 ){
      delete g_decpoint_xy;
      delete g_decpoint_phiz;
      delete g_decvec_xy;
      delete g_decvec_phiz;
      delete g_hitpoint_xy;
      delete g_hitpoint_phiz;
      delete g_hitpoint_xy_other;
      delete g_hitpoint_phiz_other;
      delete g_hitpoint_phiz_close;
      delete g_digihit_phiz;
      delete g_digihit_phiz_close;
      
      delete hist_hough_phiz;
      delete func_hough_phiz;
      delete hist_resi_phiz;
    }
    
  } // END EVENT-LOOP

  for( Int_t ibin=0; ibin<hist_eff->GetNbinsX(); ibin++ ){
    Double_t rec_eff  = ( hist_Epos->GetBinContent(ibin+1) ? hist_Nrec->GetBinContent(ibin+1)/hist_Epos->GetBinContent(ibin+1) : 0.0 );
    Double_t rec_effE = ( hist_Nrec->GetBinContent(ibin+1) ? rec_eff*sqrt((1-rec_eff)/hist_Nrec->GetBinContent(ibin+1))        : 0.0 );
    hist_eff->SetBinContent( ibin+1, rec_eff  );
    hist_eff->SetBinError  ( ibin+1, rec_effE );
  }

  TCanvas* can = new TCanvas( "can", "can", 1200, 750 );
  can->Divide(3,2);
  can->Draw();
  can->cd(1);
  hist_Epos->Draw();
  hist_Nrec->Draw("same");

  can->cd(2);
  hist_Nhit->Draw();

  can->cd(3);
  hist_Epos_Nhit->Draw("COLZ");

  can->cd(4);
  hist_eff->Draw();

  can->cd(5);
  gPad->DrawFrame( -350,-350,350,350, "Hits(x-y) [integral];X [mm];Y [mm]" );
  g_hitpoint_xy_int->Draw("Psame");

  can->cd(6);
  gPad->DrawFrame( 0.0,-250,2.0*TMath::Pi(),250, "Hits(#phi-Z) [integral];#phi;Z [mm]" );
  g_hitpoint_phiz_int->Draw("Psame");

  can->Update();

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();

  delete tree_body;
  delete tree_decay;
  delete g_orbit;
  
  return 0;
  
}

