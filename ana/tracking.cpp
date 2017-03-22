#include "setting.h"
const Int_t fl_message = 2;

// Objects
std::vector<Double_t> v_X;
std::vector<Double_t> v_Y;
std::vector<Double_t> v_Phi;
std::vector<Double_t> v_Z;
std::vector<Double_t> v_residual_PhiZ;
std::vector<Double_t> v_closePhi;
std::vector<Double_t> v_closeZ;
std::vector<Double_t> v_VaneID;
std::vector<Double_t> v_closeVaneID;
std::vector<Double_t> h_X;
std::vector<Double_t> h_Y;
std::vector<Double_t> h_Phi;
std::vector<Double_t> h_Z;

void ResetObject(){
  v_X.clear();
  v_Y.clear();
  v_Phi.clear();
  v_Z.clear();
  v_residual_PhiZ.clear();
  v_closePhi.clear();
  v_closeZ.clear();
  v_VaneID.clear();
  v_closeVaneID.clear();
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
  std::cout << "[Body ] " << tree_body ->GetEntries() << " entries" << std::endl;
  std::cout << "[Decay] " << tree_decay->GetEntries() << " entries" << std::endl;
  
  set_readbranch_body ( tree_body  );
  set_readbranch_decay( tree_decay );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Make Canvas
  TCanvas* can = new TCanvas( "can", "can", 1200, 750 );
  can->Divide(3,2);
  can->Draw();

  // Muon Orbit
  TArc* g_orbit = new TArc( 0, 0, 330 );
  g_orbit->SetFillStyle(0);

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Int_t nevt       = tree_body->GetEntries();
  Int_t cnt_signal = 0;
  for( Int_t ievt=0; ievt<nevt; ievt++ ){ // START EVENT-LOOP
    
    std::cout << "+++++++++++++++ ievt = " << ievt << " ++++++++++++++++++++" << std::endl;
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
    TGraph* g_hitpoint_xy   = new TGraph();
    TGraph* g_hitpoint_phiz = new TGraph();
    g_hitpoint_xy  ->SetMarkerColor(3);
    g_hitpoint_phiz->SetMarkerColor(3);
    g_hitpoint_xy  ->SetMarkerStyle(24);
    g_hitpoint_phiz->SetMarkerStyle(24);

    Int_t cnt_hit=0;
    std::cout << "Nvec = " << tb_pos_x->size() << std::endl;
    for( Int_t ihit=0; ihit<tb_pos_x->size(); ihit++ ){ // START HIT-LOOP
      if( tb_EachDepE  ->at(ihit)<=0    ) continue; // not zero energy-deposit
      if( tb_bodyStatus->at(ihit)!=0    ) continue; // injection hit-point (veto outgoing hit-point)
      if( tb_bodyTyp   ->at(ihit)<= 100 ) continue; // hit on vane
      if( tb_bodyTyp   ->at(ihit)>=1000 ) continue; // hit on vane
      g_hitpoint_xy  ->SetPoint( g_hitpoint_xy  ->GetN(), tb_pos_x->at(ihit),                            tb_pos_y->at(ihit) );
      g_hitpoint_phiz->SetPoint( g_hitpoint_phiz->GetN(), phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)), tb_pos_z->at(ihit) );

      // input to vector-objects
      v_X.push_back  ( tb_pos_x->at(ihit)                            );
      v_Y.push_back  ( tb_pos_y->at(ihit)                            );
      v_Z.push_back  ( tb_pos_z->at(ihit)                            );
      v_Phi.push_back( phi_uk(tb_pos_y->at(ihit),tb_pos_x->at(ihit)) );
      
      cnt_hit++;
    } // END HIT-LOOP
    std::cout << "Nhit  = " << cnt_hit        << std::endl;
    std::cout << "E(e+) = " << td_DtEnergy[1] << " MeV" << std::endl;

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
    for( Int_t ivec=0; ivec<v_residual_PhiZ.size(); ivec++ ){
      if( fabs(v_residual_PhiZ.at(ivec) - hist_resi_phiz->GetMean()) <= 3*(hist_resi_phiz->GetRMS()) ){
	v_closeZ.push_back  ( v_Z.at  (ivec) );
	v_closePhi.push_back( v_Phi.at(ivec) );
	g_hitpoint_phiz_close->SetPoint( g_hitpoint_phiz_close->GetN(), v_Phi.at(ivec),            v_Z.at(ivec) );
      }
    }

    TGraph* g_digihit_phiz       = new TGraph();
    TGraph* g_digihit_phiz_close = new TGraph();
    g_digihit_phiz      ->SetMarkerStyle(24);
    g_digihit_phiz      ->SetMarkerColor(3);
    g_digihit_phiz_close->SetMarkerColor(3);
    for( Int_t ivec=0; ivec<v_Phi.size();      ivec++ ){
      v_VaneID.push_back( GetVaneID(v_Phi.at(ivec)) );
      g_digihit_phiz->SetPoint( g_digihit_phiz->GetN(), GetVaneID(v_Phi.at(ivec)), v_Z.at(ivec) );
    }
    for( Int_t ivec=0; ivec<v_closePhi.size(); ivec++ ){
      v_closeVaneID.push_back( GetVaneID(v_closePhi.at(ivec)) );
      g_digihit_phiz_close->SetPoint( g_digihit_phiz_close->GetN(), GetVaneID(v_closePhi.at(ivec)), v_closeZ.at(ivec) );
    }

    std::cout << "fl_angle = " << fl_angle << std::endl;

    // Sorting
    std::multimap<Double_t,Double_t> tMap;
    for( Int_t ivec=0; ivec<(Int_t)v_closeVaneID.size(); ivec++ ){
      tMap.insert( std::make_pair(v_closeVaneID.at(ivec),v_closeZ.at(ivec)) );
    }
    v_closeVaneID.clear();
    v_closeZ.clear     ();
    std::multimap<Double_t,Double_t>::iterator it = tMap.begin();
    while( it != tMap.end() ){
      v_closeVaneID.push_back( (*it).first  );
      v_closeZ.push_back     ( (*it).second );
      it++;
    }
    

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Draw
    can->cd(1);
    gPad->DrawFrame(-350,-350,350,350, Form("EvtNo:%d, E(e+)=%.1f MeV, P(e+) = (%.1f, %.1f, %.1f);X [mm];Y [mm]",td_eventNum,td_DtEnergy[1],td_Dmom_x[1],td_Dmom_y[1],td_Dmom_z[1]));
    g_orbit      ->Draw("Lsame");
    g_decpoint_xy->Draw("Psame");
    g_decvec_xy  ->Draw();
    g_hitpoint_xy->Draw("Psame");

    can->cd(2);
    gPad                 ->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, Form("EvtNo:%d, E(e+)=%.1f MeV;#phi [rad];Z [mm]",td_eventNum,td_DtEnergy[1]));
    g_decpoint_phiz      ->Draw("Psame");
    g_decvec_phiz        ->Draw();
    g_hitpoint_phiz      ->Draw("Psame");
    func_hough_phiz      ->Draw("same");
    g_hitpoint_phiz_close->Draw("Psame");

    can->cd(3);
    hist_hough_phiz->Draw("COLZ");

    can->cd(4);
    hist_resi_phiz ->Draw();

    can->cd(5);
    gPad                 ->DrawFrame(0.0,-250,n_vane,250, Form("EvtNo:%d, E(e+)=%.1f MeV;Vane-ID;Z [mm]",td_eventNum,td_DtEnergy[1]));
    g_digihit_phiz      ->Draw("Psame");
    g_digihit_phiz_close->Draw("Psame");

    can->Update();
    can->WaitPrimitive();

    // Delete
    delete g_decpoint_xy;
    delete g_decpoint_phiz;
    delete g_decvec_xy;
    delete g_decvec_phiz;
    delete g_hitpoint_xy;
    delete g_hitpoint_phiz;
    delete g_hitpoint_phiz_close;

    delete hist_hough_phiz;
    delete func_hough_phiz;
    delete hist_resi_phiz;
    
  } // END EVENT-LOOP
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();
  
  return 0;
  
}

