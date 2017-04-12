#include "setting.h"

const Int_t    fl_message        = 2;
const Int_t    fl_show           = 100;
const Double_t th_show_energy    = 200.0;
const Int_t    threshold_success = 3; // Hit definition : >= threshold_success/range_success
const Int_t    range_success     = 3;
const Int_t    fl_batch          = 0; // 0(show), 1(batch), 2(batch&save)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t main( Int_t argc, Char_t** argv ){
  gROOT->SetBatch(fl_batch);
  TStyle* sty = Style(1);
  TApplication app( "app", &argc, argv );
  if( app.Argc()<2 )
    std::cerr << "Wrong input" << std::endl
	      << "Usage : " << app.Argv(0)
	      << " (char*)infilenames " << std::endl
	      << "[e.g]" << std::endl
	      << app.Argv(0) << " test.root" << std::endl
	      << std::endl, abort();
  
  Char_t* infilename = app.Argv(1);
  
  TChain* tree_body  = new TChain( "ntupleBody"  );
  TChain* tree_decay = new TChain( "ntupleDecay" );
  Int_t nfile = app.Argc()-1;
  for( Int_t ifile=0; ifile<nfile; ifile++ ){
    tree_body ->Add( app.Argv(ifile+1) );
    tree_decay->Add( app.Argv(ifile+1) );
  }
  std::cout << nfile                    << " files, "
	    << tree_body ->GetEntries() << " entries(body), "
	    << tree_decay->GetEntries() << " entries(decay)" << std::endl;
  
  set_readbranch_body ( tree_body  );
  set_readbranch_decay( tree_decay );

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Make Canvas
  TCanvas* can_1evt = new TCanvas( "can_1evt", "can_1evt", 1600, 750 );
  can_1evt->Divide(4,2);
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
  HitsArray* hits_info = new HitsArray();
  
  for( Int_t ievt=0; ievt<nevt; ievt++ ){ // START EVENT-LOOP
    //if( ievt!=669 ) continue; // tmppppp
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ) std::cout << "+++++++++++++++ ievt = " << ievt << " ++++++++++++++++++++" << std::endl;
    // read event
    hits_info ->ClearEvent();
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
    g_hitpoint_xy         ->SetMarkerColor(3);
    g_hitpoint_phiz       ->SetMarkerColor(3);
    g_hitpoint_vanez      ->SetMarkerColor(3);
    g_hitpoint_xy         ->SetMarkerStyle(24);
    g_hitpoint_phiz       ->SetMarkerStyle(24);
    g_hitpoint_vanez      ->SetMarkerStyle(24);
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
      hits_info->InputHits( cnt_hit, tb_pos_x->at(ihit), tb_pos_y->at(ihit), tb_pos_z->at(ihit), tb_ptime->at(ihit), tb_gtime->at(ihit), tb_pID->at(ihit), tb_EachDepE->at(ihit) );
      cnt_hit++;
    } // END HIT-LOOP
    hits_info->CalcOrder();
    if( fl_message && (cnt_show < fl_show || ievt==nevt-1) ){
      std::cout << "Nhit  = " << hits_info->GetNhits()    << std::endl;
      std::cout << "E(e+) = " << td_DtEnergy[1] << " MeV" << std::endl;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    for( Int_t ihit=0; ihit<hits_info->GetNhits(); ihit++ ){
      g_hitpoint_xy       ->SetPoint( g_hitpoint_xy       ->GetN(), hits_info->GetX     (ihit), hits_info->GetY(ihit) );
      g_hitpoint_phiz     ->SetPoint( g_hitpoint_phiz     ->GetN(), hits_info->GetPhi   (ihit), hits_info->GetZ(ihit) );
      g_hitpoint_vanez    ->SetPoint( g_hitpoint_vanez    ->GetN(), hits_info->GetVaneID(ihit), hits_info->GetZ(ihit) );
      g_hitpoint_xy_int   ->SetPoint( g_hitpoint_xy_int   ->GetN(), hits_info->GetX     (ihit), hits_info->GetY(ihit) );
      g_hitpoint_phiz_int ->SetPoint( g_hitpoint_phiz_int ->GetN(), hits_info->GetPhi   (ihit), hits_info->GetZ(ihit) );
      g_hitpoint_vanez_int->SetPoint( g_hitpoint_vanez_int->GetN(), hits_info->GetVaneID(ihit), hits_info->GetZ(ihit) );
      if( hits_info->GetpID(ihit)!=2 ){
	g_hitpoint_xy_other   ->SetPoint( g_hitpoint_xy_other   ->GetN(), hits_info->GetX     (ihit), hits_info->GetY(ihit) );
	g_hitpoint_phiz_other ->SetPoint( g_hitpoint_phiz_other ->GetN(), hits_info->GetPhi   (ihit), hits_info->GetZ(ihit) );
	g_hitpoint_vanez_other->SetPoint( g_hitpoint_vanez_other->GetN(), hits_info->GetVaneID(ihit), hits_info->GetZ(ihit) );
      }
      
    }
    hist_Epos     ->Fill( td_DtEnergy[1] );
    hist_Nhit     ->Fill( hits_info->GetNhits()   );
    hist_Epos_Nhit->Fill( td_DtEnergy[1], hits_info->GetNhits() );
    if( td_DtEnergy[1] > 150 && cnt_hit > 0 ) cnt_signal++;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Hough Transformation (phi-Z)
    hits_info->HoughTransform_phiz();
    hits_info->HoughFit_phiz();
    hits_info->CalcHoughResidual_phiz();

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
    
    for( Int_t ivec=0; ivec<hits_info->GetNhits(); ivec++ ){
      if( hits_info->GetClose_Hough_phiz(ivec)>=0 ){
	g_hitpoint_xy_close   ->SetPoint( g_hitpoint_xy_close   ->GetN(), hits_info->GetX     (ivec), hits_info->GetY(ivec) );
	g_hitpoint_phiz_close ->SetPoint( g_hitpoint_phiz_close ->GetN(), hits_info->GetPhi   (ivec), hits_info->GetZ(ivec) );
	g_hitpoint_vanez_close->SetPoint( g_hitpoint_vanez_close->GetN(), hits_info->GetVaneID(ivec), hits_info->GetZ(ivec) );
      }
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //clustering
    //hits_info->Print_VaneID_Order(fl_message); // tmppppp
    hits_info->Clustering(fl_message *((Int_t)(cnt_show < fl_show)));
    
    TGraph* g_hitcluster_phiz = new TGraph();
    g_hitcluster_phiz->SetMarkerColor(3);
    g_hitcluster_phiz->SetMarkerStyle(20);

    for( Int_t ivec=0; ivec<hits_info->GetNhits(); ivec++ ){
      if( hits_info->GetClusterNo( ivec )!=0 ) continue;
      g_hitcluster_phiz->SetPoint( g_hitcluster_phiz->GetN(), hits_info->GetPhi(ivec), hits_info->GetZ(ivec) );
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Judgement of success or false
    const Int_t n_primary_hit = 10;
    Int_t fl_success[n_primary_hit] = {0};

    Int_t tmp_cnt = 0;
    for( Int_t ip=0; ip<n_primary_hit; ip++ ){
      if( hits_info->GetNhits() == ip ) break;
      if( hits_info->GetClusterNo(hits_info->GetOrdergT(ip)) >= 0 ) fl_success[ip] = 1;
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
    //hits_info->Print             (fl_message*((Int_t)(cnt_show<fl_show))); // tmppppp
    //hits_info->Print_VaneID_Order(fl_message); // tmppppp
    //hits_info->Print_gT_Order    (fl_message); // tmppppp
    // Draw
    if( ((cnt_show < fl_show || ievt==nevt-1) || fl_batch==2) && td_DtEnergy[1] > th_show_energy ){
      //hits_info->Print(fl_message); // tmppppp
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
      hits_info->GetFunc_Hough_phiz(0)->Draw("same");

      can_1evt->cd(3);
      hits_info->GetHist_Hough_phiz(0)->Draw("COLZ");

      can_1evt->cd(4);
      for( Int_t iline=0; iline<hits_info->GetNHoughLines(); iline++ ){
	hits_info->GetHist_Hough_phiz_Residual(iline)->Draw( iline ? "same" : "" );
      }


      can_1evt->cd(5);
      gPad           ->DrawFrame(0.0,-250,2.0*TMath::Pi(),250, Form("EvtNo:%d, E(e+)=%.1f MeV;#phi [rad];Z [mm]",td_eventNum,td_DtEnergy[1]));
      g_decpoint_phiz->Draw("Psame");
      g_decvec_phiz  ->Draw();
      if( g_hitcluster_phiz->GetN() ) g_hitcluster_phiz->Draw("Psame");
      
      can_1evt->cd(8);
      can_1evt->cd(8)->Clear();
      tex->DrawTextNDC( 0.2,0.80, Form("EvtNo = %d", td_eventNum   ) );
      tex->DrawTextNDC( 0.2,0.75, Form("E(e+) = %.1f MeV", td_DtEnergy[1]) );
      tex->DrawTextNDC( 0.2,0.70, Form("P(e+) = (%.1f, %.1f, %.1f) [MeV]",td_Dmom_x[1],td_Dmom_y[1],td_Dmom_z[1]) );
      tex->DrawTextNDC( 0.2,0.65, Form("Nhit(  all  ) = %d",  hits_info->GetNhits()) );
      tex->DrawTextNDC( 0.2,0.50, Form("%s : %d%d%d%d%d %d%d%d%d%d",(fl_fin_success ? "Success" : "False"),
				       fl_success[0],fl_success[1],fl_success[2],fl_success[3],fl_success[4],
				       fl_success[5],fl_success[6],fl_success[7],fl_success[8],fl_success[9])
			);
      if( hits_info->GetNhits()>0 ) tex->DrawTextNDC( 0.2,0.40, Form("        (X,  Y,  Z,  phi,  Vane-ID)") );
      if( hits_info->GetNhits()>0 ) tex->DrawTextNDC( 0.2,0.35, Form("1st : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",
								     hits_info->GetX     (hits_info->GetOrdergT(0)),
								     hits_info->GetY     (hits_info->GetOrdergT(0)),
								     hits_info->GetZ     (hits_info->GetOrdergT(0)),
								     hits_info->GetPhi   (hits_info->GetOrdergT(0)),
								     hits_info->GetVaneID(hits_info->GetOrdergT(0))
								     ));
      if( hits_info->GetNhits()>1 ) tex->DrawTextNDC( 0.2,0.30, Form("2nd : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",
								     hits_info->GetX     (hits_info->GetOrdergT(1)),
								     hits_info->GetY     (hits_info->GetOrdergT(1)),
								     hits_info->GetZ     (hits_info->GetOrdergT(1)),
								     hits_info->GetPhi   (hits_info->GetOrdergT(1)),
								     hits_info->GetVaneID(hits_info->GetOrdergT(1))
								     ));
      if( hits_info->GetNhits()>2 ) tex->DrawTextNDC( 0.2,0.25, Form("3rd : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",
								     hits_info->GetX     (hits_info->GetOrdergT(2)),
								     hits_info->GetY     (hits_info->GetOrdergT(2)),
								     hits_info->GetZ     (hits_info->GetOrdergT(2)),
								     hits_info->GetPhi   (hits_info->GetOrdergT(2)),
								     hits_info->GetVaneID(hits_info->GetOrdergT(2))
								     ));
      if( hits_info->GetNhits()>3 ) tex->DrawTextNDC( 0.2,0.20, Form("4th : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",
								     hits_info->GetX     (hits_info->GetOrdergT(3)),
								     hits_info->GetY     (hits_info->GetOrdergT(3)),
								     hits_info->GetZ     (hits_info->GetOrdergT(3)),
								     hits_info->GetPhi   (hits_info->GetOrdergT(3)),
								     hits_info->GetVaneID(hits_info->GetOrdergT(3))
								     ));
      if( hits_info->GetNhits()>4 ) tex->DrawTextNDC( 0.2,0.15, Form("5th : (%4.2f, %4.2f, %4.2f, %2.2f, %3d)",
								     hits_info->GetX     (hits_info->GetOrdergT(4)),
								     hits_info->GetY     (hits_info->GetOrdergT(4)),
								     hits_info->GetZ     (hits_info->GetOrdergT(4)),
								     hits_info->GetPhi   (hits_info->GetOrdergT(4)),
								     hits_info->GetVaneID(hits_info->GetOrdergT(4))
								     ));
      
      if( ievt!=nevt-1 && !gROOT->IsBatch() ){
	can_1evt->Update();
	//hits_info->Test();
	//hits_info->Print_VaneID_Order(fl_message);
	hits_info->Print_gT_Order(fl_message);
	can_1evt->WaitPrimitive();
      }

      if( fl_batch==2 ) can_1evt->Print("pic/tracking.pdf");

      cnt_show++;
      //if( ievt>100 ) break; // tmppppp
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
      delete g_hitcluster_phiz;
      /*      
      delete hist_hough_phiz;
      delete func_hough_phiz;
      delete hist_resi_phiz;
      */
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

