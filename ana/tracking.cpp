#include "setting.h"

Int_t main( Int_t argc, Char_t** argv ){
  //gROOT->SetBatch(true);
  TStyle* sty = Style();
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
  std::cout << tree_body ->GetEntries() << std::endl;
  std::cout << tree_decay->GetEntries() << std::endl;
  
  set_readbranch_body ( tree_body  );
  set_readbranch_decay( tree_decay );


  
  TCanvas* can = new TCanvas( "can", "can", 800, 400 );
  can->Divide(2,1);
  can->Draw();
  TArc* g_orbit = new TArc( 0, 0, 330 );


  //TH2D* hist = new TH2D( "hist", "hist", 100, -400, 400, 100, -400, 400 );
  Int_t nevt = tree_body->GetEntries();
  for( Int_t ievt=0; ievt < nevt; ievt++ ){
    can->cd(1);
    gPad->DrawFrame(-400,-400,400,400, Form("EvtNo:%d, E(e+)=%.1f MeV;X [mm];Y [mm]",td_eventNum,td_DtEnergy[1]));
    g_orbit->Draw("Lsame");
    can->cd(2);
    gPad->DrawFrame(-TMath::Pi(),-800,TMath::Pi(),800, Form("EvtNo:%d, E(e+)=%.1f MeV;#phi;Z [mm]",td_eventNum,td_DtEnergy[1]));
    
    TGraph* g_decaypoint_xy   = new TGraph(); g_decaypoint_xy  ->SetMarkerColor(2);
    TGraph* g_hitpoint_xy     = new TGraph(); g_hitpoint_xy    ->SetMarkerColor(3);
    TGraph* g_decaypoint_phiz = new TGraph(); g_decaypoint_phiz->SetMarkerColor(2);
    TGraph* g_hitpoint_phiz   = new TGraph(); g_hitpoint_phiz  ->SetMarkerColor(3);
    
    tree_body ->GetEntry(ievt);
    tree_decay->GetEntry(ievt);


    g_decaypoint_xy->SetPoint( g_decaypoint_xy->GetN(), td_Dpos_x[0], td_Dpos_y[0] );
    for( Int_t ihit=0; ihit<tb_pos_x->size(); ihit++ ) g_hitpoint_xy->SetPoint( g_hitpoint_xy->GetN(), tb_pos_x->at(ihit), tb_pos_y->at(ihit) );

    g_decaypoint_phiz->SetPoint( g_decaypoint_phiz->GetN(), atan2(td_Dpos_y[0],td_Dpos_x[0]), td_Dpos_z[0] );
    for( Int_t ihit=0; ihit<tb_pos_x->size(); ihit++ ){
      g_hitpoint_phiz->SetPoint( g_hitpoint_phiz->GetN(), atan2(tb_pos_y->at(ihit),tb_pos_x->at(ihit)), tb_pos_z->at(ihit) );
    }
    

    can->cd(1);
    g_decaypoint_xy->Draw("Psame");
    g_hitpoint_xy  ->Draw("Psame");

    can->cd(2);
    g_decaypoint_phiz->Draw("Psame");
    g_hitpoint_phiz  ->Draw("Psame");

    can->Update();
    can->WaitPrimitive();

    delete g_decaypoint_xy;
    delete g_hitpoint_xy;
    delete g_decaypoint_phiz;
    delete g_hitpoint_phiz;
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  std::cout << "finish" << std::endl;
  if( !gROOT->IsBatch() ) app.Run();
  
  return 0;
  
}

