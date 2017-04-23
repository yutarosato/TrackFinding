#ifndef SETTING_H
#define SETTING_H

#include "../Util/Style.h"
#include "../Util/HitsArray.h"

#include <iostream>
#include <iomanip>
#include <map>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TChain.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1C.h>
#include <TH2C.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TArc.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>

#include <TView.h>
#include <TRandom.h>
#include <TPolyLine3D.h>
#include <TAxis3D.h>
#include <TTUBE.h>

#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

// variables for tree(body) [18 branches]
Int_t tb_eventNum;
Int_t tb_hitInfo;
std::vector<Int_t>*    tb_bodyTyp;
std::vector<Int_t>*    tb_bodyStatus;
std::vector<Int_t>*    tb_chID;
std::vector<Int_t>*    tb_pID;
std::vector<Double_t>* tb_CurrentDepE;
std::vector<Double_t>* tb_EachDepE;
std::vector<Double_t>* tb_kEnergy;
std::vector<Double_t>* tb_mom_x;
std::vector<Double_t>* tb_mom_y;
std::vector<Double_t>* tb_mom_z;
std::vector<Double_t>* tb_ptime;
std::vector<Double_t>* tb_gtime;
std::vector<Double_t>* tb_pos_x;
std::vector<Double_t>* tb_pos_y;
std::vector<Double_t>* tb_pos_z;
std::vector<Double_t>* tb_tEnergy;

// variables for tree(decay) [17 branches]
Int_t   td_eventNum;
Float_t td_Dptime  [4];
Float_t td_Dgtime  [4];
Float_t td_DkEnergy[4];
Float_t td_Dmom_x  [4];
Float_t td_Dmom_y  [4];
Float_t td_Dmom_z  [4];
Float_t td_Dmomv_x [4];
Float_t td_Dmomv_y [4];
Float_t td_Dmomv_z [4];
Float_t td_Dpol_x  [4];
Float_t td_Dpol_y  [4];
Float_t td_Dpol_z  [4];
Float_t td_Dpos_x  [4];
Float_t td_Dpos_y  [4];
Float_t td_Dpos_z  [4];
Float_t td_DtEnergy[4];

void set_readbranch_body( TChain* tree ){
  tree->SetBranchAddress( "eventNum",    &tb_eventNum    );
  tree->SetBranchAddress( "hitInfo",     &tb_hitInfo     );
  tree->SetBranchAddress( "bodyTyp",     &tb_bodyTyp     );
  tree->SetBranchAddress( "bodyStatus",  &tb_bodyStatus  );
  tree->SetBranchAddress( "chID",        &tb_chID        );
  tree->SetBranchAddress( "pID",         &tb_pID         );
  tree->SetBranchAddress( "CurrentDepE", &tb_CurrentDepE );
  tree->SetBranchAddress( "EachDepE",    &tb_EachDepE    );
  tree->SetBranchAddress( "kEnergy",     &tb_kEnergy     );
  tree->SetBranchAddress( "mom_x",       &tb_mom_x       );
  tree->SetBranchAddress( "mom_y",       &tb_mom_y       );
  tree->SetBranchAddress( "mom_z",       &tb_mom_z       );
  tree->SetBranchAddress( "ptime",       &tb_ptime       );
  tree->SetBranchAddress( "gtime",       &tb_gtime       );
  tree->SetBranchAddress( "pos_x",       &tb_pos_x       );
  tree->SetBranchAddress( "pos_y",       &tb_pos_y       );
  tree->SetBranchAddress( "pos_z",       &tb_pos_z       );
  tree->SetBranchAddress( "tEnergy",     &tb_tEnergy     );

  return;
}

void set_readbranch_decay( TChain* tree ){
  tree->SetBranchAddress( "eventNum", &td_eventNum );
  tree->SetBranchAddress( "Dptime",    td_Dptime   );
  tree->SetBranchAddress( "Dgtime",    td_Dgtime   );
  tree->SetBranchAddress( "DkEnergy",  td_DkEnergy  );
  tree->SetBranchAddress( "Dmom_x",    td_Dmom_x    );
  tree->SetBranchAddress( "Dmom_y",    td_Dmom_y    );
  tree->SetBranchAddress( "Dmom_z",    td_Dmom_z    );
  tree->SetBranchAddress( "Dmomv_x",   td_Dmomv_x   );
  tree->SetBranchAddress( "Dmomv_y",   td_Dmomv_y   );
  tree->SetBranchAddress( "Dmomv_z",   td_Dmomv_z   );
  tree->SetBranchAddress( "Dpol_x",    td_Dpol_x    );
  tree->SetBranchAddress( "Dpol_y",    td_Dpol_y    );
  tree->SetBranchAddress( "Dpol_z",    td_Dpol_z    );
  tree->SetBranchAddress( "Dpos_x",    td_Dpos_x    );
  tree->SetBranchAddress( "Dpos_y",    td_Dpos_y    );
  tree->SetBranchAddress( "Dpos_z",    td_Dpos_z    );
  tree->SetBranchAddress( "DtEnergy",  td_DtEnergy  );

  return;
}


// output tree for track fitting
Int_t t_event;
Int_t t_ntrk;
Int_t t_true_trk;
std::vector<Double_t> t_X;
std::vector<Double_t> t_Y;
std::vector<Double_t> t_Z;
std::vector<Double_t> t_gT;
std::vector<Double_t> t_pT;
std::vector<Double_t> t_pID;
std::vector<Double_t> t_EachDepE;

Double_t t_gen_X;
Double_t t_gen_Y;
Double_t t_gen_Z;
Double_t t_gen_gT;
Double_t t_gen_pT;
Double_t t_gen_PX;
Double_t t_gen_PY;
Double_t t_gen_PZ;

void set_outtree( TTree* tree ){
  tree->Branch( "event",      &t_event,    "event/I"    );
  tree->Branch( "ntrk",       &t_ntrk,     "ntrk/I"     );
  tree->Branch( "true_trk",   &t_true_trk, "true_trk/I" );
  tree->Branch( "t_X",        &t_X                      );
  tree->Branch( "t_Y",        &t_Y                      );
  tree->Branch( "t_Z",        &t_Z                      );
  tree->Branch( "t_gT",       &t_gT                     );
  tree->Branch( "t_pT",       &t_pT                     );
  tree->Branch( "t_pID",      &t_pID                    );
  tree->Branch( "t_EachDepE", &t_EachDepE               );
  tree->Branch( "gen_X",      &t_gen_X,     "gen_X/D"   );
  tree->Branch( "gen_Y",      &t_gen_Y,     "gen_Y/D"   );
  tree->Branch( "gen_Z",      &t_gen_Z,     "gen_Z/D"   );
  tree->Branch( "gen_gT",     &t_gen_gT,    "gen_gT/D"  );
  tree->Branch( "gen_pT",     &t_gen_pT,    "gen_pT/D"  );
  tree->Branch( "genPX",      &t_gen_PX,    "gen_PX/D"  );
  tree->Branch( "genPY",      &t_gen_PY,    "gen_PY/D"  );
  tree->Branch( "genPZ",      &t_gen_PZ,    "gen_PZ/D"  );

  return;
}

void clear_outtree(){
  t_X.clear();
  t_Y.clear();
  t_Z.clear();
  t_gT.clear();
  t_pT.clear();
  t_pID.clear();
  t_EachDepE.clear();
  return;
}

#endif
