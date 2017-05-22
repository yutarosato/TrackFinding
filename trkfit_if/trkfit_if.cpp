#include <iostream>
#include <iomanip>

#include <TROOT.h>
#include <TChain.h>

//****************<Variables for TTree>******************
Int_t                  t_event;    // event number
// result of track finding
Int_t                  t_ntrk;     // #candidate track in one event
Int_t                  t_true_trk; // flag for success/failure of track-finding : 1(success), 0(failure)
std::vector<Double_t>* t_X;        // hit position
std::vector<Double_t>* t_Y;        // hit position
std::vector<Double_t>* t_Z;        // hit position
std::vector<Double_t>* t_PX;       // momentum at hit position
std::vector<Double_t>* t_PY;       // momentum at hit position
std::vector<Double_t>* t_PZ;       // momentum at hit position
std::vector<Double_t>* t_gT;       // hit (global) time
std::vector<Double_t>* t_pT;       // hit (proper) time
std::vector<Double_t>* t_pID;      // hit by 1(e-), 2(e+), 3(gamma)
std::vector<Double_t>* t_EachDepE; // Energy deposit
// generator information
Double_t t_gen_X;  // muon decay position
Double_t t_gen_Y;  // muon decay position
Double_t t_gen_Z;  // muon decay position
Double_t t_gen_gT; // muon decay (global) time
Double_t t_gen_pT; // muon decay (proper) time
Double_t t_gen_PX; // positron momentum
Double_t t_gen_PY; // positron momentum
Double_t t_gen_PZ; // positron momentum

void set_readbranch( TChain* tree ){
  tree->SetBranchAddress( "event",      &t_event    );
  tree->SetBranchAddress( "ntrk",       &t_ntrk     );
  tree->SetBranchAddress( "true_trk",   &t_true_trk );
  tree->SetBranchAddress( "t_X",        &t_X        );
  tree->SetBranchAddress( "t_Y",        &t_Y        );
  tree->SetBranchAddress( "t_Z",        &t_Z        );
  tree->SetBranchAddress( "t_PX",       &t_PX       );
  tree->SetBranchAddress( "t_PY",       &t_PY       );
  tree->SetBranchAddress( "t_PZ",       &t_PZ       );
  tree->SetBranchAddress( "t_gT",       &t_gT       );
  tree->SetBranchAddress( "t_pT",       &t_pT       );
  tree->SetBranchAddress( "t_pID",      &t_pID      );
  tree->SetBranchAddress( "t_EachDepE", &t_EachDepE );
  tree->SetBranchAddress( "gen_X",      &t_gen_X    );
  tree->SetBranchAddress( "gen_Y",      &t_gen_Y    );
  tree->SetBranchAddress( "gen_Z",      &t_gen_Z    );
  tree->SetBranchAddress( "gen_gT",     &t_gen_gT   );
  tree->SetBranchAddress( "gen_pT",     &t_gen_pT   );
  tree->SetBranchAddress( "gen_PX",     &t_gen_PX   );
  tree->SetBranchAddress( "gen_PY",     &t_gen_PY   );
  tree->SetBranchAddress( "gen_PZ",     &t_gen_PZ   );

  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t main( Int_t argc, Char_t** argv ){
  const Char_t* infilename = "test.root";
  TChain* tree = new TChain( "trkfinding" );
  tree->Add( Form("%s",infilename) );
  set_readbranch( tree );
  
  std::cout << tree->GetEntries() << " candidate tracks" << std::endl;
  for( Int_t ievt=0; ievt<tree->GetEntries(); ievt++ ){ // START EVENT-LOOP
    tree->GetEntry(ievt);

    std::cout << std::endl
	      << "*********** [EvtNo = " << t_event << "] **************" << std::endl;
    std::cout << "   ntrk = "     << t_ntrk      << ", "
	      << "   true_trk = " << t_true_trk  << ", "
	      << "   Nhit = "     << t_X->size() << std::endl
	      << Form("   Muon decay position(X,Y,Z)=(%4.2f,%4.2f,%4.2f)", t_gen_X,  t_gen_Y,  t_gen_Z  ) << std::endl
	      << Form("   Positoron momentum (X,Y,Z)=(%4.2f,%4.2f,%4.2f)", t_gen_PX, t_gen_PY, t_gen_PZ ) << std::endl
	      << Form("   Muon decay time (global,proper)=(%4.2f,%4.2f)",  t_gen_gT, t_gen_pT           ) << std::endl;

    for( Int_t ihit=0; ihit<t_X->size(); ihit++ ){ // START HIT-LOOP
      std::cout << std::setw(6) << std::right << ihit
		<< Form(" : (X,Y,Z)=(%4.2f,%4.2f,%4.2f)",    t_X->at(ihit),  t_Y->at(ihit),  t_Z->at(ihit))
		<< Form(", (PX,PY,PZ)=(%4.2f,%4.2f,%4.2f)", t_PX->at(ihit), t_PY->at(ihit), t_PZ->at(ihit))
		<< ", (global time, proper time) = " << Form("(%4.2f,%4.2f)", t_gT ->at(ihit), t_pT ->at(ihit))
		<< ", gen_pID = "  << t_pID     ->at(ihit)
		<< ", EachDepE = " << t_EachDepE->at(ihit)
		<< std::endl;
    } // END HIT-LOOP
  } // END EVENT-LOOP
  
  

  return 0;
  
}

