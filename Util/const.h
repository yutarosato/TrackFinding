#ifndef CONST_H
#define CONST_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <set>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TCut.h>
#include <TLeaf.h>
#include <THStack.h>
#include <TH2C.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TFile.h>
#include <TArrow.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSpectrum.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>


const Int_t nvane = 8;
const Int_t nvane_shape = 4;
const Double_t vane_width  = 210; // [mm]
const Double_t vane_height = 700; // [mm]
const Double_t vane_inner  =  75; // [mm]

const Int_t nbase = 3;
const double r_base_point = 200; // mm

#endif
