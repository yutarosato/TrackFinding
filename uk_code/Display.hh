#ifndef DISPLAY_HH
#define DISPLAY_HH 1

#include "root.hh"
#include <vector>

using namespace std;

class Display{

public:
  Display();
  ~Display();

  void ROOT_Init(int,char**);
  void ROOT_End();
  void MakeCanvas(int,int,int);
  void MakeFrame(int,const char*,const char*);
  void MakeFrame(int,const char*,const char*,double,double,double,double);
  void MakeGraph(int,vector <double>&,vector <double>&);
  void MakeHist(int,vector <double>&);
  void Make2DHist(int,vector <double>&, vector <double>&);
  void MakePol1(int,double,double,double s=0,double e=2*TMath::Pi());
  void MakeOrbit();

  TRint *app;
  TCanvas *can[10];
  TH2F *frame[10];
  TGraph *graph[10];
  TH1F *his[10];
  TH2F *his2D[10];
  TF1 *pol1[10];
  TArc *ring, *center;

private:

};

#endif
