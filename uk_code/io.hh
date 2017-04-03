#ifndef IO_HH
#define IO_HH 1

#include "root.hh"
#include <vector>

using namespace std;

class IO{

public:
  IO();
  ~IO();

  void OpenFile();
  void CloseFile();
  void OpenFile(int);
  void CloseFile(int);
  void OpenTree();
  int GetEventNumber();
  void ClearVector();
  void ReadOneEvent(int);
  void GetEvtStrParam(vector<double>&);
  void GetParam(vector<double>&, vector<double>&, vector<double>&,
	       vector<double>&, vector<double>&,
		vector<int>&, vector<int>&, vector<double>&);
  void GetDecayPoint(double&,double&,double&);
  double GetPositronInitialEnergy();
  void ReadData(vector<double>&, vector<double>&, vector<double>&,
	       vector<double>&, vector<double>&,
	       vector<int>&, vector<int>&);
  void DigitizeVane(int, vector<double>&, vector<double>&);
  void AddGhost(int, vector<double>&, vector<double>&, 
		vector<double>&, vector<int>&,
		vector<double>&, vector<double>&);
  int GetUnit(double, double);
  void GetStrips(int, double, double, int&, int&);
  bool CheckCoinci(int,int a[][16]);
  void AddPoints(int, vector<int>&, vector<int>&, 
		 int *a, int b[][100], vector<int>&, vector<double>&);
  void ConvStripToReal(vector<int>&, vector<int>&,
		       vector<int>&, vector<double>&, vector<double>&);


  ifstream data;


private:

  int index;
  TFile *file;
  TTree *TD,*T;

  vector<double> XX, YY, ZZ, RR, PPhi, TTime;
  vector<int> IId, TThr;
  vector<double> EvStr;

  //for TD
  int eventNumD;
  float tEnergyD[4];//0=>mu+, 1=>e+
  float momv_xD[4],momv_yD[4],momv_zD[4];
  float mom_xD[4],mom_yD[4],mom_zD[4];
  float pol_xD[4],pol_yD[4],pol_zD[4];
  float ptimeD[4];//nsec
  float gtimeD[4];//nsec
  float pos_xD[4],pos_yD[4],pos_zD[4];
  int neventD;

  //for T
  int eventNum,hitInfo;
  vector<double> *tEnergy;
  vector<double> *mom_x;
  vector<double> *mom_y;
  vector<double> *mom_z;
  vector<double> *pos_x;
  vector<double> *pos_y;
  vector<double> *pos_z;
  vector<double> *EachDepE;
  vector<double> *kEnergy;
  vector<double> *CurrentDepE;
  vector<int> *bodyTyp;
  vector<int> *bodyStatus;
  vector<int> *pID;
  vector<double> *gtime;
  TBranch *btEnergy;
  TBranch *bmom_x;
  TBranch *bmom_y;
  TBranch *bmom_z;
  TBranch *bpos_x;
  TBranch *bpos_y;
  TBranch *bpos_z;
  TBranch *bbodyTyp;
  TBranch *bbodyStatus;
  TBranch *bpID;
  TBranch *bgtime;
  TBranch *bEachDepE;
  TBranch *bkEnergy;
  TBranch *bCurrentDepE;


};

#endif
