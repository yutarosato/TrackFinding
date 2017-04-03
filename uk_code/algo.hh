#ifndef ALGO_HH
#define ALGO_HH 1

#include "root.hh"
#include <vector>

using namespace std;

class Algo{

public:
  Algo();
  ~Algo();

  void HoughTransform(int,int,double*,double*,double*,double*);
  void HoughTransform(vector<double>&, vector<double>&, 
		      vector<double>&, vector<double>&);
  void OneHoughFit(TH2F*,double&,double&);
  void GetPol1Residual(TF1*, vector<double>&,vector<double>&,vector<double>&);
  void GetPol1ResidualY(TF1*, vector<double>&,vector<double>&,vector<double>&);

private:

  vector<double> afterX, afterY;

};

#endif
