#include "algo.hh"

using namespace std;

Algo::Algo(){};
Algo::~Algo(){};


void Algo::HoughTransform(int div, int N,double X[],double Y[], 
			  double R[], double Theta[]){
  for(int ii=0;ii<N;ii++){
    for(int jj=0;jj<div;jj++){
	R[ii*div+jj]=X[ii]*TMath::Cos((double)jj/(double)div*TMath::Pi())
	  +Y[ii]*TMath::Sin((double)jj/(double)div*TMath::Pi());
	Theta[ii*div+jj]=jj;
      
    }
  }


}

void Algo::HoughTransform(vector<double> &X, vector<double> &Y,
			  vector<double> &afX, vector<double> &afY){
  int div=180;
  int N=X.size();
  for(int ii=0;ii<N;ii++){
    for(int jj=0;jj<div;jj++){
      if(Y[ii]>-100 && Y[ii]<100){
	afY.push_back(X[ii]*180/TMath::Pi()*TMath::Cos((double)jj*TMath::Pi()/180)
		      +Y[ii]*TMath::Sin((double)jj*TMath::Pi()/180));
	afX.push_back(jj);
      }
    }
  }
  afterX = afX; afterY=afY;
}

void Algo::OneHoughFit(TH2F *h, double &par0, double &par1){
  int tmpX,tmpY,tmpZ;
  double Rtheta,Rr;
  h->GetMaximumBin(tmpX,tmpY,tmpZ);
  Rtheta = (tmpX+tmpX-1)/2.;
  Rr = (tmpY+tmpY-1)/2.-250;
  //cout << Rtheta << " " << Rr << endl;
  
  par0 = Rr/TMath::Sin(Rtheta*TMath::Pi()/180);
  par1 = -1/TMath::Tan(Rtheta*TMath::Pi()/180);

}

void Algo::GetPol1Residual(TF1 *tf, vector<double> &XX,vector<double> &YY,
			   vector<double> &Resi){
  double par0 = tf->GetParameter(0);
  double par1 = tf->GetParameter(1);

  int num=XX.size();
  for(int ii=0; ii<num; ii++){
    Resi.push_back(YY[ii]-(par0+par1*XX[ii]));
  }

}

void Algo::GetPol1ResidualY(TF1 *tf, vector<double> &XX,vector<double> &YY,
			   vector<double> &Resi){
  double par0 = tf->GetParameter(0);
  double par1 = tf->GetParameter(1);

  int num=XX.size();
  for(int ii=0; ii<num; ii++){
    Resi.push_back( (XX[ii]-((YY[ii]-par0)/par1)) *180/TMath::Pi());
  }

}
			   
