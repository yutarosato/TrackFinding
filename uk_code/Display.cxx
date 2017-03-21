#include "Display.hh"

using namespace std;

Display::Display(){};
Display::~Display(){};

void Display::ROOT_Init(int argc, char **argv){
  //gROOT->Reset();
  gStyle->SetOptStat(0);

  app = new TRint("app",&argc, argv);
  //cpalette();
  hesscol();
}

void Display::ROOT_End(){
  app->Run();
}

void Display::MakeCanvas(int Ncan,int x,int y){
  char can_name[20];
  sprintf(can_name,"can%d",Ncan);
  int Unit=300;
  if(x>4||y>4) Unit=200;
  can[Ncan] = new TCanvas(can_name,can_name,x*Unit,y*Unit);
  can[Ncan]->SetFillColor(10);
  if(x>1 || y>1){
    can[Ncan]->Divide(x,y);
  }
}

void Display::MakeFrame(int Nfr, const char* name,const char*title){
  frame[Nfr] = new TH2F(name,title,1000,0,6.28,1000,-250,250);
  frame[Nfr]->SetFillColor(10);
}

void Display::MakeFrame(int Nfr, const char* name,const char*title ,
			double xmin, double xmax, double ymin, double ymax){
  frame[Nfr] = new TH2F(name,title,1000,xmin,xmax,1000,ymin,ymax);
  frame[Nfr]->SetFillColor(10);
}

void Display::MakeGraph(int Ngra,vector<double> &XX,vector<double> &YY){
  int index=XX.size();
  double x[index],y[index];
  char gra_name[20];
  sprintf(gra_name,"gra%d",Ngra);
  for(int i=0;i<index;i++){
    //    cout << X[i] << " " << Y[i] << endl;
    x[i]=XX[i]; y[i]=YY[i];
  }
  graph[Ngra] = new TGraph(index,x,y);
  graph[Ngra]->SetName(gra_name);
  graph[Ngra]->SetMarkerStyle(20);
  graph[Ngra]->SetMarkerSize(.5);
}

void Display::MakeHist(int Nhis, vector<double> &XX){
  int index=XX.size();
  char his_name[20];
  sprintf(his_name,"his%d",Nhis);
  his[Nhis] = new TH1F(his_name,his_name,100,-20,20);
  for(int i=0; i<index; i++){
    his[Nhis]->Fill(XX[i]);
  }
}
  
void Display::Make2DHist(int Nhis, vector<double> &XX, vector<double> &YY){
  int index=XX.size();
  double x[index],y[index];
  char his_name[20];
  sprintf(his_name,"2Dhis%d",Nhis);
  his2D[Nhis] = new TH2F(his_name,his_name,180,0,180,500,-250,250);
  for(int i=0; i<index; i++){
    his2D[Nhis]->Fill(XX[i],YY[i]);
  }
}

void Display::MakePol1(int Np, double par0, double par1, 
		       double s, double e){
  char polname[10];
  sprintf(polname,"pol1_%d",Np);
  pol1[Np] = new TF1(polname,"[0]+[1]*x",s,e);
  pol1[Np]->SetParameter(0,par0);
  pol1[Np]->SetParameter(1,par1);
  pol1[Np]->SetLineWidth(.5);
}

void Display::MakeOrbit(){
  ring = new TArc(0,0,333);
  ring->SetLineColor(4);
  center = new TArc(0,0,70);
  center->SetLineColor(3);

}
