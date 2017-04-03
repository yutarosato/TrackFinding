{
  gROOT->Reset();
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",1000,600);
  c1->SetFillColor(10);
  c1->Draw();

  ifstream f16("plot_v16.dat");
  ifstream f24("plot_v24.dat");
  ifstream f32("plot_v32.dat");
  ifstream f40("plot_v40.dat");
  ifstream f48("plot_v48.dat");

  int num;
  double xx,yy,ex,ey;

  TGraphErrors *gg[5];
  for(int ii=0; ii<5;ii++){
    gg[ii] = new TGraphErrors();
  }
  for(int ii=0;ii<40;ii++){
    f16 >> num >> xx >> yy >> ex >> ey;
    gg[0]->SetPoint(ii,xx,yy);
    gg[0]->SetPointError(ii,ex,ey);
    f24 >> num >> xx >> yy >> ex >> ey;
    gg[1]->SetPoint(ii,xx,yy);
    gg[1]->SetPointError(ii,ex,ey);
    f32 >> num >> xx >> yy >> ex >> ey;
    gg[2]->SetPoint(ii,xx,yy);
    gg[2]->SetPointError(ii,ex,ey);
    f40 >> num >> xx >> yy >> ex >> ey;
    gg[3]->SetPoint(ii,xx,yy);
    gg[3]->SetPointError(ii,ex,ey);
    f48 >> num >> xx >> yy >> ex >> ey;
    gg[4]->SetPoint(ii,xx,yy);
    gg[4]->SetPointError(ii,ex,ey);
  }

  f16.close(); f24.close(); f32.close(); f40.close(); f48.close();

  TH2F *waku = new TH2F("waku","",1000,0,350,1000,0,1.2);
  waku->SetFillColor(10);
  waku->Draw();
  waku->SetXTitle("Energy[MeV]");
  waku->SetYTitle("Efficiency");

  gg[0]->SetLineColor(1);
  gg[1]->SetLineColor(6);
  gg[2]->SetLineColor(4);
  gg[3]->SetLineColor(3);
  gg[4]->SetLineColor(2);

  for(int ii=0; ii<5; ii++){
    gg[ii]->Draw("same,e");
  }

  

  TLine *ll = new TLine(0,0.9,350,0.9);
  ll->SetLineColor(4);
  ll->Draw("same");
  c1->SetGridy();
}
