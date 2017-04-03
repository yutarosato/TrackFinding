{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->SetFillColor(10);
  c1->Draw();


  TH1F *h = new TH1F("h","",2000,0,50000);
  h->SetFillColor(10);

  double tnum=16;
  ifstream ff("ocp_v16.dat");
  double evn,time;

  double hit[10000],tt[10000];

  for(int ii=0; ii<10000;ii++){
    hit[ii]=0; tt[ii]=0;
  }

  while(1){
    ff >> evn >> time;
    if(ff.eof()) break;
    h->Fill(time);
    //    cout << time << " " << (int)time/5 << endl;
    hit[(int)time/5]++;
  }
  ff.close();
  //  cout << hit[0]+hit[1]+hit[2] << endl;
  double er[10000];

  int all=0;
  for(int ii=0; ii<10000;ii++){
    all += hit[ii];
    er[ii]=TMath::Sqrt(hit[ii]*4);
  }
  cout << " all " << all << endl;
  for(int ii=0; ii<10000;ii++){
    tt[ii]=ii*5;
    //    hit[ii]=(double)hit[ii]*4/(1.056e8)*100;
    hit[ii]=(double)hit[ii]*4/(tnum*1152*2*2)*100;
    er[ii]=(double)er[ii]/(tnum*1152*2*2)*100;
    if(hit[ii]<1e-3){
      hit[ii]=100;
      er[ii]=0;
    }
    //    cout << hit[ii] << endl;
  }


  TGraph *gg = new TGraph(10000,tt,hit);
  gg->SetMarkerStyle(7);
  //  h->Draw();

  double ex[10000];
  for(int ii=0; ii<10000;ii++) ex[ii]=2.5;

  TGraphErrors *ge = new TGraphErrors(10000,tt,hit,ex,er);
  ge->SetMarkerStyle(7);
  TH2F *waku = new TH2F("waku","",1000,0,50000,1000,1e-3,10);
  waku->SetFillColor(10);
  waku->Draw();
  waku->SetYTitle("Occupancy [%]");
  waku->SetXTitle("Time [nsec]");
  ge->Draw("same,pe");

  TH1F *plot = new TH1F("plot","",10000,0,50000);
  plot->SetFillColor(10);
  for(int ii=0;ii<10000;ii++){
    plot->Fill(ii*5,hit[ii]);
  }
  //plot->Draw();
  double fsum=0;
  for(int ii=0;ii<10;ii++){
    fsum += hit[ii];
  }
  cout << fsum/10 << endl;

  c1->SetLogy();
  c1->SetGridy(); c1->SetGridx();
}
