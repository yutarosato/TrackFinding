void EneEff(char* name){
  gROOT->Reset();
  //gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->SetFillColor(10);
  c1->Divide(4,2);

  int enum, hitnum, allhitnum,ok,type,ok5t;
  double InEne,Dtime,pos_xD,pos_yD,pos_zD,mom_xD,mom_yD,mom_zD;
  double stime,etime,esx,esy,esz,esvx,esvy,esvz,eex,eey,eez,eevx,eevy,eevz;
  double theta;
  
  int div=10;
  ifstream ff(name);
  //ifstream ff("EvSt_t5nsec.dat");
  //ifstream ff("EvSt.dat");

  TH1F *HN = new TH1F("HN","",500,0,500);
  HN->SetFillColor(10);

  int eff150=0,eff200=0;
  int eff150all=0,eff200all=0;

  int all[40],aok[40];
  for(int ii=0; ii<40;ii++){
    all[ii]=0; aok[ii]=0;
  }
  while(1){
    ff >> enum >> hitnum >>allhitnum >> InEne >> Dtime >> pos_xD
       >> pos_yD >> pos_zD >> mom_xD >> mom_yD >> mom_zD
       >> stime >> etime >> esx >> esy >> esz >> esvx >> esvy
       >> esvz >> eex >> eey >> eez >> eevx >> eevy >> eevz
       >> theta >> ok;
    if(ff.eof()) break;
    if(allhitnum>3){
      if(InEne>150){
	eff150all++; 
	if(ok) eff150++;
      }
      if(InEne>200){
	eff200all++;
	if(ok) eff200++;
      }

    }
    

    if(allhitnum>0)HN->Fill(allhitnum);
    if(allhitnum>3){
      for(int ii=0;ii<(int)(400/div);ii++){
	if(ii*div<=InEne && (ii+1)*div>InEne){
	  all[ii]++;
	  if(ok)aok[ii]++;
	}
      }
    } 
    
  }
  ff.close();

  double er[40],rate[40];
  for(int ii=0;ii<40;ii++){
    er[ii]=-1000; rate[ii]=-1000;
  }
  for(int ii=0; ii<40; ii++){
    //    cout << all[ii] << endl;
    if(all[ii]>0){
      rate[ii]=(double)aok[ii]/all[ii];
      er[ii]=(double)TMath::Sqrt(aok[ii]*(1-rate[ii]))/all[ii];

    }
    else{
      //      cout << ii << endl;
    }
  }

  TGraphErrors *gl = new TGraphErrors();
  for(int ii=0;ii<40;ii++){
    gl->SetPoint(ii,(2*ii+1)*div/2,rate[ii]);
    gl->SetPointError(ii,div/2,er[ii]);
    cout << ii << " " << (2*ii+1)*div/2 << " " << rate[ii]
	 << " " << div/2 << " " << er[ii] << endl;
  }

  TH2F *waku = new TH2F("waku","",1000,0,350,1000,0,1.2);
  waku->SetFillColor(10);
  waku->Draw();
  waku->SetXTitle("Energy[MeV]");
  waku->SetYTitle("Efficiency");
  //gg->Draw("e");
  gl->Draw("e");

  TLine *ll = new TLine(0,0.9,350,0.9);
  ll->SetLineColor(4);
  ll->Draw("same");
  c1->SetGridy();

  TCanvas *c2 = new TCanvas("c2","",500,500);
  c2->SetFillColor(10);
  c2->Draw();
  HN->Draw();


  cout << "eff150 :" << (double)eff150/eff150all << " (" 
       << (double)TMath::Sqrt(eff150*(1-(double)eff150/eff150all))/eff150all 
       << ") " << endl;
  cout << "eff200 :" << (double)eff200/eff200all << " (" 
       << (double)TMath::Sqrt(eff200*(1-(double)eff200/eff200all))/eff200all 
       << ") " << endl;


}
