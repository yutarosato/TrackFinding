{
  gROOT->Reset();
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  gStyle->SetStatH(0.3);
  TCanvas *c1 = new TCanvas("c1","",900,600);
  c1->SetFillColor(10);
  c1->Draw();


  TH1F *h = new TH1F("h","",2000,0,50000);
  h->SetFillColor(10);

  ifstream ff("newocp_v48.dat");
  double evn,time,R,P,X,Y,Z;
  const int NN= 10000,nn= 5000;

  double hit[NN],tt[NN],newHit[NN];
  double xx[NN][nn],yy[NN][nn],zz[NN][nn],rr[NN][nn],phi[NN][nn];

  cout << "start!" << endl;
  for(int ii=0; ii<NN;ii++){
    hit[ii]=0; tt[ii]=0; newHit[ii]=0;
    //for(int jj=0; jj<nn;jj++){
    //  xx[ii][jj]=-1000; yy[ii][jj]=-1000; zz[ii][jj]=-1000;
    //  rr[ii][jj]=-1000; phi[ii][jj]=-1000;
    //}
  }
  TH1F *VZ0[48],*VZ1[48],VZ2[48],VZ3[48];
  TH1F *VR0[48],*VR1[48],VR2[48];
  char fname[30];
  for(int ii=0; ii<48;ii++){
    sprintf(fname,"vane_z0_%d",ii);
    VZ0[ii]= new TH1F(fname,"",1100,0,1100);
    sprintf(fname,"vane_z1_%d",ii);
    VZ1[ii]= new TH1F(fname,"",1100,0,1100);
    sprintf(fname,"vane_z2_%d",ii);
    VZ2[ii]= new TH1F(fname,"",1100,0,1100);
    sprintf(fname,"vane_z3_%d",ii);
    VZ3[ii]= new TH1F(fname,"",1100,0,1100);
    sprintf(fname,"vane_r0_%d",ii);    
    VR0[ii]= new TH1F(fname,"",2000,0,2000);
    sprintf(fname,"vane_r1_%d",ii);    
    VR1[ii]= new TH1F(fname,"",2000,0,2000);
    sprintf(fname,"vane_r2_%d",ii);    
    VR2[ii]= new TH1F(fname,"",2000,0,2000);
  }


  cout << "1" << endl;
  double sum=0,sum2=0;
  TGraph *tg = new TGraph();
  int tind=0;
  double vane,rdiv,zdiv;
  while(1){
    ff >> evn >> time >> R >> P >> X >> Y >> Z;
    if(ff.eof()) break;
    h->Fill(time);

    //R digitize
    if(R==290)rdiv=1099;
    else rdiv=(int)((R-70)/0.2);
    //phi digitize (vane)
    if((int)((P*180/TMath::Pi()+3.75)/7.5)==48) vane=0;
    else vane = (int)((P*180/TMath::Pi()+3.75)/7.5);
    //Z digitize
    zdiv = (int)((Z+200)/0.2);


    if((int)time/5==0 ){
      //cout << evn << " " << (int)time/5 << " " << zdiv <<
      //" " << rdiv << " " << X << " " << Y << endl;
      tg->SetPoint(tind,X,Y);
      tind++;
    }
    //Def. of Region
    // Z1->zdiv[0,499],  Z2->[500,999], Z3->[1000,1499], Z4->[1500,1999]
    // r1->rdiv[0,366], r2->[367,733], r3->[734,1099]

    xx[(int)time/5][hit[(int)time/5]]=X;
    yy[(int)time/5][hit[(int)time/5]]=Y;
    zz[(int)time/5][hit[(int)time/5]]=zdiv;
    rr[(int)time/5][hit[(int)time/5]]=rdiv;
    phi[(int)time/5][hit[(int)time/5]]=vane;
    hit[(int)time/5]++;
  }
  ff.close();
  double prer=0,prez=0,prer2=0,prez2=0;
  double preZ0[48],preZ1[48],preZ2[48],preZ3[48];
  double preR0[48],preR1[48],preR2[48];

  for(int ii=0;ii<NN;ii++){
    if(hit[ii]>0){
      for(int jj=0;jj<hit[ii];jj++){
	if(zz[ii][jj]<500){
	  if(rr[ii][jj]!=prer && zz[ii][jj]!=prez){
	    VZ0[(int)phi[ii][jj]]->Fill(rr[ii][jj]);
	    newHit[ii]++;
	  }
	  prer=rr[ii][jj]; prez=zz[ii][jj];
	}
	else if(zz[ii][jj]>=500 && zz[ii][jj]<1000){
	  if(rr[ii][jj]!=prer && zz[ii][jj]!=prez){
	    VZ1[(int)phi[ii][jj]]->Fill(rr[ii][jj]);
	    newHit[ii]++;
	  }
	  prer=rr[ii][jj]; prez=zz[ii][jj];	  
	}
	else if(zz[ii][jj]>=1000 && zz[ii][jj]<1500){
	  if(rr[ii][jj]!=prer && zz[ii][jj]!=prez){
	    VZ2[(int)phi[ii][jj]]->Fill(rr[ii][jj]);
	    newHit[ii]++;
	  }
	  prer=rr[ii][jj]; prez=zz[ii][jj];	  
	}
	else if(zz[ii][jj]>=1500){
	  if(rr[ii][jj]!=prer && zz[ii][jj]!=prez){
	    VZ3[(int)phi[ii][jj]]->Fill(rr[ii][jj]);
	    newHit[ii]++;
	  }
	  prer=rr[ii][jj]; prez=zz[ii][jj];	  
	}

	if(rr[ii][jj]<367){
	  if(rr[ii][jj]!=prer2 && zz[ii][jj]!=prez2){
	    VR0[(int)phi[ii][jj]]->Fill(zz[ii][jj]);
	    //newHit[ii]++;
	  }
	  prer2=rr[ii][jj]; prez2=zz[ii][jj];	  
	}
	else if(rr[ii][jj]>=367&&rr[ii][jj]<734){
	  if(rr[ii][jj]!=prer2 && zz[ii][jj]!=prez2){
	    VR1[(int)phi[ii][jj]]->Fill(zz[ii][jj]);
	    // newHit[ii]++;
	  }
	  prer2=rr[ii][jj]; prez2=zz[ii][jj];	  
	} 
	else if(rr[ii][jj]>=734){
	  if(rr[ii][jj]!=prer2 && zz[ii][jj]!=prez2){
	    VR2[(int)phi[ii][jj]]->Fill(zz[ii][jj]);
	    //	    newHit[ii]++;
	  }
	  prer2=rr[ii][jj]; prez2=zz[ii][jj];	  
	}
      }
    }
  }
    //Def. of Region
    // Z1->zdiv[0,499],  Z2->[500,999], Z3->[1000,1499], Z4->[1500,1999]
    // r1->rdiv[0,366], r2->[367,733], r3->[734,1099]

  TH1F *HNum[12];
  char hnumname[20];
  for(int ii=0; ii<12; ii++){
    sprintf(hnumname,"hithum%d",ii);
    HNum[ii] = new TH1F(hnumname,"",20,0,20);
    HNum[ii]->SetFillColor(10);
    HNum[ii]->SetLineWidth(2);
  }

  for(int ii=0;ii<48;ii++){
    for(int jj=1;jj<1101;jj++){
      if(jj<368){
	HNum[0]->Fill(VZ0[ii]->GetBinContent(jj));
	HNum[3]->Fill(VZ1[ii]->GetBinContent(jj));
	HNum[6]->Fill(VZ2[ii]->GetBinContent(jj));
	HNum[9]->Fill(VZ3[ii]->GetBinContent(jj));
      }
      else if(jj>=368&&jj<735){
	HNum[1]->Fill(VZ0[ii]->GetBinContent(jj));
	HNum[4]->Fill(VZ1[ii]->GetBinContent(jj));
	HNum[7]->Fill(VZ2[ii]->GetBinContent(jj));
	HNum[10]->Fill(VZ3[ii]->GetBinContent(jj));
      }
      else if(jj>=735){
	HNum[2]->Fill(VZ0[ii]->GetBinContent(jj));
	HNum[5]->Fill(VZ1[ii]->GetBinContent(jj));
	HNum[8]->Fill(VZ2[ii]->GetBinContent(jj));
	HNum[11]->Fill(VZ3[ii]->GetBinContent(jj));
      }

    }
  }



  cout << "2" << endl;
  int allhitnum=0,allhitnumred=0;

  double er[NN],ex[NN];
  for(int ii=0; ii<NN;ii++){
    allhitnum += hit[ii];
    allhitnumred += newHit[ii];
    tt[ii]=ii*5;
    er[ii]=(double)TMath::Sqrt(newHit[ii])/(48*2304*2)*100;
    newHit[ii]=(double)newHit[ii]/(48*2304*2)*100;
    ex[ii]=2.5;
    if(newHit[ii]<1e-5)newHit[ii]=1000;
  }
  cout << "All : " << allhitnum << " " << allhitnumred << endl;
  
  TGraphErrors *gg = new TGraphErrors(NN,tt,newHit,ex,er);
  gg->SetMarkerStyle(7);
  //  h->Draw();

  TH2F *waku = new TH2F("waku","",1000,0,50000,1000,1e-4,10);
  waku->SetFillColor(10);
  waku->Draw();
  waku->SetYTitle("Occupancy [%]");
  waku->SetXTitle("Time [nsec]");
  gg->Draw("same,pe");
  c1->SetLogy();
  c1->SetGridy(); c1->SetGridx();

  TCanvas *c2 = new TCanvas("c2","Zdiv",1000,800);
  c2->SetFillColor(10);
  c2->Divide(1,4);
  c2->Draw();
  c2->cd(4);  VZ0[0]->Draw();
  c2->cd(3);  VZ1[0]->Draw();
  c2->cd(2);  VZ2[0]->Draw();
  c2->cd(1);  VZ3[0]->Draw();

  TCanvas *c3 = new TCanvas("c3","Rdiv",1000,800);
  c3->SetFillColor(10);
  c3->Divide(1,3);
  c3->Draw();
  c3->cd(1);  VR0[0]->Draw();
  c3->cd(2);  VR1[0]->Draw();
  c3->cd(3);  VR2[0]->Draw();

  TCanvas *c4 = new TCanvas("c4","Hitnum",1000,900);
  c4->SetFillColor(10);
  c4->Divide(3,4);
  c4->Draw();

  char term[50];

  c4->cd(10);  HNum[0]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[0]->GetEntries());
  TF1 *pois0 = new TF1("pois0",term,0,20);
  pois0->SetParameter(0,3);
  HNum[0]->Fit("pois0","","E",0,20);

  c4->cd(11);  HNum[1]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[1]->GetEntries());
  TF1 *pois1 = new TF1("pois1",term,0,20);
  pois1->SetParameter(0,3);
  HNum[1]->Fit("pois1","","E",0,20);

  c4->cd(12);  HNum[2]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[2]->GetEntries());
  TF1 *pois2 = new TF1("pois2",term,0,20);
  pois2->SetParameter(0,3);
  HNum[2]->Fit("pois2","","E",0,20);


  c4->cd(7);  HNum[3]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[3]->GetEntries());
  TF1 *pois3 = new TF1("pois3",term,0,20);
  pois3->SetParameter(0,3);
  HNum[3]->Fit("pois3","","E",0,20);

  c4->cd(8);  HNum[4]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[4]->GetEntries());
  TF1 *pois4 = new TF1("pois4",term,0,20);
  pois4->SetParameter(0,3);
  HNum[4]->Fit("pois4","","E",0,20);

  c4->cd(9);  HNum[5]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[5]->GetEntries());
  TF1 *pois5 = new TF1("pois5",term,0,20);
  pois5->SetParameter(0,3);
  HNum[5]->Fit("pois5","","E",0,20);
  double poissum=0;
  for(int ii=0; ii<15;ii++){
    poissum += TMath::Poisson(ii,pois5->GetParameter(0));
  }
  cout << poissum << endl;



  c4->cd(4);  HNum[6]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[6]->GetEntries());
  TF1 *pois6 = new TF1("pois6",term,0,20);
  pois6->SetParameter(0,3);
  HNum[6]->Fit("pois6","","E",0,20);

  c4->cd(5);  HNum[7]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[7]->GetEntries());
  TF1 *pois7 = new TF1("pois7",term,0,20);
  pois7->SetParameter(0,3);
  HNum[7]->Fit("pois7","","E",0,20);

  c4->cd(6);  HNum[8]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[8]->GetEntries());
  TF1 *pois8 = new TF1("pois8",term,0,20);
  pois8->SetParameter(0,3);
  HNum[8]->Fit("pois8","","E",0,20);

  c4->cd(1);  HNum[9]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[9]->GetEntries());
  TF1 *pois9 = new TF1("pois9",term,0,20);
  pois9->SetParameter(0,3);
  HNum[9]->Fit("pois9","","E",0,20);

  c4->cd(2);  HNum[10]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[10]->GetEntries());
  TF1 *pois10 = new TF1("pois10",term,0,20);
  pois10->SetParameter(0,3);
  HNum[10]->Fit("pois10","","E",0,20);

  c4->cd(3);  HNum[11]->Draw("e");
  sprintf(term,"%f*TMath::Poisson(x,[0])",HNum[11]->GetEntries());
  TF1 *pois11 = new TF1("pois11",term,0,20);
  pois11->SetParameter(0,3);
  HNum[11]->Fit("pois11","","E",0,20);


  //TF1 *pois = new TF1("pois","17616*TMath::Poisson(x,[0])",0,20);
  //c4->cd(2);
  //pois->SetParameter(0,1000);
  //pois->SetParameter(0,3);
  //HNum[10]->Fit("pois","","E",0,20);

  /*
  TH1F *aa = new TH1F("aa","",30,0,30);
  aa->SetFillColor(10);
  TRandom3 rnd;
  for(int ii=0; ii<17616; ii++){
    aa->Fill(rnd.Poisson(2.8));
  }
  c4->cd(1); aa->Draw("same");
  */

  /*
  TCanvas *can[48*2];
  char canname[30];
  for(int ii=0; ii<48;ii++){
    sprintf(canname,"Z_Vane%d",ii);
    can[ii] = new TCanvas(canname,canname,1000,800);
    can[ii]->SetFillColor(10);
    can[ii]->Divide(1,4);
    //can[ii]->Draw();
    cout << ii << " Z ";
    can[ii]->cd(4); VZ0[ii]->Draw(); 
    cout << VZ0[ii]->GetSum()/1100 << " ";
    can[ii]->cd(3); VZ1[ii]->Draw();
    cout << VZ1[ii]->GetSum()/1100 << " ";
    can[ii]->cd(2); VZ2[ii]->Draw();
    cout << VZ2[ii]->GetSum()/1100 << " ";
    can[ii]->cd(1); VZ3[ii]->Draw();
    cout << VZ3[ii]->GetSum()/1100 << " ";
    //if(ii==0) can[ii]->Print("hitZred.ps(");
    //else if(ii==47) can[ii]->Print("hitZred.ps)");
    //else can[ii]->Print("hitZred.ps");

    sprintf(canname,"R_Vane%d",ii);
    can[ii+48] = new TCanvas(canname,canname,1000,800);
    can[ii+48]->SetFillColor(10);
    can[ii+48]->Divide(1,3);    
    cout << " R ";
    can[ii+48]->cd(1); VR0[ii]->Draw(); 
    cout << VR0[ii]->GetSum()/2000 << " ";
    can[ii+48]->cd(2); VR1[ii]->Draw();
    cout << VR1[ii]->GetSum()/2000 << " ";
    can[ii+48]->cd(3); VR2[ii]->Draw();
    cout << VR2[ii]->GetSum()/2000 << endl;

    //if(ii==0) can[ii+48]->Print("hitRred.ps(");
    //else if(ii==47) can[ii+48]->Print("hitRred.ps)");
    //else can[ii+48]->Print("hitRred.ps");
  }
    */

  /*
  TH2F *waku2 = new TH2F("waku2","",1000,-350,350,1000,-350,350);
  waku2->SetFillColor(10);
  waku2->SetXTitle("X[mm]");
  waku2->SetYTitle("Y[mm]");
  waku2->GetYaxis()->SetLabelSize(0.038);
  waku2->GetYaxis()->SetTitleSize(0.035);
  waku2->GetYaxis()->SetTitleOffset(1.3);
  waku2->Draw();

 TArc *ring = new TArc(0,0,333);
  ring->SetLineColor(4);
  TArc *center = new TArc(0,0,70);
  center->SetLineColor(3);
  //ring->Draw();
  //center->Draw();
  center->SetFillColor(3);
  int SiNum=48;
  double RTPC=290,Thick=0.15;
  double DetOut[48][4],DetIn[48][4]; 
  double pp;
  TLine *Det[48],*Det2[48]; 
  double pi = TMath::Pi();
  //  for(int di=0;di<SiNum;++di){
  for(int di=0;di<48;++di){
    if(SiNum==48)pp=di*7.5*pi/(double)180;
    else if(SiNum==32) pp=di*11.25*pi/(double)180;
    DetOut[di][0]=RTPC*cos(pp)-Thick*sin(pp);
    DetOut[di][1]=RTPC*sin(pp)+Thick*cos(pp);
    DetOut[di][2]=RTPC*cos(pp)+Thick*sin(pp);
    DetOut[di][3]=RTPC*sin(pp)-Thick*cos(pp);
    DetIn[di][0]=(RTPC-220)*cos(pp)-Thick*sin(pp);
    DetIn[di][1]=(RTPC-220)*sin(pp)+Thick*cos(pp);
    DetIn[di][2]=(RTPC-220)*cos(pp)+Thick*sin(pp);
    DetIn[di][3]=(RTPC-220)*sin(pp)-Thick*cos(pp);
    Det[di] = new TLine(DetOut[di][0],DetOut[di][1],DetIn[di][0],DetIn[di][1]);
    Det[di]->SetLineColor(kMagenta+3);
    Det2[di] = new TLine(DetOut[di][2],DetOut[di][3],DetIn[di][2],DetIn[di][3]);
    Det2[di]->SetLineColor(kMagenta+3);
    Det[di]->Draw("same");
    Det2[di]->Draw("same");

  }


  tg->SetMarkerStyle(20);
  tg->Draw("same,p");
  */
}
