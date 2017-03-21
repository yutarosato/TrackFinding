#include "io.hh"

using namespace std;


IO::IO(){};
IO::~IO(){};

void IO::OpenFile(){
  data.open("sample/48v10n/sample9.txt");
  if(!data){
    cout << "Fail to open a file ! " << endl;
    exit(1);
  }
}

void IO::CloseFile(){
  data.close();
}

void IO::OpenFile(int i=0){
  file = new TFile("track.root");
}

void IO::CloseFile(int i=0){
  file->Close();
}

void IO::OpenTree(){
  TD = (TTree*)file->Get("ntupleDecay");
  T = (TTree*)file->Get("ntupleBody");
  TD->ls();
  T->ls();

  //SetBranch for TD
  neventD = TD->GetEntries();
  TD->SetBranchAddress("eventNum",&eventNumD);
  TD->SetBranchAddress("DtEnergy",tEnergyD);
  TD->SetBranchAddress("Dmomv_x",momv_xD);
  TD->SetBranchAddress("Dmomv_y",momv_yD);
  TD->SetBranchAddress("Dmomv_z",momv_zD);
  TD->SetBranchAddress("Dmom_x",mom_xD);
  TD->SetBranchAddress("Dmom_y",mom_yD);
  TD->SetBranchAddress("Dmom_z",mom_zD);
  TD->SetBranchAddress("Dpol_x",pol_xD);
  TD->SetBranchAddress("Dpol_y",pol_yD);
  TD->SetBranchAddress("Dpol_z",pol_zD);
  TD->SetBranchAddress("Dptime",ptimeD);
  TD->SetBranchAddress("Dgtime",gtimeD);
  TD->SetBranchAddress("Dpos_x",pos_xD);
  TD->SetBranchAddress("Dpos_y",pos_yD);
  TD->SetBranchAddress("Dpos_z",pos_zD);

  //SetBranch for T
  T->SetBranchAddress("eventNum",&eventNum);
  T->SetBranchAddress("hitInfo",&hitInfo);
  T->SetBranchAddress("tEnergy",&tEnergy,&btEnergy);
  T->SetBranchAddress("mom_x",&mom_x,&bmom_x);
  T->SetBranchAddress("mom_y",&mom_y,&bmom_y);
  T->SetBranchAddress("mom_z",&mom_z,&bmom_z);
  T->SetBranchAddress("pos_x",&pos_x,&bpos_x);
  T->SetBranchAddress("pos_y",&pos_y,&bpos_y);
  T->SetBranchAddress("pos_z",&pos_z,&bpos_z);
  T->SetBranchAddress("bodyTyp",&bodyTyp,&bbodyTyp);
  T->SetBranchAddress("bodyStatus",&bodyStatus,&bbodyStatus);
  T->SetBranchAddress("gtime",&gtime,&bgtime);
  T->SetBranchAddress("pID",&pID,&bpID);
  T->SetBranchAddress("EachDepE",&EachDepE,&bEachDepE);
  T->SetBranchAddress("kEnergy",&kEnergy,&bkEnergy);
  T->SetBranchAddress("CurrentDepE",&CurrentDepE,&bCurrentDepE);

}

void IO::ClearVector(){
  XX.clear();  YY.clear();  ZZ.clear();
  RR.clear();  PPhi.clear(); IId.clear();
  TThr.clear(); TTime.clear();
  EvStr.clear();
}

int IO::GetEventNumber(){
  return T->GetEntries();
}

void IO::ReadOneEvent(int i){
  // i=652;
  //i=321;
  //i=721;
  //  i=185;
  //i=4704;
  TD->GetEntry(i);
  T->GetEntry(i);


  cout << i << endl;
  //cout << "---" << eventNumD << " " << eventNum << " " << tEnergyD[1]<<endl;

  int vecsize = pos_x->size();
  //cout << vecsize << " " << tEnergyD[1] << " " << gtimeD[0] << endl;
  //cout << pos_xD[0] << " " << pos_yD[0] << " " << pos_zD[0] << endl;
  //cout << mom_xD[1] << " " << mom_yD[1] << " " << mom_zD[1] << endl;
  int tflg=0;
  double stime=0,etime=0,hnum=0,onum=0;
  double esx=0,esy=0,esz=0,esvx=0,esvy=0,esvz=0;
  double eex=0,eey=0,eez=0,eevx=0,eevy=0,eevz=0;


  for(int jj=0; jj<vecsize; jj++){
    double tmp_x=pos_x->at(jj);
    double tmp_y=pos_y->at(jj);
    double tmp_z=pos_z->at(jj);
    double tmp_vx=mom_x->at(jj);
    double tmp_vy=mom_y->at(jj);
    double tmp_vz=mom_z->at(jj);
    double tmp_energy=tEnergy->at(jj);
    double tmp_time=gtime->at(jj);
    int tmp_BTyp=bodyTyp->at(jj);
    int pre_BTyp=0;
    if(jj>0) pre_BTyp=bodyTyp->at(jj-1);
    int particleID=pID->at(jj); //1=>e- 2=>e+ 3=>gamma 
    int statBody=bodyStatus->at(jj);
    double tmp_dep=EachDepE->at(jj); //deposit
    double tmp_kE=kEnergy->at(jj);
    double tmp_cdep=CurrentDepE->at(jj); // total deposit until "i"



    if(1 && 
       tmp_dep>0 &&
       statBody==0 &&
       tmp_BTyp>100 &&
       tmp_BTyp<1000
       ){
      
      onum++;
      if(particleID==2){
	hnum++;
	if(tflg==0){
	  stime=tmp_time;
	  esx = tmp_x; esy = tmp_y; esz = tmp_z;
	  esvx = tmp_vx; esvy = tmp_vy; esvz = tmp_vz;

	  tflg++;
	}
	else{
	  etime=tmp_time;
	  eex = tmp_x; eey = tmp_y; eez = tmp_z;
	  eevx = tmp_vx; eevy = tmp_vy; eevz = tmp_vz;
	}
	//cout<<tmp_time << " " <<tmp_x << " " <<tmp_y << " " <<tmp_z <<endl;
      }
      
      //added 110915
      //if(etime-stime >15) break;
      
      int id=i;
      int th=1;
      if(tEnergyD[1]<=150) th=-1;
      double tmp_r = TMath::Sqrt(tmp_x*tmp_x+tmp_y*tmp_y);
      double cosphi=tmp_x/tmp_r;
      double sinphi=tmp_y/tmp_r;
      double Phi = TMath::ACos(cosphi);
      if(sinphi<0) Phi = 2*TMath::Pi()-Phi;
    
      //cout << th << " " << tmp_time << " " << Phi << " " << tmp_z << endl;
      //      cout << tEnergyD[1] << endl;
      XX.push_back(tmp_x); YY.push_back(tmp_y); ZZ.push_back(tmp_z);
      RR.push_back(tmp_r); PPhi.push_back(Phi), IId.push_back(id);
      TThr.push_back(th); TTime.push_back(tmp_time);
    }
  }

  //cout << "pos :" << esx << " " << esy << " " << esz << endl;
  //cout << "mom :" << esvx << " " << esvy << " " << esvz << endl;
  //cout << "aa " << hnum << " " << stime << " " << etime << endl;
  double aho=mom_xD[1]*mom_xD[1]+mom_yD[1]*mom_yD[1];//+mom_zD[1]*mom_zD[1];
  aho = TMath::Sqrt(aho);
  double ahox= pos_xD[0]+ mom_xD[1]/aho;
  double ahoy=pos_yD[0]+mom_yD[1]/aho;
  //cout << TMath::Sqrt(ahox*ahox+ahoy*ahoy) << " "
  //     << TMath::Sqrt(pos_xD[0]*pos_xD[0]+pos_yD[0]*pos_yD[0]+pos_zD[0]*pos_zD[0]) << endl;

  if(onum<1){
    esx=-1000; esy=-1000; esz=-1000;
    esvx=-1000; esvy=-1000; esvz=-1000;
    eex=-1000; eey=-1000; eez=-1000;
    eevx=-1000; eevy=-1000; eevz=-1000;
    stime=-1000; etime=-1000;
  }

  EvStr.push_back(i);           //0
  EvStr.push_back(hnum);        //1
  EvStr.push_back(onum);        //2
  EvStr.push_back(tEnergyD[1]); //3
  EvStr.push_back(gtimeD[0]);   //4
  EvStr.push_back(pos_xD[0]);   //5
  EvStr.push_back(pos_yD[0]);   //6
  EvStr.push_back(pos_zD[0]);
  EvStr.push_back(mom_xD[1]);
  EvStr.push_back(mom_yD[1]);
  EvStr.push_back(mom_zD[1]);   //10
  EvStr.push_back(stime);       //11
  EvStr.push_back(etime);
  EvStr.push_back(esx);
  EvStr.push_back(esy);
  EvStr.push_back(esz);        //15
  EvStr.push_back(esvx);       //16
  EvStr.push_back(esvy);
  EvStr.push_back(esvz);
  EvStr.push_back(eex);
  EvStr.push_back(eey);        //20
  EvStr.push_back(eez);        //21
  EvStr.push_back(eevx);
  EvStr.push_back(eevy);
  EvStr.push_back(eevz);


}

void IO::GetEvtStrParam(vector<double> &EvSt){
  EvSt = EvStr;
}

void IO::GetParam(vector<double> &X, vector<double> &Y, vector<double> &Z,
		 vector<double> &R, vector<double> &Phi,
		  vector<int> &Id, vector<int> &Thr, vector<double> &Time){
  X = XX; Y=YY; Z=ZZ; R=RR; Phi=PPhi; Id=IId; Thr=TThr; Time=TTime;

}

void IO::GetDecayPoint(double &x,double &y,double &z){
  x = pos_xD[0];   y = pos_yD[0];   z = pos_zD[0]; 
  //cout << "aa " << x << " " << y << " " << z << endl;
}

double IO::GetPositronInitialEnergy(){
  return tEnergyD[1];
}


void IO::ReadData(vector<double> &X, vector<double> &Y, vector<double> &Z,
		 vector<double> &R, vector<double> &Phi,
		 vector<int> &Id, vector<int> &Thr){
  int id,th;
  double x,y,z,r,phi,time,cE,iE;
  while(1){
    data >> id >> th >> time >> x >> y >> z >> r >> phi >> cE >> iE;
    if(data.eof())break;

    phi = TMath::ATan(y/x);
    if(x<0) phi=phi+TMath::Pi();
    else if(x>=0&&y<0) phi = phi+TMath::Pi()*2;

    r = x*x+y*y;
    r = TMath::Sqrt(r);

    X.push_back(x); Y.push_back(y); Z.push_back(z);
    R.push_back(r); Phi.push_back(phi); Id.push_back(id);
    Thr.push_back(th);
    //X[index]=X[index]*360/6.28;
  }
}

void IO::DigitizeVane(int vane, vector<double> &Phi, vector<double> &VaneNum){
  int ind=Phi.size();
  double tmpVane;
  int VNum;
  for(int i=0;i<ind;i++){
    tmpVane=(Phi[i]/(2*TMath::Pi()/vane));
    VNum = (int)tmpVane;
    if(tmpVane-VNum >=0.5){
      VNum++;
      if(VNum==vane) VNum=0;
    }
    VaneNum.push_back(VNum);
  }
}

void IO::AddGhost(int vane, vector<double> &R, vector<double> &Z, 
		  vector<double> &VaneNum, vector<int> &UnitNum,
		  vector<double> &AfterR, vector<double> &AfterZ){
  // int index=R.size();
  int unitn;
  int stripR=0,stripZ=0;
  int VaneUnit[vane][16],EvID[vane][100];
  int VaneEvN[vane];
  vector<int> StripNumX, StripNumY;
  bool Coinc[vane];
  //cout << "index " << index <<endl;
  for(int i=0;i<vane;i++){
    VaneEvN[i]=0;
    for(int j=0;j<16;j++){
      VaneUnit[i][j]=0;
    }
    for(int k=0;k<100;k++){
      EvID[i][k]=0;
    }
  }
  for(int i=0; i<index;i++){
    unitn = GetUnit(R[i],Z[i]);
    UnitNum.push_back(unitn);
    GetStrips(unitn,R[i],Z[i],stripR,stripZ);
    StripNumX.push_back(stripR);
    StripNumY.push_back(stripZ);
    //cout << stripR << " " << stripZ << endl;
    VaneUnit[(int)VaneNum[i]][unitn]++;
    EvID[(int)VaneNum[i]][VaneEvN[(int)VaneNum[i]]]=i;
    VaneEvN[(int)VaneNum[i]]++;
  }
  int index2=StripNumX.size();

//cout << "index2  " << index2 <<endl;
  for(int i=0;i<vane;i++){
    Coinc[i]=CheckCoinci(i,VaneUnit);
    if(Coinc[i]){
      AddPoints(i,StripNumX,StripNumY,VaneEvN,EvID,UnitNum,VaneNum);
    }
    //    cout << i << " " << Coinc[i] << endl;
  }
  int index3=StripNumX.size();
  //cout << "index3  " << index3 <<endl;

  ConvStripToReal(UnitNum, StripNumX, StripNumY, AfterR, AfterZ);

  int index4=AfterR.size();


  //cout << "index4  " << index4 <<endl;

}

int IO::GetUnit(double r, double z){
  //R[75,295]
  //Z[-200,200]

  int x,y,unum;
  x = ((int)r-75)/55; 
  y = -1*((int)z-200)/100;
  unum = x+4*y;
  return unum;
}

void IO::GetStrips(int uni, double r, double z, int &sR, int &sZ){
  //R[75,295]
  //Z[-200,200]
  double tmpx=0.,tmpy=0.;
  tmpx=r-75;
  tmpy=-1*(z-200);
  if(uni%4==1) tmpx = tmpx-55;
  else if(uni%4==2) tmpx = tmpx-55*2;
  else if(uni%4==3) tmpx = tmpx-55*3;
  if(uni/4==1) tmpy = tmpy-100;
  else if(uni/4==2) tmpy = tmpy-200;
  else if(uni/4==3) tmpy = tmpy-300;

  sR = (int)(tmpx/0.2);
  sZ = (int)(tmpy/0.2);
  //  cout << r << " " << tmpx << " " << z << " " << tmpy << endl;
  //if(sR>274 || sZ>499) cout << sR << " " << sZ << endl;
}

bool IO::CheckCoinci(int vnum, int VU[][16]){
  for(int i=0;i<16;i++){
    if(VU[vnum][i]>1) return true;
  }
  return false;
}

void IO::AddPoints(int vaneN, vector<int> &R, vector<int> &Z,
		   int ID[], int Ev[][100], vector<int> &Uni, vector<double> &VN){
  int tmpUnit[16];
  for(int i=0; i<16;i++) tmpUnit[i]=0;

  for(int i=0; i<ID[vaneN];i++){
    tmpUnit[Uni[Ev[vaneN][i]]]++;
    //cout << Ev[vaneN][i] << " " << Uni[Ev[vaneN][i]] << " "
    //	 <<R[Ev[vaneN][i]] << " " << Z[Ev[vaneN][i]] << endl;
  }
  int tmpX[50],tmpY[50];
  for(int i=0; i<50;i++){
    tmpX[i]=0; tmpY[i]=0;
  }
  int tmpindex=0;
  for(int i=0; i<16;i++){
    if(tmpUnit[i]>1){
      tmpindex=0;
      for(int j=0; j<ID[vaneN];j++){
	if(Uni[Ev[vaneN][j]] ==i){
	  tmpX[tmpindex]=R[Ev[vaneN][j]];
	  tmpY[tmpindex]=Z[Ev[vaneN][j]];
	  tmpindex++;
	}
      }
      for(int j=0;j<tmpindex;j++){
	for(int k=0;k<tmpindex;k++){
	  if(j!=k){
	    VN.push_back(vaneN);
	    Uni.push_back(i);
	    R.push_back(tmpX[j]);
	    Z.push_back(tmpY[k]);
	    //cout << i << " " << tmpX[j] << " " << tmpY[k] << endl;
	  }
	}
      }
    }
  } 
}

void IO::ConvStripToReal(vector<int> &UnitNum, 
			 vector<int> &StripNumX, vector<int> &StripNumY, 
			 vector<double> &AfterR, vector<double> &AfterZ){
  int evnum=UnitNum.size();
  double tmpX,tmpY;
  for(int i=0; i<evnum;i++){
    tmpX = (double)StripNumX[i]*0.2+0.1;
    tmpY = (double)StripNumY[i]*0.2+0.1;
    if(UnitNum[i]%4==1) tmpX=tmpX+55;
    else if(UnitNum[i]%4==2) tmpX=tmpX+55*2;
    else if(UnitNum[i]%4==3) tmpX=tmpX+55*3;
    tmpX = tmpX+75;
    if(UnitNum[i]/4==1)tmpY=tmpY+100;
    else if(UnitNum[i]/4==2) tmpY=tmpY+200;
    else if(UnitNum[i]/4==3) tmpY=tmpY+300;
    tmpY = -1*tmpY+200;

    AfterR.push_back(tmpX);
    AfterZ.push_back(tmpY);
  }

}
