#include "root.hh"
#include "io.hh"
#include "algo.hh"
#include "Display.hh"
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

//-Class
IO *io;
Algo *algo;
Display *disp;

//-Function
void Initial();
void ClearVector();

//-Variables
double Rawx[1000],Rawy[1000];
double HoughR[1000*1000],HoughT[1000*1000];

//int vane=16;
//int vane=32;
//int vane=48;
//int vane=24;
int vane=40;

vector<double> X, Y, Z, R, Phi, VaneNum, AfterR, AfterZ, AfterX, AfterY;
vector<double> dX, dY, dZ, redVaneNum, Time;
vector<int> Id, Thr, UnitNum;

//added in ver.0.12
vector<double> Zcp, Phicp;

//added in event structure
vector<double> ESt;
double ETheta=-1000;

//-Main
int main(int argc, char **argv){

  //--Initialize
  disp = new Display();
  io = new IO();
  algo = new Algo();

  disp->ROOT_Init(argc,argv);
  Initial();

  disp->MakeCanvas(0,4,3);
  disp->MakeFrame(0,"frame","#phi-Z");
  disp->MakeFrame(1,"frame2","X-Y",-350,350,-350,350);
  disp->MakeFrame(2,"frame3","VaneNum-Z",0,vane,-250,250);

  disp->MakeOrbit();

  disp->frame[0]->SetXTitle("#phi [rad]");
  disp->frame[0]->SetYTitle("Z [mm]");
  disp->frame[0]->GetXaxis()->SetTitleSize(0.045);
  disp->frame[0]->GetYaxis()->SetTitleSize(0.045);
  disp->frame[0]->GetXaxis()->SetLabelSize(0.04);
  disp->frame[0]->GetYaxis()->SetLabelSize(0.035);
  disp->frame[0]->GetYaxis()->SetTitleOffset(1.1);

  disp->frame[1]->SetXTitle("X [mm]");
  disp->frame[1]->SetYTitle("Y [mm]");
  disp->frame[1]->GetXaxis()->SetTitleSize(0.045);
  disp->frame[1]->GetYaxis()->SetTitleSize(0.045);
  disp->frame[1]->GetXaxis()->SetLabelSize(0.04);
  disp->frame[1]->GetYaxis()->SetLabelSize(0.035);
  disp->frame[1]->GetYaxis()->SetTitleOffset(1.1);

  disp->frame[2]->SetXTitle("Vane #");
  disp->frame[2]->SetYTitle("Z [mm]");
  disp->frame[2]->GetXaxis()->SetTitleSize(0.045);
  disp->frame[2]->GetYaxis()->SetTitleSize(0.045);
  disp->frame[2]->GetXaxis()->SetLabelSize(0.04);
  disp->frame[2]->GetYaxis()->SetLabelSize(0.035);
  disp->frame[2]->GetYaxis()->SetTitleOffset(1.1);

  //added in event structure
  disp->MakeCanvas(1,2,2);
  disp->can[1]->Draw();
  disp->can[1]->cd(1);
  disp->frame[1]->Draw();
  disp->ring->Draw("same");
  disp->center->Draw("same");
  disp->can[1]->cd(2);
  disp->frame[0]->Draw();
  disp->can[1]->cd(3);
  disp->frame[2]->Draw();
  // --added in event structure

  disp->can[0]->Draw();
  disp->can[0]->cd(2);//
  disp->frame[0]->Draw();
  disp->can[0]->cd(4);
  disp->frame[0]->Draw();
  disp->can[0]->cd(6);
  disp->frame[0]->Draw();
  disp->can[0]->cd(1);
  disp->frame[1]->Draw();
  disp->ring->Draw("same");
  disp->center->Draw("same");

  disp->can[0]->cd(7);
  disp->frame[2]->Draw();
  disp->can[0]->cd(8);
  disp->frame[2]->Draw();
  disp->can[0]->cd(9);
  disp->frame[2]->Draw();
  disp->can[0]->cd(10);
  disp->frame[2]->Draw();


  //disp->can[0]->Print("EDisp.ps(");
  //    disp->can[1]->Print("EDisp.ps(");

  io->OpenFile(0);
  io->OpenTree();
  int delflag=0,delflag2=0,delflag3=0,delflag4=0;

  int OKcount=0,OKtrue=0;
  int counter0=0,counter1=0, counter2=0,counter3=0;
  int truecnt=0, truecnt2=0, truecnt3=0,truecnt4=0;
  int temp[500],temp2[500];
  int tempnum=0,tempnum2=0;
  int trueflag=0, truefirstpointflag=0;
  int evnum = io->GetEventNumber();

  double dpos[3];

  TLatex *texNum,*texEne,*texDet;


  for(int i=0;i<evnum;i++){ //loop
  //  for(int i=0;i<evnum;i=i+500){ //loop
  // for(int i=0;i<30;i++){ //loop
    ETheta=-1000;
    counter0++;
    ClearVector();
    io->ClearVector();
    if(i>0){
      texNum->Delete();
      texEne->Delete();
    }

    if(i>0&&delflag==1){
      disp->graph[0]->Delete();
      disp->graph[1]->Delete();
      disp->graph[3]->Delete();
      disp->his2D[0]->Delete();
      disp->pol1[0]->Delete();
      disp->his[0]->Delete();
      if(delflag2==1){
	disp->graph[2]->Delete();
	disp->graph[4]->Delete();
	disp->graph[5]->Delete();
	if(delflag3==1){
	  disp->graph[6]->Delete();
	  disp->graph[7]->Delete();
	  if(delflag4==1){
	    texDet->Delete();
	  }
	}
      }

      delflag=0; delflag2=0;delflag3=0; delflag4=0;
    }
    
    trueflag=0;
    io->ReadOneEvent(i);
    io->GetParam(X,Y,Z,R,Phi,Id,Thr,Time);
    io->GetDecayPoint(dpos[0],dpos[1],dpos[2]);
    dX.push_back(dpos[0]);
    dY.push_back(dpos[1]);
    dZ.push_back(dpos[2]);
    double dE = io->GetPositronInitialEnergy();

   ofstream fop("newocp.dat",ios::app);
    for(int ii=0; ii<Time.size();ii++){
      fop << i << " " << Time[ii] << " " << R[ii] << " " << Phi[ii] 
          << " " << X[ii] << " " << Y[ii] << " " << Z[ii] << endl;
    }
    fop.close();

    /*
    ofstream fop("ocp.dat",ios::app);
    for(int ii=0; ii<Time.size();ii++){
      fop << i << " " << Time[ii] << endl;
    }
    fop.close();
    */

    char texNumc[40],texEnec[40];
    sprintf(texNumc,"Event Number : %d",i);
    sprintf(texEnec,"Positron Energy : %3.1f [MeV]",dE);
    texNum = new TLatex(0.05,0.8,texNumc);
    texEne = new TLatex(0.05,0.6,texEnec);
    texNum->SetTextSize(0.07);
    texEne->SetTextSize(0.07);
    if(dE>150){
      texEne->SetTextColor(2);
    }
    else{
      texEne->SetTextColor(4);
    }

    if(X.size()>0 && dE>150){ 
      truecnt++;
    }
    if(X.size()>3){// 

      //ofstream tmpout("Zfirst.dat",ios::app);
      //tmpout << Z[0] << " " << dE << endl;
      //tmpout.close();

      counter1++;
      if(Thr[0]>0) truecnt3++;
      delflag=1;
      disp->MakeGraph(0,Phi,Z);
      disp->MakeGraph(1,X,Y);
      disp->MakeGraph(3,dX,dY);
      disp->can[0]->cd(2);//
      disp->graph[0]->Draw("p,same");
      disp->can[0]->cd(1);
      disp->graph[1]->Draw("p,same");
      disp->graph[3]->SetMarkerColor(2);
      disp->graph[3]->Draw("p,same");
   
      disp->can[1]->cd(2);//
      disp->graph[0]->Draw("p,same");
      disp->can[1]->cd(1);
      disp->graph[1]->Draw("p,same");
 

      algo->HoughTransform(Phi,Z,AfterX,AfterY);
     
      disp->Make2DHist(0,AfterX,AfterY);

      disp->can[0]->cd(3);
      disp->his2D[0]->SetXTitle("#theta [#circ]");
      disp->his2D[0]->SetYTitle("r");
      disp->his2D[0]->GetXaxis()->SetTitleSize(0.045);
      disp->his2D[0]->GetXaxis()->SetTitleOffset(1.1);
      disp->his2D[0]->GetYaxis()->SetTitleOffset(1.3);
      disp->his2D[0]->Draw("colz");
      double par0,par1;
      algo->OneHoughFit(disp->his2D[0],par0,par1);
      
      disp->can[0]->cd(4);
      disp->graph[0]->Draw("p,same");
      disp->MakePol1(0,par0,par1*180/TMath::Pi());
      disp->pol1[0]->Draw("same");
      
      vector<double> Resi;
      algo->GetPol1Residual(disp->pol1[0],Phi,Z,Resi);
      disp->MakeHist(0,Resi);
      disp->can[0]->cd(5);
      disp->his[0]->SetXTitle("residual[mm]");
      disp->his[0]->Draw();
      
      vector<double> redZ,redPhi;
      for(int ic=0; ic<(int)Resi.size();ic++){
	if(fabs(Resi[ic]-disp->his[0]->GetMean())<=3*(disp->his[0]->GetRMS())){
	  redZ.push_back(Z[ic]);
	  redPhi.push_back(Phi[ic]);
	}
      }
      double angleflag=0;
      angleflag = TMath::ATan(par1*180/TMath::Pi())*180/TMath::Pi();


      ////added in ver.0.12

      if(fabs(angleflag)>=89.9){
	std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
	Resi.clear();
	algo->GetPol1ResidualY(disp->pol1[0],Phi,Z,Resi);
	disp->his[0]->Delete();
	disp->MakeHist(0,Resi);
	disp->can[0]->cd(5);
	disp->his[0]->SetXTitle("residual[degree]");
	disp->his[0]->Draw();

	for(int ic=0; ic<(int)Resi.size();ic++){
	  if(fabs(Resi[ic]-disp->his[0]->GetMean())<=3*(disp->his[0]->GetRMS())){
	  }
	  else{
	    Zcp.push_back(Z[ic]);
	    Phicp.push_back(Phi[ic]);
	  }
	}

	
	AfterX.clear(); AfterY.clear();
	redZ.clear(); redPhi.clear();
	Resi.clear();
	
	disp->his2D[0]->Delete();
	disp->pol1[0]->Delete();
	disp->his[0]->Delete();
	algo->HoughTransform(Phicp,Zcp,AfterX,AfterY);

	disp->Make2DHist(0,AfterX,AfterY);

	disp->can[0]->cd(3);
	disp->his2D[0]->SetXTitle("#theta [#circ]");
	disp->his2D[0]->SetYTitle("r");
	disp->his2D[0]->GetXaxis()->SetTitleSize(0.045);
	disp->his2D[0]->GetXaxis()->SetTitleOffset(1.1);
	disp->his2D[0]->GetYaxis()->SetTitleOffset(1.3);
	disp->his2D[0]->Draw("colz");
 
	algo->OneHoughFit(disp->his2D[0],par0,par1);
	//disp->graph[0]->Delete();
	//disp->MakeGraph(0,Phicp,Zcp);
	disp->can[0]->cd(4);
	disp->graph[0]->Draw("p,same");
	disp->MakePol1(0,par0,par1*180/TMath::Pi());
	disp->pol1[0]->Draw("same");

	algo->GetPol1Residual(disp->pol1[0],Phi,Z,Resi);
	disp->MakeHist(0,Resi);
	disp->can[0]->cd(5);
	disp->his[0]->SetXTitle("residual[mm]");
	disp->his[0]->Draw();
	
	for(int ic=0; ic<(int)Resi.size();ic++){
	  if(fabs(Resi[ic]-disp->his[0]->GetMean())<=3*(disp->his[0]->GetRMS())){
	    redZ.push_back(Z[ic]);
	    redPhi.push_back(Phi[ic]);
	  }
	}
	
	angleflag = TMath::ATan(par1*180/TMath::Pi())*180/TMath::Pi();
	//cout << angleflag << endl;
	
      }
      ////ver.0.12 fin

      ETheta = angleflag;

      if(redZ.size()>0 && fabs(angleflag)<89.9){
	counter2++;
	delflag2=1;
	if(Thr[0]>0){
	  truecnt2++;
	}
	disp->MakeGraph(2,redPhi,redZ);
	disp->can[0]->cd(6);
	disp->graph[2]->Draw("p,same");

	//Digitize
	io->DigitizeVane(vane,Phi,VaneNum);
	io->DigitizeVane(vane,redPhi,redVaneNum);
	disp->MakeGraph(4,VaneNum,Z);
	disp->MakeGraph(5,redVaneNum,redZ);
	disp->can[0]->cd(7);
	disp->graph[4]->Draw("p,same");
	disp->can[0]->cd(8);
	disp->graph[5]->Draw("p,same");


	//Sorting
	multimap<double,double> tMap;
	for(int ii=0; ii<(int)redVaneNum.size(); ii++){
	  tMap.insert(make_pair(redVaneNum[ii],redZ[ii]));
	}
	redVaneNum.clear(); redZ.clear();
	multimap<double,double>::iterator it = tMap.begin();
	while(it != tMap.end()){
	  redVaneNum.push_back((*it).first);
	  redZ.push_back((*it).second);
	  it++;
	}
	/*
	cout << "--------------" << endl;
	for(int ii=0;ii<(int)redVaneNum.size();ii++){
	  cout << redVaneNum[ii] << " " << redZ[ii] << endl;
	}
	*/
	//cout << "--------------" << endl;
	//Search overlap
	//-get divided groups
	int group[100],groupS[100],groupE[100];
	int groupnum=0;
	int groupflag=0;

	for(int ii=0;ii<100;ii++){
	  group[ii]=0; groupS[ii]=0; groupE[ii]=0;
	}

	for(int ii=(int)redVaneNum[0]; ii<=(int)redVaneNum[redVaneNum.size()-1];ii++){
	  if(tMap.count(ii)!=1){
	    if(groupflag==0 && ii>(int)redVaneNum[0]){
	      groupE[groupnum]=ii-1;
	    }
	    groupflag++;
	  }
	  else if(tMap.count(ii)==1){
	    if(groupflag>0){
	      if(group[groupnum]>0)
		groupnum++;
	      if(group[groupnum]==0){
		groupS[groupnum]=ii;
		group[groupnum]++;
	      }
	    }
	    else{
	      if(group[groupnum]==0){
		groupS[groupnum]=ii;
	      }
	      group[groupnum]++;
	      if(ii== (int)redVaneNum[redVaneNum.size()-1]){
		groupE[groupnum]=ii;
	      }
	    }
	    groupflag=0;
	  }
	  else{
	    group[groupnum]++;
	  }
	}

	//-get group with maximum number
	int tmpmax=0,maxgroup=0;
	if(group[0]>0){
	  for(int ii=0; ii<groupnum+1;ii++){
	    if(tmpmax<group[ii]){
	      tmpmax=group[ii];
	      maxgroup=ii;
	    }
	  }
	  //	  cout << maxgroup << " " << tmpmax << endl;
	}

	//-cut group with small number & clustering
	vector<double> clsVaneNum,clsZ;
	if(tmpmax>1){//def(v012)->3 changed in ver012_2
	  counter3++;
	  if(Thr[0]>0) truecnt4++;
	  delflag3=1;
	  for(int ii=0; ii<(int)redVaneNum.size(); ii++){
	    if(redVaneNum[ii]>=groupS[maxgroup]&&
	       redVaneNum[ii]<=groupE[maxgroup]){
	      clsVaneNum.push_back(redVaneNum[ii]);
	      clsZ.push_back(redZ[ii]);
	    }
	  }
	  disp->MakeGraph(6,clsVaneNum,clsZ);
	  disp->can[0]->cd(9);
	  disp->graph[6]->Draw("p,same");

	  //-steering
	  //-back
	  double hashi,mae,hashiY,maeY,exval;
	  double subt,minsub,minsubY;
	  int overii[1000*VaneNum.size()];
	  int overiicnt=0;
	  for(int ii=0;ii<(int)VaneNum.size()*1000;ii++){
	    overii[ii]=VaneNum.size();
	  }
	  int mplflag=0; //added in v012_3
	  while(1){
	    minsub=1000;
	    hashi = clsVaneNum[clsVaneNum.size()-1];
	    hashiY = clsZ[clsVaneNum.size()-1];
	    maeY = clsZ[clsVaneNum.size()-2];
	    exval = hashiY + hashiY - maeY;
	    //cout << hashi << endl;
	    if(hashi==vane-1) hashi=hashi-vane;
	    for(int ii=0;ii<(int)VaneNum.size();ii++){
	      if(VaneNum[ii] == hashi+1){
		subt = TMath::Abs(exval-Z[ii]);
		if(subt<minsub){
		  minsub=subt;
		  minsubY=Z[ii];
		  for(int jj=0; jj<overiicnt;jj++){
		    if(overii[jj]==ii)minsub=1000;
		  }
		  overii[overiicnt]=ii;
		  overiicnt++;
		}
	      }
	    }
	    //cout << minsub << endl;
	    if(minsub>10){
	      // for bug file
	      int tmpfl=0;
	      mplflag++; //added in v012_3
	      for(int ii=0;ii<(int)VaneNum.size();ii++){
		//tmpfl=0;
		if(VaneNum[ii] == hashi+2){
		  subt = TMath::Abs(exval-Z[ii]);
		  if(subt<minsub){
		    tmpfl=1;
		    minsub=subt;
		    minsubY=Z[ii];
		    for(int jj=0; jj<overiicnt;jj++){
		      if(overii[jj]==ii)minsub=100;
		    }
		    overii[overiicnt]=ii;
		    overiicnt++;
		  }
		}
	      } // ---for bug file
	      if(tmpfl==0 || mplflag>300){ //added in v012_3
		break;
	      }
	      else{//added in v012_3
		hashi=hashi+1;
		//cout << subt << endl;
	      }//--added in v012_3
	    }
	    clsVaneNum.push_back(hashi+1);
	    clsZ.push_back(minsubY);
	  }
	  for(int ii=0;ii<(int)VaneNum.size()*1000;ii++){
	    overii[ii]=VaneNum.size();
	  }
	  overiicnt=0;
	  //-front
	  while(1){
	    minsub=1000;
	    hashi = clsVaneNum[0];
	    hashiY = clsZ[0];
	    maeY = clsZ[1];
	    exval = hashiY + hashiY - maeY;
	    if(hashi==0) hashi=hashi+vane;
	    for(int ii=0;ii<(int)VaneNum.size();ii++){
	      if(VaneNum[ii]==hashi-1){
		subt=TMath::Abs(exval-Z[ii]);
		if(subt<minsub){
		  minsub=subt;
		  minsubY=Z[ii];
		  for(int jj=0; jj<overiicnt;jj++){
		    if(overii[jj]==ii)minsub=1000;
		  }
		  overii[overiicnt]=ii;
		  overiicnt++;
		}
	      }
	    }
	    if(minsub>10){
	      // for bug file
	      int tmpfl=0;
	      for(int ii=0;ii<(int)VaneNum.size();ii++){
		tmpfl=0; //added in v012_3
		if(VaneNum[ii]==hashi-2){
		  subt=TMath::Abs(exval-Z[ii]);
		  if(subt<minsub){
		    tmpfl=1;
		    minsub=subt;
		    minsubY=Z[ii];
		    for(int jj=0; jj<overiicnt;jj++){
		      if(overii[jj]==ii)minsub=100;//changed in ver012_3
		    }
		    overii[overiicnt]=ii;
		    overiicnt++;
		  }
		}
	      } // ---for bug file

	      if(tmpfl==0){
		break;
	      }
	      else{//added in v012_3
		hashi=hashi-1;
		//cout << subt << endl;
	      }//--added in v012_3
	    }
	    clsVaneNum.insert(clsVaneNum.begin(),hashi-1);
	    clsZ.insert(clsZ.begin(),minsubY);
	  }

	  disp->MakeGraph(7,clsVaneNum,clsZ);
	  disp->can[0]->cd(10);
	  disp->graph[7]->Draw("p,same");
	  disp->can[1]->cd(3);
	  disp->graph[7]->Draw("p,same");

	  double mintime[3],mintimei[3],minZ[3],minVane[3];
	  for(int ii=0;ii<3;ii++){
	    mintime[ii]=1000000;
	    mintimei[ii]=0;
	    minZ[ii]=0;
	    minVane[ii]=0;
	  }
		
	  //GetTrue3hits   //removed in v012_3
	  for(int ii=0;ii<(int)Z.size();ii++){
	    if(mintime[0]>Time[ii]){
	      mintime[0]=Time[ii];
	      mintimei[0]=ii;
	    }
	  }
	  for(int ii=0;ii<(int)Z.size();ii++){
	    if(mintime[1]>Time[ii] && Time[ii]>mintime[0]){
	      mintime[1]=Time[ii];
	      mintimei[1]=ii;
	    }
	  }
	  for(int ii=0;ii<(int)Z.size();ii++){
	    if(mintime[2]>Time[ii] && Time[ii]>mintime[1]){
	      mintime[2]=Time[ii];
	      mintimei[2]=ii;
	    }
	  }
	  for(int ii=0;ii<3;ii++){
	    minZ[ii]= Z[mintimei[ii]];
	    minVane[ii]= VaneNum[mintimei[ii]];
	  }
	  /*
	  cout << minZ[0] << " " << minZ[1] << " " << minZ[2] <<endl;
	  cout << Z[0] << " " << Z[1] << " " << Z[2] 
	       << " " << Z[3] << " " << Z[4] << endl;

	  for(int ii=0; ii<(int)clsZ.size();ii++){
	    cout << clsZ[ii] << endl;
	  }
	  */

	  //Compare  //changed in v012_3
	  trueflag=0;
	  truefirstpointflag=0;
	  /*
	  for(int ii=0;ii<(int)clsZ.size();ii++){
	    for(int jj=0;jj<3;jj++){
	      if(clsVaneNum[ii]==minVane[jj] && clsZ[ii]==minZ[jj]){
		trueflag++;
	      }
	    }
	  }
	  */
	  for(int ii=0;ii<(int)clsZ.size();ii++){
	    for(int jj=0;jj<4;jj++){
	      if(clsVaneNum[ii]==VaneNum[jj] && clsZ[ii]==Z[jj]){
		trueflag++;
	      }
	    }
	    if(clsVaneNum[ii]==VaneNum[0] &&clsZ[ii]==Z[0])
	      truefirstpointflag=1;
	  }  
	  

	  //cout << trueflag << endl;
	  if(trueflag>=3 && truefirstpointflag==1){//changed in v012_3
	    //ofstream fobt("true_num.dat",ios::app);
	    //fobt << i << " " << dE << endl;
	    //fobt.close();
	    //cout << "----OK----" << endl;	  
	    OKcount++;
	    texDet = new TLatex(0.2,0.3,"Detected!");
	    disp->can[0]->cd(12);
	    texDet->SetTextSize(0.1);
	    texDet->Draw("same");
	    delflag4=1;
	    if(Thr[0]>0) OKtrue++;
	  }
	
	}	//-End -cut group with small number & clustering



	//	cout << "--------------" << endl;

	
      }
      else{
	//	if(Thr[0]>0) cout << "----------" << endl;
      }
      disp->can[0]->cd(12);
      texNum->Draw("same");
      texEne->Draw("same");
      disp->can[0]->Update();
      //disp->can[0]->Print("EDisp.ps");
      //if(dE>150 && trueflag<3) disp->can[0]->Print("EDisp.ps");
      //int tmp; scanf("%d",&tmp);
    }
    else{
      
      if(dE>150 && X.size()>0){
	temp[tempnum]=i;
	tempnum++;
	disp->can[0]->cd(1);
	disp->MakeGraph(3,dX,dY);
	disp->graph[3]->SetMarkerColor(2);
	disp->graph[3]->Draw("p,same");
	disp->MakeGraph(0,Phi,Z);
	disp->MakeGraph(1,X,Y);
	disp->graph[1]->Draw("p,same");
	disp->can[0]->cd(2);//
	disp->graph[0]->Draw("p,same");
	disp->can[0]->cd(12);
	texNum->Draw("same");
	texEne->Draw("same");

	disp->can[0]->Update();


	//disp->can[0]->Print("EDisp.ps");
	//int tmp; scanf("%d",&tmp);
	disp->graph[0]->Delete();
	disp->graph[1]->Delete();
	disp->graph[3]->Delete();

      }
      //else if(dE>150 && X.size()<=0){
      else{
	temp2[tempnum2]=i;
	tempnum2++;
	disp->can[0]->cd(1);
	disp->MakeGraph(3,dX,dY);
	disp->graph[3]->SetMarkerColor(2);
	disp->graph[3]->Draw("p,same");
	disp->can[0]->cd(12);
	texNum->Draw("same");
	texEne->Draw("same");
	disp->can[0]->Update();


	//disp->can[0]->Print("EDisp.ps");
	//int tmp; scanf("%d",&tmp);
	disp->graph[3]->Delete();
      }
      
    }

    //added in event structure

    disp->MakeGraph(9,dX,dY);
    disp->can[1]->cd(1);
    disp->graph[9]->SetMarkerColor(2);
    disp->graph[9]->Draw("p,same");

    io->GetEvtStrParam(ESt);
    ESt.push_back(ETheta);

    TLatex *es_evn,*es_det,*es_poshit,*es_hit,*es_Dene,*es_Dtime;
    TLatex *es_Dpos, *es_Dmom,*es_Time,*es_dTime;
    TLatex *es_Fpos, *es_Fmom, *es_Lpos, *es_Lmom;
    TLatex *es_theta;
    char esevn[30],esdet[30],esposhit[30],eshit[30],esDene[30],esDtime[30];
    char esDpos[40], esDmom[40],esTime[40],esdTime[40];
    char esFpos[40],esFmom[40],esLpos[40],esLmom[40];
    char estheta[40];

    double Dabs= ESt[8]*ESt[8]+ESt[9]*ESt[9]+ESt[10]*ESt[10];
    Dabs = TMath::Sqrt(Dabs);
    double Dmomx=ESt[8]/Dabs;
    double Dmomy=ESt[9]/Dabs;
    double Dmomz=ESt[10]/Dabs;

    double Fabs= ESt[16]*ESt[16]+ESt[17]*ESt[17]+ESt[18]*ESt[18];
    Fabs = TMath::Sqrt(Fabs);
    double Fmomx=ESt[16]/Fabs;
    double Fmomy=ESt[17]/Fabs;
    double Fmomz=ESt[18]/Fabs;

    double Labs= ESt[22]*ESt[22]+ESt[23]*ESt[23]+ESt[24]*ESt[24];
    Labs = TMath::Sqrt(Labs);
    double Lmomx=ESt[22]/Labs;
    double Lmomy=ESt[23]/Labs;
    double Lmomz=ESt[24]/Labs;


    double deltaT=ESt[12]-ESt[11];

    sprintf(esevn,"Event # : %d",(int)ESt[0]);
    sprintf(esDene,"Initial energy of e+: %3.1f MeV",ESt[3]);
    sprintf(esDtime,"Decay time: %3.1f ns",ESt[4]);
    sprintf(esDpos,"Decay position : (%3.1f, %3.1f, %3.1f)",
	    ESt[5],ESt[6],ESt[7]);
    sprintf(esDmom,"Initial vector of e+: (%3.3f, %3.3f, %3.3f)",
	    Dmomx,Dmomy,Dmomz);
    sprintf(esposhit,"Hit # in vanes (e+) : %d",(int)ESt[1]);
    sprintf(eshit,", (all) : %d",(int)ESt[2]);
    sprintf(esTime,"Time (first): %3.1f ns, (last): %3.1f ns"
	    ,ESt[11],ESt[12]);
    sprintf(esdTime,"   ( #DeltaT : %3.1f ns )",deltaT);
    sprintf(esFpos,"Position :First (%3.1f, %3.1f, %3.1f)",
	    ESt[13],ESt[14],ESt[15]);
    sprintf(esLpos,"Last (%3.1f, %3.1f, %3.1f)"
	    ,ESt[19],ESt[20],ESt[21]);
    sprintf(esFmom,"Vector :First (%3.3f, %3.3f, %3.3f)",
	    Fmomx,Fmomy,Fmomz);
    sprintf(esLmom," Last (%3.3f, %3.3f, %3.3f)",Lmomx,Lmomy,Lmomz);
    sprintf(estheta,"#theta of Hough fit line : %3.1f deg.",ESt[25]);

    es_evn = new TLatex(0.05,0.92,esevn);
    es_Dene = new TLatex(0.05,0.84,esDene);
    es_Dtime = new TLatex(0.05,0.76,esDtime);
    es_Dpos = new TLatex(0.05,0.68,esDpos);
    es_Dmom = new TLatex(0.05,0.6,esDmom);
    es_poshit = new TLatex(0.05,0.52,esposhit);
    es_hit = new TLatex(0.55,0.52,eshit);
    es_Time = new TLatex(0.05,0.44,esTime);
    es_dTime = new TLatex(0.05,0.38,esdTime);
    es_Fpos = new TLatex(0.05,0.3,esFpos);
    es_Lpos = new TLatex(0.25,0.25,esLpos);
    es_Fmom = new TLatex(0.05,0.17,esFmom);
    es_Lmom = new TLatex(0.2,0.12,esLmom);
    es_theta = new TLatex(0.05,0.04,estheta);


    es_evn->SetTextSize(0.05);
    es_Dene->SetTextSize(0.05);
    if(ESt[3]>150) es_Dene->SetTextColor(2);
    else es_Dene->SetTextColor(4);
    es_Dtime->SetTextSize(0.05);
    es_Dpos->SetTextSize(0.05);
    es_Dmom->SetTextSize(0.05);
    es_poshit->SetTextSize(0.05);
    es_hit->SetTextSize(0.05);
    es_Time->SetTextSize(0.05);
    es_dTime->SetTextSize(0.05);
    es_Fpos->SetTextSize(0.05);
    es_Lpos->SetTextSize(0.05);
    es_Fmom->SetTextSize(0.05);
    es_Lmom->SetTextSize(0.05);
    es_theta->SetTextSize(0.05);


    disp->can[1]->cd(4);
    es_evn->Draw();
    es_Dene->Draw();
    es_Dtime->Draw();
    es_Dpos->Draw();
    es_Dmom->Draw();
    es_poshit->Draw();
    es_hit->Draw();
    es_Time->Draw();
    es_dTime->Draw();
    es_Fpos->Draw();
    es_Lpos->Draw();
    es_Fmom->Draw();
    es_Lmom->Draw();
    es_theta->Draw();


    if(delflag4==1){
      sprintf(esdet,"Detected!!");
      es_det = new TLatex(0.5,0.92,esdet);
      es_det->SetTextSize(0.05);
      es_det->SetTextColor(616);
      es_det->Draw();
    }

    ofstream fstrct("EvSt.dat",ios::app);
    for(int ii=0;ii<26;ii++){
      fstrct << ESt[ii] << " ";
    }
    double oknot=0;
    if(delflag4==1) oknot=1;
    ESt.push_back(oknot);
    fstrct << ESt[26] << endl;

    fstrct.close();

    //disp->can[1]->Print("EDisp.ps");
    disp->can[1]->Update();
    disp->can[1]->WaitPrimitive();
    //int tmp; scanf("%d",&tmp);

    disp->graph[9]->Delete();
    es_evn->Delete();
    es_Dene->Delete();
    es_Dtime->Delete();
    es_Dpos->Delete();
    es_Dmom->Delete();
    es_poshit->Delete();
    es_hit->Delete();
    es_Time->Delete();
    es_dTime->Delete();
    es_Fpos->Delete();
    es_Lpos->Delete();
    es_Fmom->Delete();
    es_Lmom->Delete();
    es_theta->Delete();


    if(delflag4==1){
      es_det->Delete();
    }

  }

  cout << endl;
  cout << "all                       : " << counter0 << endl;
  cout << "Nhit cut                 : " << counter1 << " (" 
       << truecnt3 << ")" << endl;
  cout << "Hough & angle cut (true) : " << counter2 << " (" 
       << truecnt2 << ")" << endl;
  cout << "Clustering cut (true) : " << counter3 << " (" 
       << truecnt4 << ")" << endl;
  cout << "Final (true) : " << OKcount << " (" 
       << OKtrue << ")" << endl;
  cout << endl;
  cout << "True event (E>150MeV) :" << truecnt << endl;

  //disp->can[0]->Print("EDisp.ps)");
  //disp->can[1]->Print("EDisp.ps)");

  ofstream fout("Counter.dat");
  fout << "all                       : " << counter0 << endl;
  fout << "Nhit cut                 : " << counter1 << " (" 
       << truecnt3 << ")" << endl;
  fout << "Hough & angle cut (true) : " << counter2 << " (" 
       << truecnt2 << ")" << endl;
  fout << "Clustering cut (true) : " << counter3 << " (" 
       << truecnt4 << ")" << endl;
  fout << "Final (true) : " << OKcount << " (" 
       << OKtrue << ")" << endl;
  fout << endl;
  fout << "True event (E>150MeV) :" << truecnt << endl;
  fout.close();

  /*
  ofstream fout("tmp.dat");
  for(int i=0; i<tempnum; i++){
    fout << temp[i] << endl;
  }
  fout.close();
  ofstream fout2("tmp2.dat");
  for(int i=0; i<tempnum2; i++){
    fout2 << temp2[i] << endl;
  }
  fout2.close();
  */

  //  disp->MakeGraph(0,Phi,Z);



  //--End
 
  //  io->CloseFile();
  disp->ROOT_End();  
  delete disp; delete io; delete algo;
  
 
}

void Initial(){

  cout << "////////////////////////////////////////" << endl;
  cout << "//                                    //" << endl;
  cout << "//       Track Finding Tool           //" << endl;
  cout << "//               For Muon g-2         //" << endl;
  cout << "//                        ver.0       //" << endl;
  cout << "//                                    //" << endl;
  cout << "//                 writen by K. Ueno  //" << endl;
  cout << "//                                    //" << endl;
  cout << "////////////////////////////////////////" << endl;

}

void ClearVector(){
  X.clear(); Y.clear(); Z.clear(); R.clear(); Phi.clear();
  VaneNum.clear(); AfterR.clear(); AfterZ.clear(); AfterX.clear(); AfterY.clear();
  Id.clear(); Thr.clear(); UnitNum.clear();
  dX.clear(); dY.clear(); dZ.clear();
  redVaneNum.clear();
  //added in ver.0.12
  Zcp.clear(); Phicp.clear();
  //added in event structure
  ESt.clear();

}
