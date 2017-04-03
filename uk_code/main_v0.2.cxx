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

int vane=48;
//int vane=32;

vector<double> X, Y, Z, R, Phi, VaneNum, AfterR, AfterZ, AfterX, AfterY;
vector<double> dX, dY, dZ, redVaneNum, Time;
vector<int> Id, Thr, UnitNum;


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


  //  disp->can[0]->Print("EDisp.ps(");

  io->OpenFile(0);
  io->OpenTree();
  int delflag=0,delflag2=0,delflag3=0,delflag4=0;

  int OKcount=0,OKtrue=0;
  int counter0=0,counter1=0, counter2=0,counter3=0;
  int truecnt=0, truecnt2=0, truecnt3=0,truecnt4=0;
  int temp[500],temp2[500];
  int tempnum=0,tempnum2=0;

  int evnum = io->GetEventNumber();

  double dpos[3];

  TLatex *texNum,*texEne,*texDet;


  for(int i=0;i<evnum;i++){ //loop
  //for(int i=0;i<120;i++){ //loop
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

    io->ReadOneEvent(i);
    io->GetParam(X,Y,Z,R,Phi,Id,Thr,Time);
    io->GetDecayPoint(dpos[0],dpos[1],dpos[2]);
    dX.push_back(dpos[0]);
    dY.push_back(dpos[1]);
    dZ.push_back(dpos[2]);
    double dE = io->GetPositronInitialEnergy();
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
      // cout << disp->his[0]->GetRMS() << endl;
      
      vector<double> redZ,redPhi;
      for(int i=0; i<(int)Resi.size();i++){
	if(fabs(Resi[i]-disp->his[0]->GetMean())<=3*(disp->his[0]->GetRMS())){
	  redZ.push_back(Z[i]);
	  redPhi.push_back(Phi[i]);
	}
      }
      double angleflag=0;
      angleflag = TMath::ATan(par1*180/TMath::Pi())*180/TMath::Pi();
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
	if(tmpmax>3){
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
	  int overii[10*VaneNum.size()];
	  int overiicnt=0;
	  for(int ii=0;ii<(int)VaneNum.size()*10;ii++){
	    overii[ii]=VaneNum.size();
	  }

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
	    if(minsub>10)break;
	    clsVaneNum.push_back(hashi+1);
	    clsZ.push_back(minsubY);
	  }

	  for(int ii=0;ii<(int)VaneNum.size()*10;ii++){
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
	    if(minsub>10)break;
	    clsVaneNum.insert(clsVaneNum.begin(),hashi-1);
	    clsZ.insert(clsZ.begin(),minsubY);
	  }

	  disp->MakeGraph(7,clsVaneNum,clsZ);
	  disp->can[0]->cd(10);
	  disp->graph[7]->Draw("p,same");

	  double mintime[3],mintimei[3],minZ[3],minVane[3];
	  for(int ii=0;ii<3;ii++){
	    mintime[ii]=1000000;
	    mintimei[ii]=0;
	    minZ[ii]=0;
	    minVane[ii]=0;
	  }

	  //GetTrue3hits
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

	  //Compare
	  int trueflag=0;
	  for(int ii=0;ii<(int)clsZ.size();ii++){
	    for(int jj=0;jj<3;jj++){
	      if(clsVaneNum[ii]==minVane[jj] && clsZ[ii]==minZ[jj]){
		trueflag++;
	      }
	    }
	  }

	  if(trueflag==3){
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
      int tmp; scanf("%d",&tmp);
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
	int tmp; scanf("%d",&tmp);
	disp->graph[0]->Delete();
	disp->graph[1]->Delete();
	disp->graph[3]->Delete();

      }
      else if(dE>150 && X.size()<=0){
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
	int tmp; scanf("%d",&tmp);
	disp->graph[3]->Delete();
      }
      
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
}
