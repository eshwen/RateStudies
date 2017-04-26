
#include <iostream>
#include <fstream>
#include <utility>
#include <tuple>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"


using namespace std;
void scan(int targetRate){

  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias14Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events/";
  

  //ZeroBias sample L1
  TFile *file = TFile::Open(Form("%srateL1_tausscaleToLumi18E33.root", directory.Data()),"read");
  if (file == NULL) cout<<"File not found"<<endl;

  int xcut = 0;
  int ycut = 0;
  ofstream fOut(Form("%srate%dkHz_etTaus_etTauTau_1d8E34.txt",directory.Data(),targetRate));

  TH2D * hInput = (TH2D*) file->Get("Rate_DiTau2D");
  if (hInput == NULL) cout<<"Histogram not found"<<endl;
  for(int j = 1 ; j <= hInput->GetNbinsY() ; ++j){

    for(int  i = 1 ; i <= hInput->GetNbinsX() ; ++i)
      {
	xcut= 0;
	ycut= 0;

	if(hInput->GetBinContent(i,j)<=(float)targetRate){
	  xcut= hInput->GetXaxis()->GetBinLowEdge(i);
	  ycut= hInput->GetYaxis()->GetBinLowEdge(j);
	  cout<<"pT_tau "<<ycut<<" GeV; pT_tautau "<<xcut<<" GeV; rate "<<hInput->GetBinContent(i,j)<<endl;	 
	  fOut<<ycut<<"\t"<<xcut<<"\t"<<hInput->GetBinContent(i,j)<<endl;	 
	  break;
	}
	
      }

  }
  fOut.close();
}



void scanTot(int targetRate){

  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias14Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events/";
  

  //  ZeroBias sample L1
  TFile *file = TFile::Open(Form("%srateL1_tausscaleToLumi20E33.root", directory.Data()),"read");
  if (file == NULL) cout<<"File not found"<<endl;

  int xcut = 0;
  int ycut = 0;
  ofstream fOut(Form("%srateTot%dkHz_etTaus_etTauTau_2E34.txt",directory.Data(),targetRate));

  TH2D * hInputRate = (TH2D*) file->Get("Rate_DiTau2D");
  if (hInputRate == NULL) cout<<"Histogram not found"<<endl;
  TH2D * hInputPure[8];
  cout<<"high boost"<<endl;
  for(int xx = 0; xx<8; xx++){
    hInputPure[xx] = (TH2D*) file->Get(Form("PureRate_DiTau2D_wrt%d",30+xx));
    if (hInputPure[xx] == NULL) cout<<"Histogram not found"<<endl;
  }
    fOut<<"pT_tau [GeV]\tpT_tautau [GeV]\tDitauXXer cut [GeV]\tRate DitauXXer [kHz]\tpure rate wrt DiTauXXer [kHz]\tTotal rate (of the OR)[kHz]"<<endl;
  for(int xx = 0; xx<8; xx++){      
    xcut= 70;
    ycut= 25;
    // cout <<i<< endl;
    double totRate =hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30)) + hInputPure[xx]->GetBinContent(hInputRate->GetXaxis()->FindBin(xcut),hInputRate->GetYaxis()->FindBin(ycut));

    //if((totRate<=(float)targetRate)){
    cout<<"pT_tau "<<ycut<<" GeV; pT_tautau "<<xcut<<" GeV; rate ditau "<<hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30))<<" kHz; pure rate wrt "<<xx+30<<" "<<hInputPure[xx]->GetBinContent(hInputRate->GetXaxis()->FindBin(xcut),hInputRate->GetYaxis()->FindBin(ycut))<<" kHz; total rate"<<totRate<<endl;

    fOut<<ycut<<"\t"<<xcut<<"\t"<<xx+30<<"\t"<<hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30))<<"\t"<<hInputPure[xx]->GetBinContent(hInputRate->GetXaxis()->FindBin(xcut),hInputRate->GetYaxis()->FindBin(ycut))<<"\t"<<totRate<<endl;
	   
    //}
  }
  cout<<"* * * * *"<<endl;
  cout<<"moderate boost"<<endl;
  int flag = 0;
  for(int j = 1 ; j <= hInputRate->GetYaxis()->FindBin(35.) ; ++j){
    for(int xx = 0; xx<8; xx++){      
      for(int  i = 1 ; i <= hInputRate->GetNbinsX() ; ++i)
	{
	  flag = 0;
	  xcut= 0;
	  ycut= 0;
	  // cout <<i<< endl;
	  double totRate =hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30)) + hInputPure[xx]->GetBinContent(i,j);
	    xcut= hInputRate->GetXaxis()->GetBinLowEdge(i);
	    ycut= hInputRate->GetYaxis()->GetBinLowEdge(j);	  
	    if((totRate<=(float)targetRate)&&(xx+30>=31) && (xx+30<=35) && (xcut<=90) && (xcut>=70) && (ycut<=28) && (ycut>=25)){
	    if(ycut<xx+30){
	      cout<<"pT_tau "<<ycut<<" GeV; pT_tautau "<<xcut<<" GeV; rate ditau "<<hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30))<<" kHz; pure rate wrt "<<xx+30<<" "<<hInputPure[xx]->GetBinContent(hInputRate->GetXaxis()->FindBin(xcut),hInputRate->GetYaxis()->FindBin(ycut))<<" kHz"<<endl;
	      fOut<<ycut<<"\t"<<xcut<<"\t"<<xx+30<<"\t"<<hInputRate->GetBinContent(hInputRate->GetXaxis()->FindBin(0.),hInputRate->GetYaxis()->FindBin(xx+30))<<"\t"<<hInputPure[xx]->GetBinContent(hInputRate->GetXaxis()->FindBin(xcut),hInputRate->GetYaxis()->FindBin(ycut))<<"\t"<<totRate<<endl;
	      flag = 1;
	      break;
	    }
	    
	    }
	    if (flag == 1) break;
	}
if (flag == 1) break;
    }
    cout<<"* * *"<<endl;
  }
  fOut.close();
  cout<<"Saved in "<<directory.Data()<<"rateTot"<<targetRate<<"kHz_etTaus_etTauTau_2E34.txt"<<endl;
}

