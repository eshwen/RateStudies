
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
  TFile *file = TFile::Open(Form("%srateL1_taus.root", directory.Data()),"read");
  if (file == NULL) cout<<"File not found"<<endl;

  int xcut = 0;
  int ycut = 0;
  ofstream fOut(Form("%srate%dkHz_etTaus_ptTauTau.txt",directory.Data(),targetRate));

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
	  fOut<<ycut<<"\t"<<xcut<<endl;	 
	  break;
	}
	
      }

  }
  fOut.close();
}
