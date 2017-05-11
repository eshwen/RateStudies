#include <sstream>
#include <fstream>
#include <map>
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
#include "interface/object.h"
#include "interface/utils.h"

using namespace std;


int main(int argc, char** argv){


  


  //ZeroBias sample L1
  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias26Apr2017HighPU_ibx0_BunchTrain0-5_2016H9Nov_-1Events/";

  //Takashi's emulated
  TString fileList = "fileLists/L1NtuplesL1Menu2017ZeroBiasBunchTrainsX.txt";
   
    


  TChain * emuInput;
  TChain * unInput;
  TChain * lInput; 
  

  emuInput = new TChain ("l1UpgradeEmuTree/L1UpgradeTree");
  appendFromFileList(emuInput, fileList);
  unInput = new TChain ("l1UpgradeTree/L1UpgradeTree");
  appendFromFileList(unInput, fileList);
  lInput = new TChain ("l1EventTree/L1EventTree");
  appendFromFileList(lInput, fileList);
  
  TString fOutNameVBF;
  
  
  
  fOutNameVBF = directory+"L1NTuplesTakPlots.root";
  
  
  
  lInput->SetMakeClass(1);
  emuInput->SetMakeClass(1);
  unInput->SetMakeClass(1);  
  
  Int_t lumi ;
  UShort_t emu_jetN ;
  std::vector<float> emu_jetEt;
  std::vector<float> emu_jetEta;
  std::vector<float> emu_jetPhi;

  UShort_t un_jetN ;
  std::vector<float> un_jetEt;
  std::vector<float> un_jetEta;
  std::vector<float> un_jetPhi;

  
  
  
  // set branch and variables
  TBranch *b_lumi ;

  TBranch *b_emu_jetN ;
  TBranch *b_emu_jetEt;
  TBranch *b_emu_jetEta;
  TBranch *b_emu_jetPhi;

  TBranch *b_un_jetN ;
  TBranch *b_un_jetEt;
  TBranch *b_un_jetEta;
  TBranch *b_un_jetPhi;
  
  

  

  lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
  
  emuInput ->SetBranchAddress("nJets", &emu_jetN , &b_emu_jetN);
  emuInput ->SetBranchAddress("jetEta", &emu_jetEta, &b_emu_jetEta);
  emuInput ->SetBranchAddress("jetPhi", &emu_jetPhi, &b_emu_jetPhi);
  emuInput ->SetBranchAddress("jetEt", &emu_jetEt, &b_emu_jetEt);
  
  unInput ->SetBranchAddress("nJets", &un_jetN , &b_un_jetN);
  unInput ->SetBranchAddress("jetEta", &un_jetEta, &b_un_jetEta);
  unInput ->SetBranchAddress("jetPhi", &un_jetPhi, &b_un_jetPhi);
  unInput ->SetBranchAddress("jetEt", &un_jetEt, &b_un_jetEt);
  
  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");
  
 
  TH1D* EnScaleB = new TH1D ("EnScaleB", "enscaleB", 50,-1,1);
  TH1D* EnScaleE = new TH1D ("EnScaleE", "enscaleE", 50,-1,1); 
  
  
  // analyze data    
  
  std::vector<object> jet20un;
  std::vector<object> jet20emu;    
  std::vector< std::tuple<double,int,int> > mjj_pass;
  
   bool L1_DoubleJet_80_20_20j20_mjj400 = false;
   
   
   
   for (Long64_t iEv = 0 ;true ; ++iEv){
     lumi =0;			   
     
     emu_jetN =0;		   
     emu_jetEt.clear(); 
     emu_jetEta.clear();
     emu_jetPhi.clear();
     
     un_jetN =0;		   
     un_jetEt.clear(); 
     un_jetEta.clear();
     un_jetPhi.clear();
     int got = 0;
     
     got =  lInput->GetEntry(iEv);
     emuInput->GetEntry(iEv);
     unInput->GetEntry(iEv);
     
     if (got == 0) break;
     
     if (iEv%1000 == 0) cout << iEv << endl;
     
          
     jet20un.clear();
     jet20emu.clear();
     mjj_pass.clear();
     
     if(lumi<56 || lumi>69) continue;
   
     L1_DoubleJet_80_20_20j20_mjj400 = false;
     
     
     
     

     for (long int iL1 = 0; iL1 < un_jetEt.size(); iL1++){ //loop on jets emulated
       // selections
       double jetPtUn  = un_jetEt.at(iL1);
       
       if(jetPtUn>20.) jet20un.push_back(object(un_jetEt.at(iL1),un_jetEta.at(iL1),un_jetPhi.at(iL1),-999)) ;
     }
     
     for (long int iL1emu = 0; iL1emu < emu_jetEt.size(); iL1emu++){ //loop on jets unpacked
	double jetPtEmu  = emu_jetEt.at(iL1emu);
	if(jetPtEmu>20.) jet20emu.push_back(object(emu_jetEt.at(iL1emu),emu_jetEta.at(iL1emu),emu_jetPhi.at(iL1emu),-999)) ;
      }


 
       
    if (jet20un.size() >= 2){
      for (int iJet = 0; iJet <jet20un.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet20un.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet20un[iJet].Et(),
			      jet20un[iJet].Eta(),
			      jet20un[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet20un[kJet].Et(),
			      jet20un[kJet].Eta(),
			      jet20un[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_pass.push_back(make_tuple(jetPair.M(),iJet,kJet));

	  }
	  
	}
	
      }

      std::sort(mjj_pass.begin(),mjj_pass.end());

      if(jet20un[0].Et()>80 && std::get<0>(*(mjj_pass.rbegin()))>400)     L1_DoubleJet_80_20_20j20_mjj400 = true;
      
    }
    float res = 0;
    for (long int iL1 = 0; iL1 < un_jetEt.size(); iL1++){ //loop on jets unpacked
      object::SortByDeltaR DeltaRJetUn(jet20un[iL1].Phi(),jet20un[iL1].Eta());
      std::sort(jet20emu.begin(),jet20emu.end(),DeltaRJetUn);
      if(jet20emu.size()>0){
	if(jet20emu[0].DeltaR(jet20un[iL1])<0.3){
	res = (jet20un[iL1].Et()-jet20emu[0].Et())/jet20un[iL1].Et();
	if(fabs(jet20un[iL1].Eta())<3){
	    EnScaleB->Fill(res);
	}else{
	  EnScaleE->Fill(res);
	}
	
      }
      }
    }
	  
	  
   }
   
   fOutVBF->cd();
   fOutVBF -> Write();
}

