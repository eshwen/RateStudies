
// evaluates rate  - compile with c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
// and launch the executable 

#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

using namespace std;

int main(int argc, char** argv){
  //double scaleRunI = 1./1.713;
  double scaleRunI = 1.075/1.713;
  bool doScaleRunI = false;
  
  cout << "Do scale: " << doScaleRunI << endl;
  cout << "Scale factor: " << scaleRunI << endl;
  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias25Gen2016HighPU_ibx0_100000Events/";
  //ZeroBias sample L1
  TFile *file = TFile::Open(Form("%sL1ntuples_0.root", directory.Data()),"read");
  TTree * tInput = (TTree*) file->Get("L1Tree/L1Tree");
  
  TString fOutName = "rateL1_";
  fOutName = directory+fOutName;
  if (doScaleRunI) fOutName += "Scaled.root";
  else fOutName += "NotScaled.root";
    
  
  Int_t stage2_tauN ;
  std::vector<double>* stage2_tauEt=0;
  std::vector<double>* stage2_tauEta=0;
  std::vector<double>* stage2_tauPhi=0;
  
  
  Int_t stage2_jetN ;
  std::vector<Float_t> *stage2_jetEt=0;
  std::vector<Float_t> *stage2_jetEta=0;
  std::vector<Float_t> *stage2_jetPhi=0;
  
  // set branch and variables
  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  
  
  tInput ->SetBranchAddress("stage2_tauN", &stage2_tauN , &b_stage2_tauN );
  tInput ->SetBranchAddress("stage2_tauEta", &stage2_tauEta, &b_stage2_tauEta);
  tInput ->SetBranchAddress("stage2_tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
  tInput ->SetBranchAddress("stage2_tauEt", &stage2_tauEt, &b_stage2_tauEt);
  tInput ->SetBranchAddress("stage2_jetN", &stage2_jetN , &b_stage2_jetN);
  tInput ->SetBranchAddress("stage2_jetEta", &stage2_jetEta, &b_stage2_jetEta);
  tInput ->SetBranchAddress("stage2_jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
  tInput ->SetBranchAddress("stage2_jetEt", &stage2_jetEt, &b_stage2_jetEt);
 
  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOut = new TFile (Form("%s",fOutName.Data()), "recreate");
  
  // the tree histograms are filed in identical way, structure is simply copied from the Stage 2 rate Evaluator
  TH1D* LeadTauPt_Pass = new TH1D ("LeadTauPt_Pass", "LeadTauPt_Pass", 50, 0, 100);
  TH1D* SingleTauPt_Pass = new TH1D ("SingleTauPt_Pass", "SingleTauPt_Pass", 50, 0, 100);
  TH1D* SubleadTauPt_Pass = new TH1D ("SubleadTauPt_Pass", "SubleadTauPt_Pass", 50, 0, 100);
  TH1D* DiTauPt_Pass = new TH1D ("DiTauPt_Pass", "DiTauPt_Pass", 50, 0, 100);
  TH2D* DiTau2D_Pass = new TH2D ("DiTau2D_Pass", "DiTau2D_Pass", 50, 0, 100, 50, 0 , 100);
  TH1D* DiTauSingleJetPt_Pass = new TH1D ("DiTauSingleJetPt_Pass", "DiTauSingleJetPt_Pass", 50, 0, 100);
  TH2D* DiTauSingleJet2D_Pass = new TH2D ("DiTauSingleJet2D", "single tau plus jet", 50, 0, 100, 10, 0, 5);
  TH2D* SingleTauSingleJet2D_Pass = new TH2D ("SingleTauSingleJet2D", "double tau const", 50, 0, 100, 10, 0, 5);

  TH1D* RelativeRate_singleTau = new TH1D ("RelativeRate_singleTau", "Relative rate - single tau", 50, 0, 100);
  TH1D* RelativeRate_diTau = new TH1D ("RelativeRate_diTau", "Relative rate - double tau", 50, 0, 100);
  TH1D* RelativeRate_diTauSublead = new TH1D ("RelativeRate_diTauSublead", "Relative rate - double tau", 50, 0, 100);
  TH2D* RelativeRate_diTau2D = new TH2D ("RelativeRate_diTau2D", "Relative rate - double tau", 50, 0, 100, 50, 0 ,100);
  TH1D* RelativeRate_diTauSingleJet = new TH1D ("RelativeRate_diTauSingleJet", "Relative rate - double tau plus jet", 50, 0, 100);
  TH1D* RelativeRate_diTauRef = new TH1D ("RelativeRate_diTauRef", "Relative rate - double tau const", 50, 0, 100);
  TH2D* RelativeRate_diTauSingleJet2D = new TH2D ("RelativeRate_diTauSingleJet2D", "Relative rate - single tau plus jet", 50, 0, 100, 50, 0, 100);
  TH2D* RelativeRate_SingleTauSingleJet2D = new TH2D ("RelativeRate_SingleTauSingleJet2D", "Relative rate - double tau const", 50, 0, 100, 50, 0, 100);

  

  
  
  // analyze data    
  long int nEvents = tInput->GetEntries(); 
  std::vector<double> ptTau_pass; 
  std::vector<double> ptJet_pass;
  double singletau = 0;
  double ditau = 0;  
  // loop on all events
  for (long int iEv = 0; iEv < nEvents; iEv++){
    
    tInput->GetEntry(iEv);
    
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
    
    // clear pt vector for this event 
    ptTau_pass.clear();
    ptJet_pass.clear();
    
    
    for (long int iL1 = 0; iL1 < stage2_tauN; iL1++){
      // selections
      double abseta = TMath::Abs(stage2_tauEta->at(iL1) );
      double tauPt  = stage2_tauEt->at(iL1);
      
      if (doScaleRunI) tauPt = scaleRunI*tauPt;
      
      if (abseta < 2.1){	      
	ptTau_pass.push_back (tauPt);
      }
    }
    
    for (long int iL1 = 0; iL1 < stage2_jetN; iL1++){
      // selections
      double abseta = TMath::Abs( (*stage2_jetEta)[iL1] );
      double jetPt  = (*stage2_jetPhi)[iL1];
      
      if (doScaleRunI) jetPt = scaleRunI*jetPt;
      
      if (abseta < 2.4){
	ptJet_pass.push_back (jetPt);
      }
    }
    
    
    // now that all taus in the event are analyzed, fill the histograms of lead and sublead pt (if missing, fill with a -1)
    std::sort (ptTau_pass.begin(), ptTau_pass.end());
    std::sort (ptJet_pass.begin(), ptJet_pass.end());
    
    
    if (ptTau_pass.size() >= 2 )
      {
	ditau += 1;
	singletau += 1;
	    LeadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
	    SingleTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
	    SubleadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin() + 1) );
	    DiTauPt_Pass -> Fill (*(ptTau_pass.rbegin()) +*(ptTau_pass.rbegin() + 1) );
	    DiTau2D_Pass -> Fill (*(ptTau_pass.rbegin()),*(ptTau_pass.rbegin() + 1) );
	    if (ptJet_pass.size()>=1){
	      DiTauSingleJetPt_Pass->Fill( *(ptJet_pass.rbegin()) );
	      DiTauSingleJet2D_Pass->Fill( *(ptJet_pass.rbegin()), *(ptTau_pass.rbegin()) +*(ptTau_pass.rbegin() + 1) );
	    }else{
	      DiTauSingleJetPt_Pass->Fill( -1 );
	      DiTauSingleJet2D_Pass->Fill( -1,-1 );
	    }
      } 
    
    else if (ptTau_pass.size() == 1){
      singletau +=1;
      LeadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
      SingleTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
      SubleadTauPt_Pass -> Fill ( -1 );
      DiTauPt_Pass -> Fill ( -1); 
      DiTauSingleJetPt_Pass->Fill( -1 );
      DiTau2D_Pass->Fill(-1,-1);
      if (ptJet_pass.size()>=1){
	SingleTauSingleJet2D_Pass->Fill( *(ptJet_pass.rbegin()), *(ptTau_pass.rbegin()) );
      }else{
	SingleTauSingleJet2D_Pass->Fill( -1,-1 );
      }       
    } 
    
    else
      {
	LeadTauPt_Pass -> Fill ( -1 ); 
	SingleTauPt_Pass -> Fill ( -1 ); 
	SubleadTauPt_Pass -> Fill ( -1 );        
	DiTauPt_Pass -> Fill ( -1 );        
	DiTauSingleJetPt_Pass->Fill( -1 );
	DiTauSingleJet2D_Pass->Fill( -1,-1 );
	SingleTauSingleJet2D_Pass->Fill( -1,-1 );
	DiTau2D_Pass->Fill(-1,-1);
      }       
  }
  
  // compute rate plots
  cout << "Computing rates..." << endl; 
  
  for (int i = 1; i <= 50; i++){
    // int Tot_singleTau = LeadTauPt_Pass -> Integral (0, 51);
    int Tot_diTauSublead = SubleadTauPt_Pass -> Integral (1,51);
    
    int Tot_singleTau = SingleTauPt_Pass -> Integral (1,51);
    int Tot_diTau = DiTauPt_Pass -> Integral (1,51);        
    int Tot_diTauSingleJet = DiTauSingleJetPt_Pass -> Integral (1,51);        
    int Tot_diTauSingleJet2D = DiTauSingleJet2D_Pass -> Integral (1,51);       
    int Tot_SingleTauSingleJet2D = SingleTauSingleJet2D_Pass -> Integral (1,51);        
    int Tot_diTau2D = DiTau2D_Pass -> Integral (1,51);        
    
    //  double relRateSingleTau = 1.*(LeadTauPt_Pass->Integral(i, 51))/Tot_singleTau;
    //  double relRateDiTau = 1.*(SubleadTauPt_Pass->Integral(i, 51))/Tot_diTau;
    
    double relRateSingleTau = 1.*(SingleTauPt_Pass->Integral(i, 51))/Tot_singleTau;
    double relRateDiTau = 1.*(DiTauPt_Pass->Integral(i, 51))/Tot_diTau;
    double relRateDiTauSingleJet = 1.*(DiTauSingleJetPt_Pass->Integral(i, 51))/Tot_diTauSingleJet;
    double relRateDiTauSublead =  1.*(SubleadTauPt_Pass->Integral(i, 51))/Tot_diTauSublead;
    
    RelativeRate_singleTau -> SetBinContent (i, relRateSingleTau);
    RelativeRate_diTau -> SetBinContent (i, relRateDiTau);        
    RelativeRate_diTauSublead -> SetBinContent (i, relRateDiTauSublead);        
    
    RelativeRate_diTauSingleJet -> SetBinContent (i, relRateDiTauSingleJet);        
    for (int j = 1; j<=50;j++ ){
      double relRateDiTauSingleJet2D = 1.*(DiTauSingleJet2D_Pass->Integral(i, 51,j,51))/Tot_diTauSingleJet2D;
      double relRateSingleTauSingleJet2D = 1.*(SingleTauSingleJet2D_Pass->Integral(i, 51,j,51))/Tot_SingleTauSingleJet2D;
      double relRateDiTau2D = 1.*(DiTau2D_Pass->Integral(i, 51,j,51))/Tot_diTau2D;
      RelativeRate_SingleTauSingleJet2D -> SetBinContent (i,j,relRateSingleTauSingleJet2D  );
      RelativeRate_diTauSingleJet2D -> SetBinContent (i, j,relRateDiTauSingleJet2D);        
      RelativeRate_diTau2D -> SetBinContent (i, j, relRateDiTau2D);        
    }
  }
  cout<<"Fraction of events with at least 1 tau: "<<singletau/nEvents<<endl;
  cout<<"Fraction of events with at least 2 taus: "<<ditau/nEvents<<endl;   
  fOut -> Write();
}
