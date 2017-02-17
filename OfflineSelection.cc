
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


using namespace std;

bool SortByPt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }

int main(int argc, char** argv){

   
   TString LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_VBFHToTauTau_-1Events_0Skipped_1487239220.05/";
   TFile *fOut = TFile::Open(Form("%sVBFHTauTau_tau25tau25_jet30.root", LLRdir.Data()),"recreate");
   TTree *tOutput = new TTree("HTauTauTreeSel","HTauTauTreeSel");
   TFile *file = TFile::Open(Form("%sHTauTauAnalysis_total.root", LLRdir.Data()),"read");
   if (file == NULL) cout<<"File not found"<<endl;


   TTree * tInput = (TTree*) file->Get("HTauTauTree/HTauTauTree");
   

   ULong64_t       EventNumber;
   Int_t           RunNumber;
   Int_t           lumi;
   std::vector<int>     *particleType=0;
   std::vector<float>   *daughters_px=0;
   std::vector<float>   *daughters_py=0;
   std::vector<float>   *daughters_pz=0;
   std::vector<float>   *daughters_e=0;
   std::vector<Long64_t> *tauID=0;
   Int_t           JetsNumber;
   std::vector<float>   *jets_px=0;
   std::vector<float>   *jets_py=0;
   std::vector<float>   *jets_pz=0;
   std::vector<float>   *jets_e=0;
   Int_t stage2_tauN;
   std::vector<Float_t>* stage2_tauEt=0;
   std::vector<Float_t>* stage2_tauEta=0;
   std::vector<Float_t>* stage2_tauPhi=0;
   std::vector<int> *stage2_tauIso=0;  
   
   Int_t stage2_jetN;
   std::vector<Float_t> *stage2_jetEt=0;
   std::vector<Float_t> *stage2_jetEta=0;
   std::vector<Float_t> *stage2_jetPhi=0;
  

  
  


   TBranch        *b_EventNumber;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_particleType;   //!
   TBranch        *b_daughters_px;   //!
   TBranch        *b_daughters_py;   //!
   TBranch        *b_daughters_pz;   //!
   TBranch        *b_daughters_e;   //!
   TBranch        *b_tauID;   //!
   TBranch        *b_JetsNumber;   //!
   TBranch        *b_jets_px;   //!
   TBranch        *b_jets_py;   //!
   TBranch        *b_jets_pz;   //!
   TBranch        *b_jets_e;   //!
  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;
  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;

   
   tInput->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   tInput->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   tInput->SetBranchAddress("lumi", &lumi, &b_lumi);
   tInput->SetBranchAddress("particleType", &particleType, &b_particleType);
   tInput->SetBranchAddress("daughters_px", &daughters_px, &b_daughters_px);
   tInput->SetBranchAddress("daughters_py", &daughters_py, &b_daughters_py);
   tInput->SetBranchAddress("daughters_pz", &daughters_pz, &b_daughters_pz);
   tInput->SetBranchAddress("daughters_e", &daughters_e, &b_daughters_e);
   tInput->SetBranchAddress("tauID", &tauID, &b_tauID);
   tInput->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
   tInput->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
   tInput->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
   tInput->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
   tInput->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
   tInput ->SetBranchAddress("stage2_tauN", &stage2_tauN , &b_stage2_tauN );
   tInput ->SetBranchAddress("stage2_tauEta", &stage2_tauEta, &b_stage2_tauEta);
   tInput ->SetBranchAddress("stage2_tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
   tInput ->SetBranchAddress("stage2_tauEt", &stage2_tauEt, &b_stage2_tauEt);
   tInput ->SetBranchAddress("stage2_tauIso", &stage2_tauIso, &b_stage2_tauIso);
   tInput ->SetBranchAddress("stage2_jetN", &stage2_jetN , &b_stage2_jetN);
   tInput ->SetBranchAddress("stage2_jetEta", &stage2_jetEta, &b_stage2_jetEta);
   tInput ->SetBranchAddress("stage2_jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
   tInput ->SetBranchAddress("stage2_jetEt", &stage2_jetEt, &b_stage2_jetEt);   

   tOutput->Branch("EventNumber", &EventNumber);
   tOutput->Branch("RunNumber", &RunNumber);
   tOutput->Branch("lumi", &lumi);
   tOutput->Branch("particleType", &particleType);
   tOutput->Branch("daughters_px", &daughters_px);
   tOutput->Branch("daughters_py", &daughters_py);
   tOutput->Branch("daughters_pz", &daughters_pz);
   tOutput->Branch("daughters_e", &daughters_e);
   tOutput->Branch("tauID", &tauID);
   tOutput->Branch("JetsNumber", &JetsNumber);
   tOutput->Branch("jets_px", &jets_px);
   tOutput->Branch("jets_py", &jets_py);
   tOutput->Branch("jets_pz", &jets_pz);
   tOutput->Branch("jets_e", &jets_e);
   tOutput ->Branch("stage2_tauN", &stage2_tauN);
   tOutput ->Branch("stage2_tauEta", &stage2_tauEta);
   tOutput ->Branch("stage2_tauPhi", &stage2_tauPhi);
   tOutput ->Branch("stage2_tauEt", &stage2_tauEt);
   tOutput ->Branch("stage2_tauIso", &stage2_tauIso);
   tOutput ->Branch("stage2_jetN", &stage2_jetN);
   tOutput ->Branch("stage2_jetEta", &stage2_jetEta);
   tOutput ->Branch("stage2_jetPhi", &stage2_jetPhi);
   tOutput ->Branch("stage2_jetEt", &stage2_jetEt);   


   std::vector <TLorentzVector> tau;

   
   long int nEvents = tInput->GetEntries(); 
   
     for (long int iEv = 0; iEv < nEvents; iEv++){
     tau.clear();
     file->cd();
     tInput->GetEntry(iEv);
     if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
     
     for(int i = 0; i<daughters_px->size(); i++ ){
       if (particleType->at(i) ==2){
	 
	 TLorentzVector tlv_tau;
	 tlv_tau.SetPxPyPzE(
			    daughters_px->at(i),
			    daughters_py->at(i),
			    daughters_pz->at(i),
			    daughters_e->at(i)
			    );
	 if(tlv_tau.Pt()>25 && fabs(tlv_tau.Eta())<2.1 && ((tauID->at(i)>>5)&1) && ((tauID->at(i)>>7)&1) && ((tauID->at(i)>>5)&11)){
	 tau.push_back(tlv_tau);
	 }
       }
     }
     std::sort(tau.begin(),tau.end(),SortByPt);
     if(tau.size()>=2){     
       TLorentzVector tlv_tauPair = tau[0]+tau[1];
       for(int j=0; j<JetsNumber;j++){
	 TLorentzVector jet;
	 jet.SetPxPyPzE(
			jets_px->at(j),
			jets_py->at(j),
			jets_pz->at(j),
			jets_e->at(j)
			) ;
	 //	   if(tlv_tauPair.Pt()>100 && jet.Pt()>30){
	   if(jet.Pt()>30){
	     
	     tOutput->Fill();
	     break;
	   }
       }       
     }

     }
     fOut->cd();
     fOut->Write();
     fOut->Close();
}



