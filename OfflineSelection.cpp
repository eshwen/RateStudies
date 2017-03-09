
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
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "interface/object.h"


using namespace std;

bool SortByPt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }

int main(int argc, char** argv){

  int pt_tautau_cut = 0;
  int pt_tau_cut = 20;
  int pt_jet_cut = 30;
  TString LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_VBFHToTauTau_-1Events_0Skipped_1487239220.05/";
  TFile *file = TFile::Open(Form("%sHTauTauAnalysis_total.root", LLRdir.Data()),"read");
  if (file == NULL) cout<<"File not found"<<endl;
  TFile *fPlotsVBF = TFile::Open(Form("%sVBFHTauTauSignal_VBFseed.root", LLRdir.Data()),"recreate");

  fPlotsVBF->cd() ;  
  TH1D* PtTauTau = new TH1D ("PtTauTau", "", 100, 0, 400);
  TH1D* MTauTau = new TH1D ("MTauTau", "", 100, 0, 400);
  TH1D* EtaJet1 = new TH1D ("EtaJet1", "", 50, -5, 5);
  TH1D* EtaJet2 = new TH1D ("EtaJet2", "", 50, -5, 5);

  TH2D* Mjj_corr = new TH2D ("Mjj_corr", "", 50, 0, 800,50,0,800);
  TH1D* Mjj_res = new TH1D ("Mjj_res", "", 50, -1, 1);

  TH1D* Mjj_VBF_pass = new TH1D ("Mjj_VBF_pass", "", 40, 400, 2000);
  TH1D* Mjj_VBF = new TH1D ("Mjj_VBF", "", 40, 400, 2000);
  TH1D* Mjj_VBFMatched_pass = new TH1D ("Mjj_VBFMatched_pass", "", 40, 400, 2000);
  TH1D* Mjj_VBF_tot = new TH1D ("Mjj_VBF_tot", "", 40, 400, 2000);
  TH1D* Mjj_ditau = new TH1D ("Mjj_ditau", "",40 ,400, 2000);
  TH1D* Mjj_ditau_VBF = new TH1D ("Mjj_ditau_VBF", "",40 ,400, 2000);
  TH1D* Mjj_ditau_noVBF = new TH1D ("Mjj_ditau_noVBF", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBF = new TH1D ("Mjj_noditau_VBF", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbySeed = new TH1D ("Mjj_noditau_VBFbySeed", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbySel = new TH1D ("Mjj_noditau_VBFbySel", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbyBoth = new TH1D ("Mjj_noditau_VBFbyBoth", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFsel = new TH1D ("Mjj_noditau_VBFsel", "",40 ,400, 2000);
  TH1D* Mjj_or = new TH1D ("Mjj_or", "",40 ,400, 2000);

  TH1D* PtTau_VBF = new TH1D ("PtTau_VBF", "", 40, 20, 120);
  TH1D* PtTau_ditau = new TH1D ("PtTau_ditau", "", 40, 20, 120);
  TH1D* PtTau_or = new TH1D ("PtTau_or", "",40, 20, 120);
  TH1D* PtTau_ditau_noVBF = new TH1D ("PtTau_ditau_noVBF", "",40, 20, 120);
  TH1D* PtTau_noditau_VBF = new TH1D ("PtTau_noditau_VBF", "",40, 20, 120);

  TGraphAsymmErrors* Mjj_TurnOn = new  TGraphAsymmErrors();  
  Mjj_TurnOn->SetName("Mjj_TurnOn_VBF");
  TGraphAsymmErrors* Mjj_TurnOnMatched = new  TGraphAsymmErrors();  
  Mjj_TurnOnMatched->SetName("Mjj_TurnOn_VBFMatched");

  TFile *fPlotsDiTauPt = TFile::Open(Form("%sVBFHTauTauSignal_DiTauPtseed.root", LLRdir.Data()),"recreate");
  TH1D* PtTauPair_DiTauPt = new TH1D ("PtTauPair_DiTauPt", "", 40, 0, 400); 
  TH1D* PtTauPair_DiTau = new TH1D ("PtTauPair_DiTau", "",40 ,0, 400);
  TH1D* PtTauPair_noDiTau_DiTauPt = new TH1D ("PtTauPair_noDiTau_DiTauPt", "",40 ,0, 400);
  TH1D* PtTauPair_or = new TH1D ("PtTauPair_or", "",40, 0, 400);
  TH1D* PtTauPair_DiTau_noDiTauPt = new TH1D ("PtTauPair_DiTau_noDiTauPt", "",40, 0, 400);
  TH1D* PtTauPair_noDiTau_BoostbySeed = new TH1D ("PtTauPair_noDiTau_BoostbySeed", "",40,0 ,400);
  TH1D* PtTauPair_noDiTau_BoostbySel = new TH1D ("PtTauPair_noDiTau_BoostbySel", "",40,0 ,400);
  TH1D* PtTauPair_noDiTau_BoostbyBoth = new TH1D ("PtTauPair_noDiTau_BoostbyBoth", "",40,0 ,400);

  TH1D* BPtTau_DiTauPair = new TH1D ("BPtTau_DiTauPair", "", 40, 20, 120);
  TH1D* BPtTau_DiTau = new TH1D ("BPtTau_DiTau", "", 40, 20, 120);
  TH1D* BPtTau_or = new TH1D ("BPtTau_or", "",40, 20, 120);
  TH1D* BPtTau_DiTau_noDiTauPair = new TH1D ("BPtTau_DiTau_noDiTauPair", "",40, 20, 120);
  TH1D* BPtTau_noDiTau_DiTauPair = new TH1D ("BPtTau_noDiTau_DiTauPair", "",40, 20, 120);

  file->cd();
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

  
  std::vector <TLorentzVector> OffTau;
  std::vector <TLorentzVector> OffJet;
  std::vector <TLorentzVector> OffJetNoOverlap;

   
  std::vector<double> ptTau_pass; 
  std::vector<double> ptJet_pass;
  std::vector<object> tau;
  std::vector<object> tauNoOverlap;
  std::vector<object> tau25noOverlap;
  std::vector<object> jet30;   
  std::vector<object> jet30_sortByDeltaR;   
  std::vector<object> jet30noOverlap;   

  std::vector< std::tuple<double,int,int> > et_ditau_pass; //et of tau pair  
  std::vector< std::tuple<double,int,int> > m_ditau_pass; //m of tau pair  
  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int> > mjj_off_passditau; //VBF
  std::vector< std::tuple<double,int,int> > mjj_off_passVBF; //VBF
  std::vector< std::tuple<double,int,int> > ptpair_off; //ditauPt
  bool MatchingTau = false; 
  bool MatchingJet = false; 
  bool L1_DoubleIsoTau25er_Jet50= false;
  bool L1_DoubleIsoTau32er= false;
  bool L1_DoubleIsoTau25er_PtTauTau70 = false;
  bool L1_DoubleJet_90_30_Mj30j30_620 = false;
  
  
  double N_selXtrigger = 0;
  double N_seedXtrigger = 0;
  double N_selDiTauPt = 0;
  double N_seedDiTauPt = 0;
  double N_selDiTau = 0;
  double N_seedDiTau = 0;
  double N_seedXtrigger_noseedDiTau = 0;
  double N_seedDiTauPt_noseedDiTau = 0;
  double N_selDiTau_noVBF = 0;
  double N_selDiTau_VBF = 0;
  double N_noDiTau_VBF = 0;
  double N_DiTau_noDiTauPair = 0;
  double N_noDiTau_DiTauPair = 0;
  double N_DiTau_DiTauPair = 0;


  long int nEvents = tInput->GetEntries(); 

  for (long int iEv = 0; iEv < nEvents; iEv++){
    OffTau.clear();
    OffJetNoOverlap.clear();
    OffJet.clear();
      
    jet30.clear();
    jet30_sortByDeltaR.clear();
    jet30noOverlap.clear();
    tauNoOverlap.clear();
    tau25noOverlap.clear();
    mjj_pass.clear();
    mjj_off_passditau.clear();
    mjj_off_passVBF.clear();
    et_ditau_pass.clear();
    m_ditau_pass.clear();
    ptpair_off.clear();
    tau.clear(); 
    MatchingTau = false;
    MatchingJet = false;
    L1_DoubleIsoTau25er_PtTauTau70 =false;
    L1_DoubleIsoTau25er_Jet50= false;
    L1_DoubleIsoTau32er= false;
    L1_DoubleJet_90_30_Mj30j30_620 = false;


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
	if(tlv_tau.Pt()>pt_tau_cut && fabs(tlv_tau.Eta())<2.1 && ((tauID->at(i)>>5)&1) && ((tauID->at(i)>>7)&1) && ((tauID->at(i)>>5)&11)){
	  for(int iTau=0; iTau<stage2_tauN; iTau++){//matching
	    if(stage2_tauEt->at(iTau)>8){
	      TLorentzVector tlv_L1tau;
	      tlv_L1tau.SetPtEtaPhiM(stage2_tauEt->at(iTau),
				     stage2_tauEta->at(iTau),
				     stage2_tauPhi->at(iTau),
				     0.);        
	      if(tlv_tau.DeltaR(tlv_L1tau)<0.3){
		MatchingTau = true;
		break;	     
	      }
	    }	
	  }
	  if (MatchingTau = true) OffTau.push_back(tlv_tau);
	}
      }
    }
    std::sort(OffTau.begin(),OffTau.end(),SortByPt);
    if(OffTau.size()>=2){ //this holds for all the selections under study     
      TLorentzVector tlv_tauPair = OffTau[0]+OffTau[1];
      for(int j=0; j<JetsNumber;j++){
	TLorentzVector tlv_jet;
	tlv_jet.SetPxPyPzE(
			   jets_px->at(j),
			   jets_py->at(j),
			   jets_pz->at(j),
			   jets_e->at(j)
			   ) ;
	
	if(tlv_jet.Pt()>pt_jet_cut){ 
	  for(int iJet=0; iJet<stage2_jetN; iJet++){//matching with L1
	    if(stage2_jetEt->at(iJet)>8){
	      TLorentzVector tlv_L1jet;
	      tlv_L1jet.SetPtEtaPhiM(stage2_jetEt->at(iJet),
				     stage2_jetEta->at(iJet),
				     stage2_jetPhi->at(iJet),
				     0.);        
	      if(tlv_jet.DeltaR(tlv_L1jet)<0.3){
		MatchingJet = true;
		break;	     
	      }
	    }	
	  }
	  if(MatchingJet) OffJet.push_back(tlv_jet);
	  if(tlv_jet.DeltaR(OffTau[0])>0.3 && tlv_jet.DeltaR(OffTau[1])>0.3){ //overlap removal offline
	    if (MatchingJet)  OffJetNoOverlap.push_back(tlv_jet);
	  }
	}
      }
      std::sort(OffJetNoOverlap.begin(),OffJetNoOverlap.end(),SortByPt);
      std::sort(OffJet.begin(),OffJet.end(),SortByPt);
    
      ///// now all the offline vectors are filled
       
      //L1 objects
      

	for (long int iL1 = 0; iL1 < stage2_tauN; iL1++){ //loop on taus
	  // selections
	  double tauEta  = stage2_tauEta->at(iL1);
	  double tauIso  = stage2_tauIso->at(iL1);
	  if(tauIso>0.5 && fabs(tauEta)<2.1) tau.push_back(object(stage2_tauEt->at(iL1),stage2_tauEta->at(iL1),stage2_tauPhi->at(iL1),stage2_tauIso->at(iL1))) ;
	}
	
	for (long int iL1 = 0; iL1 < stage2_jetN; iL1++){ //loop on jets
	  // selections
	  double jetPt  = (*stage2_jetEt)[iL1];
	  ptJet_pass.push_back (jetPt);
	  if(jetPt>30.)	jet30.push_back(object(stage2_jetEt->at(iL1),stage2_jetEta->at(iL1),stage2_jetPhi->at(iL1),-999)) ;
	  if(jetPt>30.)	jet30_sortByDeltaR.push_back(object(stage2_jetEt->at(iL1),stage2_jetEta->at(iL1),stage2_jetPhi->at(iL1),-999)) ;
	}
	
	std::sort (jet30.begin(),jet30.end());
	
	std::sort (tau.begin(),tau.end());
	
	//overap removal L1
       	tauNoOverlap.push_back(object(tau[0].Et(),tau[0].Eta(),tau[0].Phi(),tau[0].Iso())) ;      
	tau25noOverlap.push_back(object(tau[0].Et(),tau[0].Eta(),tau[0].Phi(),tau[0].Iso())) ;      
	for (int iTau =1;iTau<tau.size();iTau++){
	  if (tau[iTau].DeltaR(tau[0])>0.2) tauNoOverlap.push_back(object(tau[iTau].Et(),tau[iTau].Eta(),tau[iTau].Phi(),tau[iTau].Iso())) ;      
	  if (tau[iTau].DeltaR(tau[0])>0.2 && tau[iTau].Et()>25) tau25noOverlap.push_back(object(tau[iTau].Et(),tau[iTau].Eta(),tau[iTau].Phi(),tau[iTau].Iso())) ;      
	}
	for (int iJet =0;iJet<jet30.size();iJet++){
	  if ((jet30[iJet].DeltaR(tau[0])>0.2)&&(jet30[iJet].DeltaR(tau[1])>0.2)) jet30noOverlap.push_back(object(jet30[iJet].Et(),jet30[iJet].Eta(),jet30[iJet].Phi(),-999)) ;      
	}
       
	std::sort (jet30noOverlap.begin(),jet30noOverlap.end());//not necessary
	std::sort (tauNoOverlap.begin(),tauNoOverlap.end());//not necessary
	std::sort (tau25noOverlap.begin(),tau25noOverlap.end());//not necessary
	
	//matched L1 jets
	if(jet30_sortByDeltaR.size()>1 && OffJetNoOverlap.size()>1){
	  object::SortByDeltaR DeltaRL1JetOffJet1(OffJetNoOverlap[0].Phi(),OffJetNoOverlap[0].Eta());
	  object::SortByDeltaR DeltaRL1JetOffJet2(OffJetNoOverlap[1].Phi(),OffJetNoOverlap[1].Eta());	  
	  
	  std::sort(jet30_sortByDeltaR.begin(),jet30_sortByDeltaR.end(),DeltaRL1JetOffJet1);
	  std::sort(jet30_sortByDeltaR.begin()+1,jet30_sortByDeltaR.end(),DeltaRL1JetOffJet2);
	  }
	
	
	//ditau+ pttautau
	if (tau25noOverlap.size() >= 2){
	  int Ntau = tau25noOverlap.size();
	  if (Ntau > 5) Ntau=5; 
	  for (int iTau = 0; iTau <Ntau; iTau++){      
	    for (int kTau = 0; kTau <Ntau; kTau++){      
	      if (kTau!=iTau) {
		TLorentzVector itau;
		itau.SetPtEtaPhiM(
				  tau25noOverlap[iTau].Et(),
				  tau25noOverlap[iTau].Eta(),
				  tau25noOverlap[iTau].Phi(),
				  0.);
		TLorentzVector ktau;
		ktau.SetPtEtaPhiM(
				  tau25noOverlap[kTau].Et(),
				  tau25noOverlap[kTau].Eta(),
				  tau25noOverlap[kTau].Phi(),
				  0.);
		TLorentzVector tauPair = itau+ktau;
		et_ditau_pass.push_back(make_tuple(tauPair.Et(),iTau,kTau));
		m_ditau_pass.push_back(make_tuple(tauPair.M(),iTau,kTau));
	      }
	    }
	  }
	  std::sort(et_ditau_pass.begin(),et_ditau_pass.end());
	  std::sort(m_ditau_pass.begin(),m_ditau_pass.end());	       
	}
   
	     
      //VBF seed
      if (jet30.size() >= 2){
	for (int iJet = 0; iJet <jet30.size(); iJet++){      
	  for (int kJet = 0; kJet <jet30.size(); kJet++){      
	    if (kJet!=iJet) {
	     	      TLorentzVector ijet;
	      ijet.SetPtEtaPhiM(
				jet30[iJet].Et(),
				jet30[iJet].Eta(),
				jet30[iJet].Phi(),
				0.);
		TLorentzVector kjet;
		kjet.SetPtEtaPhiM(
				  jet30[kJet].Et(),
				  jet30[kJet].Eta(),
				  jet30[kJet].Phi(),
				  0.);
		TLorentzVector jetPair = ijet+kjet;
		mjj_pass.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    }
	    
	  }
	  
	}
	std::sort(mjj_pass.begin(),mjj_pass.end());
	if(std::get<0>(*(mjj_pass.rbegin()))>620 && jet30[0].Et()>90) L1_DoubleJet_90_30_Mj30j30_620 = true;      
      }
      //// all the L1 objects are filled

      if(tauNoOverlap.size()>=2){
	//DiTau seed
	if (tauNoOverlap[1].Et()>32) L1_DoubleIsoTau32er = true;
	//Xtrigger seed	       
	if(jet30noOverlap.size()>0){
	  if (tauNoOverlap[1].Et()>25 && jet30noOverlap[0].Et()>50) L1_DoubleIsoTau25er_Jet50 = true;
	}    
	
	//DiTau + PtTauTau seed 
	if (tau25noOverlap[1].Et()>25){
	  if(std::get<0>(*(et_ditau_pass.rbegin()))>70) L1_DoubleIsoTau25er_PtTauTau70 = true;  
	}
	     
      }
      
      
      //acceptance       
      bool selDiTau = false;
      bool selXtrigger = false;
      bool selDiTauPt = false;
      bool selVBF = false;
      bool selDiTau_forVBF = false;
      
      
      
      //plots
      TLorentzVector tlv_jetPair;
      TLorentzVector tlv_L1jetPair;
      
      if(OffJetNoOverlap.size()>=2){ 
	tlv_jetPair = OffJetNoOverlap[0] + OffJetNoOverlap[1];
	for (int iJet = 0; iJet <OffJetNoOverlap.size(); iJet++){ 
	  if(OffJetNoOverlap[iJet].Pt()<30) continue;     
	  for (int kJet = 0; kJet <OffJetNoOverlap.size(); kJet++){      
	  if(OffJetNoOverlap[kJet].Pt()<30) continue; 
	  if (OffJetNoOverlap[kJet].Pt()<30 ||OffJetNoOverlap[iJet].Pt()<30) cout<<"argh"<<endl;
	    if (kJet!=iJet) {
	      TLorentzVector ijet;
	      ijet.SetPtEtaPhiM(
				OffJetNoOverlap[iJet].Et(),
				OffJetNoOverlap[iJet].Eta(),
				OffJetNoOverlap[iJet].Phi(),
				0.);
	      TLorentzVector kjet;
	      kjet.SetPtEtaPhiM(
				OffJetNoOverlap[kJet].Et(),
				OffJetNoOverlap[kJet].Eta(),
				OffJetNoOverlap[kJet].Phi(),
				0.);
	      TLorentzVector jetPair = ijet+kjet;
	      mjj_off_passVBF.push_back(make_tuple(jetPair.M(),iJet,kJet));
	  if (OffJetNoOverlap[kJet].Pt()>35 && OffJetNoOverlap[iJet].Pt()>35) mjj_off_passditau.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    }
	    
	  }
	}
	std::sort(mjj_off_passditau.begin(),mjj_off_passditau.end());
	std::sort(mjj_off_passVBF.begin(),mjj_off_passVBF.end());
	//for acceptance L1_DoubleJet_90_30_Mj30j30_620 //VBF
	if(OffJet.size()>1){
	  if(((OffJet[0].Pt()>100)||(OffTau[0].Pt()>100)) && OffJet[1].Pt()>30 &&  OffTau[1].Pt()>20){
	    selVBF = true;
	    if( L1_DoubleJet_90_30_Mj30j30_620 ) {
	      Mjj_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));	
	      PtTau_VBF->Fill(OffTau[1].Pt());
	    }
	    if (L1_DoubleIsoTau32er) Mjj_ditau_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));	    
	  }
	  if(OffJet[1].Pt()>35 && OffTau[1].Pt()>35 && std::get<0>(*(mjj_off_passditau.rbegin()))>400){
	    selDiTau_forVBF = true;
	    if(L1_DoubleIsoTau32er ){
	      PtTau_ditau->Fill(OffTau[1].Pt());
	      Mjj_ditau ->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
	      if(!(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF &&L1_DoubleJet_90_30_Mj30j30_620 )){
		N_selDiTau_noVBF += 1;	     
		Mjj_ditau_noVBF->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
		PtTau_ditau_noVBF->Fill(OffTau[1].Pt());
	      }
	    }	            
	  }
	  if(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF &&L1_DoubleJet_90_30_Mj30j30_620){
	    if(!(selDiTau_forVBF && L1_DoubleIsoTau32er)){
	      Mjj_noditau_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	      PtTau_noditau_VBF->Fill(OffTau[1].Pt())	;
	      N_noDiTau_VBF += 1;
	      //		cout<<"tau 1 pt = "<< OffTau[0].Pt()<<"; tau 2 pt = "<<OffTau[1].Pt()<<endl;	     
	      //		cout<<"jet 1 pt = "<< OffJet[0].Pt()<<"; jet 2 pt = "<<OffJet[1].Pt()<<endl;	     	
	    }
	    if(!(L1_DoubleIsoTau32er) && selDiTau_forVBF)  Mjj_noditau_VBFbySeed ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    if(!(selDiTau_forVBF) && !(L1_DoubleIsoTau32er)) Mjj_noditau_VBFbyBoth ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    if(!(selDiTau_forVBF) && L1_DoubleIsoTau32er) Mjj_noditau_VBFbySel ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    
	    if(!L1_DoubleIsoTau32er)      Mjj_noditau_VBFsel ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	  }

	  if((std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620) && (selDiTau_forVBF && L1_DoubleIsoTau32er)){
	    N_selDiTau_VBF+=1;	  
	  }
	  if((std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620) || (selDiTau_forVBF && L1_DoubleIsoTau32er)){	 
	    if((selDiTau_forVBF && L1_DoubleIsoTau32er)&&!(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620)) {
	      Mjj_or->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
	    }else{
	      Mjj_or->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    }
	    PtTau_or->Fill(OffTau[1].Pt());
	  }
          if(jet30_sortByDeltaR.size()>1){
	    TLorentzVector MatchL1jet1;
	    MatchL1jet1.SetPtEtaPhiM(
				     jet30_sortByDeltaR[0].Et(),
				     jet30_sortByDeltaR[0].Eta(),
				     jet30_sortByDeltaR[0].Phi(),
				     0.);
	    TLorentzVector MatchL1jet2;
	    MatchL1jet2.SetPtEtaPhiM(
				     jet30_sortByDeltaR[1].Et(),
				     jet30_sortByDeltaR[1].Eta(),
				     jet30_sortByDeltaR[1].Phi(),
				     0.);
	    
	    if(MatchL1jet1.DeltaR(OffJetNoOverlap[0])<0.3 && MatchL1jet2.DeltaR(OffJetNoOverlap[1])<0.3){
	      tlv_L1jetPair = MatchL1jet1+MatchL1jet2;
	      
	      Mjj_corr->Fill(tlv_L1jetPair.M(),tlv_jetPair.M());	  
	      Mjj_res->Fill((tlv_jetPair.M()-tlv_L1jetPair.M())/tlv_jetPair.M());	   

	      //VBF selection for turnon 
	      
	       if(((OffJet[0].Pt()>120)||(OffTau[0].Pt()>120)) && OffJet[1].Pt()>40 && OffTau[1].Pt()>20){
		
		Mjj_VBF_tot ->Fill(tlv_jetPair.M());
		if( L1_DoubleJet_90_30_Mj30j30_620 ){
		  Mjj_VBF_pass ->Fill(tlv_jetPair.M());
		}
		
		if(MatchL1jet1.Et()>90 && tlv_L1jetPair.M()>620)  Mjj_VBFMatched_pass ->Fill(tlv_jetPair.M());	     
	      }
	    }
	  }
	}
        
      if(OffJetNoOverlap.size()>0){	  


	//for acceptance L1_DoubleIsoTau25er_Jet50 //Xtrigger 
		
	if(OffJetNoOverlap[0].Pt()>60 && OffTau[0].Pt()>28 && OffTau[1].Pt()>28){
	  selXtrigger = true;
	  if (L1_DoubleIsoTau25er_Jet50) N_seedXtrigger +=1;		
	  N_selXtrigger+=1;
	}      
	if(selXtrigger &&L1_DoubleIsoTau25er_Jet50){
	  if(!(selDiTau &&L1_DoubleIsoTau32er)) N_seedXtrigger_noseedDiTau +=1;
	}
	int Ntau1 = 0;
	int Ntau2 = 0;
	//for acceptance L1_DoubleIsoTau25er_PtTauTau70 //DiTauPt 
	for (int iTau = 0; iTau <OffTau.size(); iTau++){ 
	  if(OffTau[iTau].Pt()<28) continue;
	  Ntau1++;
	  if(Ntau1>5) continue;
	  for (int kTau = 0; kTau <OffTau.size(); kTau++){      
	    if(OffTau[kTau].Pt()<28) continue; 
	    Ntau2++;
	    if(Ntau2>5) continue;
	    if (OffTau[kTau].Pt()<28 ||OffTau[iTau].Pt()<28) cout<<"argh"<<endl;
	    if (kTau!=iTau) {
	      
	      TLorentzVector itau;
	      itau.SetPtEtaPhiM(
				OffTau[iTau].Et(),
				OffTau[iTau].Eta(),
				OffTau[iTau].Phi(),
				0.);
	      TLorentzVector ktau;
	      ktau.SetPtEtaPhiM(
				OffTau[kTau].Et(),
				OffTau[kTau].Eta(),
				OffTau[kTau].Phi(),
				0.);
	      TLorentzVector tauPair = itau+ktau;
	      ptpair_off.push_back(make_tuple(tauPair.Pt(),iTau,kTau));
	    }    
	  }
	}
	std::sort(ptpair_off.begin(),ptpair_off.end());
	if(OffJetNoOverlap[0].Pt()>30 && OffTau[1].Pt()>35){
	  if(std::get<0>(*(ptpair_off.rbegin()))>100){
	  selDiTau = true;
	  N_selDiTau+=1;
	  if(L1_DoubleIsoTau32er){ 
	    N_seedDiTau +=1;

	  }	 
	  }     
	}
	if(selDiTau && L1_DoubleIsoTau32er){
	  PtTauPair_DiTau->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  BPtTau_DiTau->Fill(OffTau[1].Pt());
	}
	if(OffJetNoOverlap[0].Pt()>30 && OffTau[1].Pt()>28){
	  PtTauPair_DiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  BPtTau_DiTauPair->Fill(OffTau[1].Pt());
	  if(std::get<0>(*(ptpair_off.rbegin()))>100){
	    selDiTauPt = true;
	    if (L1_DoubleIsoTau25er_PtTauTau70) N_seedDiTauPt +=1;
	    N_selDiTauPt+=1;
	  }
	}
	if((selDiTauPt &&L1_DoubleIsoTau25er_PtTauTau70)||(selDiTau &&L1_DoubleIsoTau32er)) {
	  PtTauPair_or->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  BPtTau_or->Fill(OffTau[1].Pt());
	}
	if((selDiTauPt &&L1_DoubleIsoTau25er_PtTauTau70)&&(selDiTau &&L1_DoubleIsoTau32er)) N_DiTau_DiTauPair+=1;
	if(selDiTauPt &&L1_DoubleIsoTau25er_PtTauTau70){
	  if(!(selDiTau && L1_DoubleIsoTau32er)){
	    N_seedDiTauPt_noseedDiTau +=1;
	    N_noDiTau_DiTauPair+=1;
	    PtTauPair_noDiTau_DiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));	   
	    BPtTau_noDiTau_DiTauPair->Fill(OffTau[1].Pt());
	  }
	  if(!(selDiTau)&&L1_DoubleIsoTau32er)  PtTauPair_noDiTau_BoostbySel->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  if((selDiTau)&&!(L1_DoubleIsoTau32er))  PtTauPair_noDiTau_BoostbySeed->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  if(!(selDiTau)&&!(L1_DoubleIsoTau32er))  PtTauPair_noDiTau_BoostbyBoth->Fill(std::get<0>(*(ptpair_off.rbegin())));
	}
		     
	if((selDiTau &&L1_DoubleIsoTau32er)&&!(selDiTauPt &&L1_DoubleIsoTau25er_PtTauTau70))  {
	  PtTauPair_DiTau_noDiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));
	  BPtTau_DiTau_noDiTauPair->Fill(OffTau[1].Pt());
	  N_DiTau_noDiTauPair+=1;
	  
	}       
	}      
      }
      }
    }
   
    fPlotsVBF->cd();
  Mjj_TurnOn->BayesDivide(Mjj_VBF_pass, Mjj_VBF_tot);
  Mjj_TurnOn->Write();
  Mjj_TurnOnMatched->BayesDivide(Mjj_VBFMatched_pass, Mjj_VBF_tot);
  Mjj_TurnOnMatched->Write();
  
 
   fPlotsVBF->Write();
  fPlotsVBF->Close();
  fPlotsDiTauPt->cd();
  fPlotsDiTauPt->Write();
  fPlotsDiTauPt->Close();
  cout<<"Acceptance L1_DoubleIsoTau32er                                  "<<N_seedDiTau/N_selDiTau<<endl;
  cout<<"Events passing L1_DoubleIsoTau32er and selection                "<<N_seedDiTau<<endl;  
  cout<<"Acceptance L1_DoubleIsoTau25er_Jet50                            "<<N_seedXtrigger/N_selXtrigger<<endl;
  cout<<"Events passing L1_DoubleIsoTau25er_Jet50 and selection          "<<N_seedXtrigger<<endl;  
  cout<<"Gain                                                            "<<(N_seedXtrigger_noseedDiTau/N_selXtrigger)/(N_seedDiTau/N_selDiTau)<<endl;
  
  cout<<"Acceptance L1_DoubleIsoTau32er                                  "<<N_seedDiTau/N_selDiTau<<endl;
  cout<<"Events passing L1_DoubleIsoTau32er and selection                "<<N_seedDiTau<<endl;  
  cout<<"Acceptance L1_DoubleIsoTau25er_PtTauTau70                       "<<N_seedDiTauPt/N_selDiTauPt<<endl;
  cout<<"Events passing L1_DoubleIsoTau25er_PtTauTau70 and selection     "<<N_seedDiTauPt<<endl;  
  cout<<"Gain                                                            "<<(N_seedDiTauPt_noseedDiTau/N_selDiTauPt)/(N_seedDiTau/N_selDiTau)<<endl;
  cout<<endl;
  cout<<"VBF"<<endl;
  cout<<"Only DiTau "<<N_selDiTau_noVBF<<endl;
  cout<<"DiTau and VFB "<<N_selDiTau_VBF<<endl;
  cout<<"Only VBF "<<N_noDiTau_VBF<<endl;

  cout<<"DiTauPair"<<endl;
  cout<<"Only DiTau "<<N_DiTau_noDiTauPair<<endl;
  cout<<"DiTau and DiTauPair "<<N_DiTau_DiTauPair<<endl;
  cout<<"Only DiTauPair"<<N_noDiTau_DiTauPair<<endl;
  }



