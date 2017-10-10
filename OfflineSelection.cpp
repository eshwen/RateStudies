
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


  bool VBFHSignal= false;
  bool ggHSignal= false;
  bool DYJetsToLLSignal= false;
  bool HHbbTauTauSignal= false;

  
  if(atoi(argv[4])==1) VBFHSignal = true;
  if(atoi(argv[5])==1) ggHSignal = true;
  if(atoi(argv[6])==1) DYJetsToLLSignal = true;
  if(atoi(argv[7])==1) HHbbTauTauSignal= true; 
     
     
  int lowbin = 0;
  double Boost_OnlinePTtaucut = atof(argv[1]);
  double Boost_OnlinePTpaircut = atof(argv[2]);
  double Boost_OffPTtaucut = Boost_OnlinePTtaucut+3;
  double Boost_OffPTpaircut = Boost_OnlinePTpaircut+6;
  if(atof(argv[11])>1) Boost_OffPTpaircut = atof(argv[11]);
  if(atoi(argv[10])==1) Boost_OffPTpaircut = 100;
  double Boost_OffMjjcut = atof(argv[8]);
  double Boost_OffDeltaEtajjcut = atof(argv[9]);
  int Boost_OffjetsN = 1;
  if(Boost_OffMjjcut>1) Boost_OffjetsN = 2;
  
  double DiTau_OnlinePTtaucut = atof(argv[3]);
  double DiTau_OffPTtaucut = DiTau_OnlinePTtaucut+3;
  double DiTau_OffPTpaircut = 0;
  if(atoi(argv[10])==1) DiTau_OffPTpaircut = Boost_OffPTpaircut;
  double DiTau_OffMjjcut = Boost_OffMjjcut;
  double DiTau_OffDeltaEtajjcut = Boost_OffDeltaEtajjcut;
  int DiTau_OffjetsN = Boost_OffjetsN;

  double DiTauOR_OnlinePTtaucut = atof(argv[12]);
  double DiTauOR_OffPTtaucut = DiTauOR_OnlinePTtaucut+3;
  double DiTauOR_OffPTpaircut = 0;
  if(atoi(argv[10])==1) DiTauOR_OffPTpaircut = Boost_OffPTpaircut;
  double DiTauOR_OffMjjcut = Boost_OffMjjcut;
  double DiTauOR_OffDeltaEtajjcut = Boost_OffDeltaEtajjcut;
  int DiTauOR_OffjetsN = Boost_OffjetsN;
  
  int pt_tautau_cut = 0;
  int pt_tau_cut = 20;
  int pt_jet_cut = 30;
  TString LLRdir;
  if(VBFHSignal)          LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_VBFHToTauTau_-1Events_0Skipped/";
  if(ggHSignal)           LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_ggHToTauTau_-1Events_0Skipped/";
  if(HHbbTauTauSignal)    LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_ggRadionToHHTo2b2Tau_500_-1Events_0Skipped/";
  if(DYJetsToLLSignal)    LLRdir = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_DYJetsToLL_1500000Events_0Skipped/";
  TFile *file = TFile::Open(Form("%sHTauTauAnalysis_total.root", LLRdir.Data()),"read");
  if (file == NULL) cout<<"File not found"<<endl;

  cout<<"Signal sample: "<<LLRdir<<endl;
  
  cout<<"DiTau Seed: DoubleIsoTau"<<fixed << setprecision(0)<<DiTau_OnlinePTtaucut<<"er"<<endl;
  cout<<"BoostedDiTau Seed: DoubleIsoTau"<<fixed << setprecision(0)<<Boost_OnlinePTtaucut<<"er_PtTauTau"<<fixed << setprecision(0)<<Boost_OnlinePTpaircut<<endl;

  cout<<"Test OR of: DoubleIsoTau"<<fixed <<setprecision(0)<<DiTauOR_OnlinePTtaucut<<"er + DoubleIsoTau"<<Boost_OnlinePTtaucut<<"er_PtTauTau"<<fixed << setprecision(0)<<Boost_OnlinePTpaircut<<endl;

  cout<<endl;
  cout<<"DiTau OfflineSelection:"<<endl;
  cout<<"subleading tau PT>"<<fixed << setprecision(0)<<DiTau_OffPTtaucut<<" GeV"<<endl;
  cout<<"tau pair PT>"<<fixed << setprecision(0)<<DiTau_OffPTpaircut<<" GeV"<<endl;
  cout<<"at least "<<fixed << setprecision(0)<<DiTau_OffjetsN<<" jets with PT>30 GeV";
  if (DiTau_OffjetsN>1) cout<<" with mjj>"<<fixed << setprecision(0)<<DiTau_OffMjjcut<<" GeV, |DeltaEta(j1,j2)|>"<<fixed << setprecision(1)<<DiTau_OffDeltaEtajjcut;
  cout<<endl;

  cout<<endl;
  cout<<"Boost OfflineSelection:"<<endl;
  cout<<"subleading tau PT>"<<fixed << setprecision(0)<<Boost_OffPTtaucut<<" GeV"<<endl;
  cout<<"tau pair PT>"<<fixed << setprecision(0)<<Boost_OffPTpaircut<<" GeV"<<endl;
  cout<<"at least "<<fixed << setprecision(0)<<Boost_OffjetsN<<" jets with PT>30 GeV";
  if (Boost_OffjetsN>1) cout<<" with mjj>"<<fixed << setprecision(0)<<Boost_OffMjjcut<<" GeV, |DeltaEta(j1,j2)|>"<<fixed << setprecision(1)<<Boost_OffDeltaEtajjcut;
  cout<<endl;

  
  cout<<endl;
  cout<<"DiTau (in OR with boost) OfflineSelection:"<<endl;
  cout<<"subleading tau PT>"<<fixed << setprecision(0)<<DiTauOR_OffPTtaucut<<" GeV"<<endl;
  cout<<"tau pair PT>"<<fixed << setprecision(0)<<DiTauOR_OffPTpaircut<<" GeV"<<endl;
  cout<<"at least "<<fixed << setprecision(0)<<DiTauOR_OffjetsN<<" jets with PT>30 GeV";
  if (DiTauOR_OffjetsN>1) cout<<" with mjj>"<<fixed << setprecision(0)<<DiTauOR_OffMjjcut<<" GeV, |DeltaEta(j1,j2)|>"<<fixed << setprecision(1)<<DiTauOR_OffDeltaEtajjcut;
  cout<<endl;


  TFile *fPlotsVBF;
  if(VBFHSignal)       fPlotsVBF = TFile::Open(Form("%sVBFHTauTauSignal_VBFseed.root", LLRdir.Data()),"recreate");
  if(ggHSignal)        fPlotsVBF = TFile::Open(Form("%sggHTauTauSignal_VBFseed.root", LLRdir.Data()),"recreate");
  if(DYJetsToLLSignal) fPlotsVBF = TFile::Open(Form("%sDYJetsToLL_VBFseed.root", LLRdir.Data()),"recreate");
  if(HHbbTauTauSignal) fPlotsVBF = TFile::Open(Form("%sHbbTauTauSignal_VBFseed.root", LLRdir.Data()),"recreate");
  fPlotsVBF->cd() ;  

  TH1D* EtaJet1 = new TH1D ("EtaJet1", "", 50, -5, 5);
  TH1D* EtaJet2 = new TH1D ("EtaJet2", "", 50, -5, 5);
  TH1D* DeltaEta_VBF = new TH1D ("DeltaEta_VBF", "", 50, 0, 5);
  TH2D* Mjj_corr = new TH2D ("Mjj_corr", "", 50, 0, 800,50,0,800);
  TH1D* Mjj_res = new TH1D ("Mjj_res", "", 50, -1, 1);

  TH1D* Mjj_VBF_pass = new TH1D ("Mjj_VBF_pass", "", 40, 400, 2000);
  TH1D* Mjj_VBF = new TH1D ("Mjj_VBF", "", 100, 0, 2000);
  TH1D* Mjj_VBFMatched_pass = new TH1D ("Mjj_VBFMatched_pass", "", 40, 400, 2000);
  TH1D* Mjj_VBF_tot = new TH1D ("Mjj_VBF_tot", "", 40, 400, 2000);
  TH1D* Mjj_ditau = new TH1D ("Mjj_DiTau", "",40 ,400, 2000);
  TH1D* Mjj_ditau_VBF = new TH1D ("Mjj_DiTau_VBF", "",40 ,400, 2000);
  TH1D* Mjj_ditau_noVBF = new TH1D ("Mjj_DiTau_noVBF", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBF = new TH1D ("Mjj_noDiTau_VBF", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbySeed = new TH1D ("Mjj_noDiTau_NewbySeed", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbySel = new TH1D ("Mjj_noDiTau_NewbySel", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFbyBoth = new TH1D ("Mjj_noDiTau_NewbyBoth", "",40 ,400, 2000);
  TH1D* Mjj_noditau_VBFsel = new TH1D ("Mjj_noDiTau_VBFsel", "",40 ,400, 2000);
  TH1D* Mjj_or = new TH1D ("Mjj_or", "",40 ,400, 2000);
  TH1D* Mjj_and = new TH1D ("Mjj_and", "",40 ,400, 2000);
  
  TH1D* PtTau_VBF = new TH1D ("PtTau_VBF", "", 40, 20, 120);
  TH1D* PtTau_ditau = new TH1D ("PtTau_DiTau", "", 40, 20, 120);
  TH1D* PtTau_or = new TH1D ("PtTau_DiTau_VBF", "",40, 20, 120);
  TH1D* PtTau_and = new TH1D ("PtTau_and", "",40, 20, 120);
  TH1D* PtTau_ditau_noVBF = new TH1D ("PtTau_DiTau_noVBF", "",40, 20, 120);
  TH1D* PtTau_noditau_VBF = new TH1D ("PtTau_noDiTau_VBF", "",40, 20, 120);

  TGraphAsymmErrors* Mjj_TurnOn = new  TGraphAsymmErrors();  
  Mjj_TurnOn->SetName("Mjj_TurnOn_VBF");
  TGraphAsymmErrors* Mjj_TurnOnMatched = new  TGraphAsymmErrors();  
  Mjj_TurnOnMatched->SetName("Mjj_TurnOn_VBFMatched");

  TFile *fPlotsDiTauPt;
  if(VBFHSignal)      fPlotsDiTauPt = TFile::Open(Form("%sVBFHTauTauSignal_DiTauPtseed_%.0lf_%.0lf.root", LLRdir.Data(),Boost_OnlinePTtaucut,Boost_OnlinePTpaircut),"recreate");
  if(ggHSignal)       fPlotsDiTauPt = TFile::Open(Form("%sggHTauTauSignal_DiTauPtseed_%.0lf_%.0lf.root", LLRdir.Data(),Boost_OnlinePTtaucut,Boost_OnlinePTpaircut),"recreate");
  if(DYJetsToLLSignal)fPlotsDiTauPt = TFile::Open(Form("%sDYJetsToLL_DiTauPtseed_%.0lf_%.0lf.root", LLRdir.Data(),Boost_OnlinePTtaucut,Boost_OnlinePTpaircut),"recreate");
  if(HHbbTauTauSignal)fPlotsDiTauPt = TFile::Open(Form("%sHbbTauTauSignal_DiTauPtseed_%.0lf_%.0lf.root", LLRdir.Data(),Boost_OnlinePTtaucut,Boost_OnlinePTpaircut),"recreate");  
	    
  TH1D* PtTauPair_DiTauPt = new TH1D ("PtTauPair_DiTauPt", "", 40, lowbin, 400); 
  TH1D* PtTauPair_DiTau = new TH1D ("PtTauPair_DiTau", "",40 ,lowbin, 400);
  TH1D* PtTauPair_noDiTau_DiTauPt = new TH1D ("PtTauPair_noDiTau_DiTauPt", "",40 ,lowbin, 400);
  TH1D* PtTauPair_or = new TH1D ("PtTauPair_or", "",40,lowbin, 400);
  TH1D* PtTauPair_and = new TH1D ("PtTauPair_and", "",40,lowbin, 400);
  TH1D* PtTauPair_DiTau_noDiTauPt = new TH1D ("PtTauPair_DiTau_noDiTauPt", "",40, lowbin, 400);
  TH1D* PtTauPair_noDiTau_BoostbySeed = new TH1D ("PtTauPair_noDiTau_NewbySeed", "",40,lowbin ,400);
  TH1D* PtTauPair_noDiTau_BoostbySel = new TH1D ("PtTauPair_noDiTau_NewbySel", "",40,lowbin ,400);
  TH1D* PtTauPair_noDiTau_BoostbyBoth = new TH1D ("PtTauPair_noDiTau_NewbyBoth", "",40,lowbin ,400);

  TH1D* Boost_PtTau_DiTauPair = new TH1D ("Boost_PtTau_DiTauPt", "", 40, 20, 120);
  TH1D* Boost_PtTau_DiTau = new TH1D ("Boost_PtTau_DiTau", "", 40, 20, 120);
  TH1D* Boost_PtTau_and = new TH1D ("Boost_PtTau_and", "", 40, 20, 120);
  TH1D* Boost_PtTau_or = new TH1D ("Boost_PtTau_or", "",40, 20, 120);
  TH1D* Boost_PtTau_DiTau_noDiTauPair = new TH1D ("Boost_PtTau_DiTau_noDiTauPt", "",40, 20, 120);
  TH1D* Boost_PtTau_noDiTau_DiTauPair = new TH1D ("Boost_PtTau_noDiTau_DiTauPt", "",40, 20, 120);

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
  std::vector< std::tuple<double,int,int> > ptpair_off_passditau; //ditauPt
  std::vector< std::tuple<double,int,int> > ptpair_off_passditauOR; //ditauPt
  bool MatchingTau = false; 
  bool MatchingJet = false; 
  bool L1_DoubleIsoTau25er_Jet50= false;
  bool L1_DoubleIsoTauYYer= false;
  bool L1_DoubleIsoTauXXer= false;
  bool L1_DoubleIsoTauXXer_PtTauTauYY = false;
  bool L1_DoubleJet_90_30_Mj30j30_620 = false;
  
  
  double N_selXtrigger = 0;
  double N_seedXtrigger = 0;
  double N_selDiTauPt = 0;
  double N_seedOR = 0;
  double N_selDiTau = 0;
  double N_selOR = 0;
  double N_seedDiTau = 0;
  double N_seedDiTauPt = 0;
  double N_seedXtrigger_noseedDiTau = 0;
  double N_seedDiTauPt_noseedDiTau = 0;
  double N_selDiTau_noVBF = 0;
  double N_selDiTau_VBF = 0;
  double N_noDiTau_VBF = 0;
  double N_DiTau_noDiTauPair = 0;
  double N_noDiTau_DiTauPair = 0;
  double N_DiTau_DiTauPair = 0;
  double N_DiTau_noOR = 0;
  double N_noDiTau_OR = 0;
  double N_DiTau_OR = 0;


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
    ptpair_off_passditau.clear();
    ptpair_off_passditauOR.clear();
    tau.clear(); 
    MatchingTau = false;
    MatchingJet = false;
    L1_DoubleIsoTauXXer_PtTauTauYY =false;
    L1_DoubleIsoTau25er_Jet50= false;
    L1_DoubleIsoTauXXer= false;
    L1_DoubleIsoTauYYer= false;
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
	if(jetPt>30.)	jet30.push_back(object(stage2_jetEt->at(iL1),stage2_jetEta->at(iL1),stage2_jetPhi->at(iL1),-999)) ; //CHECK
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
       //	if(std::get<0>(*(mjj_pass.rbegin()))>630 && jet30[0].Et()>100) L1_DoubleJet_90_30_Mj30j30_620 = true;      //CHECK

      }

      //// all the L1 objects are filled

      if(tauNoOverlap.size()>=2){
	//DiTau seeds

	if (tauNoOverlap[1].Et()>DiTau_OnlinePTtaucut) L1_DoubleIsoTauXXer = true;
	if (tauNoOverlap[1].Et()>DiTauOR_OnlinePTtaucut) L1_DoubleIsoTauYYer = true;

	//Xtrigger seed	       
	if(jet30noOverlap.size()>0){
	  if (tauNoOverlap[1].Et()>25 && jet30noOverlap[0].Et()>50) L1_DoubleIsoTau25er_Jet50 = true;
	}    
	
	//DiTau + PtTauTau seed 
	if (tau25noOverlap[1].Et()>Boost_OnlinePTtaucut){
	  if(std::get<0>(*(et_ditau_pass.rbegin()))>Boost_OnlinePTpaircut) L1_DoubleIsoTauXXer_PtTauTauYY = true;  
	}
	     
      }
      
      
      //acceptance       
      bool selDiTau = false;
      bool selOR = false;
      bool DiTauORBoost = false;
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
				OffJetNoOverlap[iJet].Pt(),
				OffJetNoOverlap[iJet].Eta(),
				OffJetNoOverlap[iJet].Phi(),
				OffJetNoOverlap[iJet].M()
				);
	      TLorentzVector kjet;
	      kjet.SetPtEtaPhiM(
				OffJetNoOverlap[kJet].Pt(),
				OffJetNoOverlap[kJet].Eta(),
				OffJetNoOverlap[kJet].Phi(),
				OffJetNoOverlap[kJet].M()
				);
	      TLorentzVector jetPair = ijet+kjet;
	      mjj_off_passVBF.push_back(make_tuple(jetPair.M(),iJet,kJet));
	      if (OffJetNoOverlap[kJet].Pt()>DiTau_OffPTtaucut && OffJetNoOverlap[iJet].Pt()>DiTau_OffPTtaucut) mjj_off_passditau.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    }
	    
	  }
	}
	std::sort(mjj_off_passditau.begin(),mjj_off_passditau.end());
	std::sort(mjj_off_passVBF.begin(),mjj_off_passVBF.end());
	//for acceptance L1_DoubleJet_90_30_Mj30j30_620 //VBF

	if(OffJet.size()>1){
	   if(((OffJet[0].Pt()>100)||(OffTau[0].Pt()>100)) && OffJet[1].Pt()>30 &&  OffTau[1].Pt()>20){
	  //	     	  if(((OffJet[0].Pt()>110)||(OffTau[0].Pt()>110)) && OffJet[1].Pt()>35 &&  OffTau[1].Pt()>20){ //CHECK
	     Mjj_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));	
	     selVBF = true;
	    if( L1_DoubleJet_90_30_Mj30j30_620 ) {
	  
	      PtTau_VBF->Fill(OffTau[1].Pt());
	      DeltaEta_VBF->Fill(fabs(OffJet[0].Eta()-OffJet[1].Eta()));
	    }
	    if (L1_DoubleIsoTauXXer) Mjj_ditau_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));	    
	  }
	  
	  if(mjj_off_passditau.size()>0){ 
	    if(OffJet[1].Pt()>DiTau_OffPTtaucut && OffTau[1].Pt()>DiTau_OffPTtaucut && std::get<0>(*(mjj_off_passditau.rbegin()))>400){
	      selDiTau_forVBF = true;
	      if(L1_DoubleIsoTauXXer ){
		PtTau_ditau->Fill(OffTau[1].Pt());
	      
		Mjj_ditau ->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
		if(!(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF &&L1_DoubleJet_90_30_Mj30j30_620 )){
		  N_selDiTau_noVBF += 1;	     
		  Mjj_ditau_noVBF->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
		  PtTau_ditau_noVBF->Fill(OffTau[1].Pt());
		
		}
	      }	            
	    }
	  }
	  if(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF &&L1_DoubleJet_90_30_Mj30j30_620){
	    if(!(selDiTau_forVBF && L1_DoubleIsoTauXXer)){
	      Mjj_noditau_VBF ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	      PtTau_noditau_VBF->Fill(OffTau[1].Pt())	;
	      N_noDiTau_VBF += 1;
	    }
	    if(!(L1_DoubleIsoTauXXer) && selDiTau_forVBF)  Mjj_noditau_VBFbySeed ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    if(!(selDiTau_forVBF) && !(L1_DoubleIsoTauXXer)) Mjj_noditau_VBFbyBoth ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    if(!(selDiTau_forVBF) && L1_DoubleIsoTauXXer) Mjj_noditau_VBFbySel ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    
	    if(!L1_DoubleIsoTauXXer)      Mjj_noditau_VBFsel ->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	  }

	  if((std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620) && (selDiTau_forVBF && L1_DoubleIsoTauXXer)){
	    N_selDiTau_VBF+=1;
	    Mjj_and->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    PtTau_and->Fill(OffTau[1].Pt());
	  }
	  if((std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620) || (selDiTau_forVBF && L1_DoubleIsoTauXXer)){	 
	    if((selDiTau_forVBF && L1_DoubleIsoTauXXer)&&!(std::get<0>(*(mjj_off_passVBF.rbegin()))>700 && selVBF && L1_DoubleJet_90_30_Mj30j30_620)) {
	      Mjj_or->Fill(std::get<0>(*(mjj_off_passditau.rbegin())));
	    }else{
	      Mjj_or->Fill(std::get<0>(*(mjj_off_passVBF.rbegin())));
	    }
	    PtTau_or->Fill(OffTau[1].Pt());
	  }

          if(jet30_sortByDeltaR.size()>1 && OffJetNoOverlap.size()>1){
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
	    if(!(selDiTau &&L1_DoubleIsoTauXXer)) N_seedXtrigger_noseedDiTau +=1;
	  }
	  int Ntau = OffTau.size();
	  if (OffTau.size()> 5) Ntau = 5;
	  //for acceptance L1_DoubleIsoTau25er_PtTauTau70 //DiTauPt 
	  for (int iTau = 0; iTau <Ntau; iTau++){ 

	    if(OffTau[iTau].Pt()<Boost_OffPTtaucut) continue;
	    // if(OffTau[iTau].Pt()<28) continue;
	    
	    for (int kTau = 0; kTau <OffTau.size(); kTau++){      
	      if(OffTau[kTau].Pt()<Boost_OffPTtaucut) continue;
	      // if(OffTau[kTau].Pt()<28) continue; 
	      
	      if (OffTau[kTau].Pt()<Boost_OffPTtaucut ||OffTau[iTau].Pt()<Boost_OffPTtaucut) cout<<"argh"<<endl;
	      if (kTau!=iTau) {
	      
		TLorentzVector itau;
		itau.SetPtEtaPhiM(
				  OffTau[iTau].Et(),
				  OffTau[iTau].Eta(),
				  OffTau[iTau].Phi(),
				  OffTau[iTau].M());
			       
		TLorentzVector ktau;
		ktau.SetPtEtaPhiM(
				  OffTau[kTau].Et(),
				  OffTau[kTau].Eta(),
				  OffTau[kTau].Phi(),
				  OffTau[kTau].M()
				  );
		TLorentzVector tauPair = itau+ktau;
		ptpair_off.push_back(make_tuple(tauPair.Pt(),iTau,kTau));
		if (OffTau[kTau].Pt()>DiTau_OffPTtaucut && OffTau[iTau].Pt()>DiTau_OffPTtaucut) ptpair_off_passditau.push_back(make_tuple(tauPair.Pt(),iTau,kTau));
		if (OffTau[kTau].Pt()>DiTauOR_OffPTtaucut && OffTau[iTau].Pt()>DiTauOR_OffPTtaucut) ptpair_off_passditauOR.push_back(make_tuple(tauPair.Pt(),iTau,kTau));
		
	      }    
	    }
	  }

	  std::sort(ptpair_off.begin(),ptpair_off.end());
	  std::sort(ptpair_off_passditau.begin(),ptpair_off_passditau.end());
	  std::sort(ptpair_off_passditauOR.begin(),ptpair_off_passditauOR.end());
      
	  if(OffJetNoOverlap[0].Pt()>30 && OffTau[1].Pt()>DiTau_OffPTtaucut){
	    //if(ptpair_off_passditau.size()>0){
	    if(std::get<0>(*(ptpair_off_passditau.rbegin()))>DiTau_OffPTpaircut){
	      if(DiTau_OffjetsN==2){
		if(OffJetNoOverlap.size()>=DiTau_OffjetsN){
		  if(OffJetNoOverlap[1].Pt()>30){
		    double DeltaEtajjOff = fabs(OffJetNoOverlap[1].Eta() -OffJetNoOverlap[0].Eta());
		    double MjjOff = (OffJetNoOverlap[1] + OffJetNoOverlap[0]).M();
		    if(MjjOff>DiTau_OffMjjcut && DeltaEtajjOff>DiTau_OffDeltaEtajjcut){
		      selDiTau = true;
		      N_selDiTau+=1;
		      if(L1_DoubleIsoTauXXer){ 
			N_seedDiTau +=1;
		      }
		    }
		  }
		}
	      }else{
		selDiTau = true;
		N_selDiTau+=1;
		if(L1_DoubleIsoTauXXer){ 
		  N_seedDiTau +=1;
		}	 
	      }
	    }
	    //  }
	    
	  }

	  if(OffJetNoOverlap[0].Pt()>30 && OffTau[1].Pt()>DiTauOR_OffPTtaucut){
	    //if(ptpair_off_passditau.size()>0){
	    if(std::get<0>(*(ptpair_off_passditauOR.rbegin()))>DiTauOR_OffPTpaircut){
	      if(DiTauOR_OffjetsN==2){
		if(OffJetNoOverlap.size()>=DiTauOR_OffjetsN){
		  if(OffJetNoOverlap[1].Pt()>30){
		    double DeltaEtajjOff = fabs(OffJetNoOverlap[1].Eta() -OffJetNoOverlap[0].Eta());
		    double MjjOff = (OffJetNoOverlap[1] + OffJetNoOverlap[0]).M();
		    if(MjjOff>DiTauOR_OffMjjcut && DeltaEtajjOff>DiTauOR_OffDeltaEtajjcut){
		      selOR = true;
		      N_selOR+=1;
		      if(L1_DoubleIsoTauYYer){ 
			//  N_seedOR +=1;
		      }
		    }
		  }
		}
	      }else{
		selOR = true;
		N_selOR+=1;
		if(L1_DoubleIsoTauYYer){ 
		  // N_seedOR +=1;
		}	 
	      }
	    }
	    //  }
	    
	  }


	  
	  if(selDiTau && L1_DoubleIsoTauXXer){
	    // if(ptpair_off_passditau.size()>0)
	    PtTauPair_DiTau->Fill(std::get<0>(*(ptpair_off_passditau.rbegin())));
	    Boost_PtTau_DiTau->Fill(OffTau[1].Pt());
	  }
	  if(OffJetNoOverlap[0].Pt()>30 && OffTau[1].Pt()>Boost_OffPTtaucut){
	    
	    
	    if(Boost_OffjetsN==2){
	      if(OffJetNoOverlap.size()>=Boost_OffjetsN){
		if(OffJetNoOverlap[1].Pt()>30){
		  double DeltaEtajjOff = fabs(OffJetNoOverlap[1].Eta() - OffJetNoOverlap[0].Eta());
		  double MjjOff = (OffJetNoOverlap[1] + OffJetNoOverlap[0]).M();
		  if(MjjOff>Boost_OffMjjcut && DeltaEtajjOff>Boost_OffDeltaEtajjcut){
		    PtTauPair_DiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));
		    Boost_PtTau_DiTauPair->Fill(OffTau[1].Pt());
		    if(std::get<0>(*(ptpair_off.rbegin()))>Boost_OffPTpaircut){
		      selDiTauPt = true;
		      if (L1_DoubleIsoTauXXer_PtTauTauYY) N_seedDiTauPt +=1;
		      N_selDiTauPt+=1;
		    }
		  }
		}
		
	      }
	      
	      
	    }else{
	      PtTauPair_DiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));
	      Boost_PtTau_DiTauPair->Fill(OffTau[1].Pt());
	      if(std::get<0>(*(ptpair_off.rbegin()))>Boost_OffPTpaircut){
		selDiTauPt = true;
		if (L1_DoubleIsoTauXXer_PtTauTauYY) N_seedDiTauPt +=1;
		N_selDiTauPt+=1;
	      } 
	    }
		
	      
	  }
	   
	
	  if((selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY)||(selDiTau &&L1_DoubleIsoTauXXer)) {

	    if(!(selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY)&&(selDiTau &&L1_DoubleIsoTauXXer)) {
	      if(ptpair_off_passditau.size()>0)	PtTauPair_or->Fill(std::get<0>(*(ptpair_off_passditau.rbegin())));
	    }else{
	      if(ptpair_off.size()>0) PtTauPair_or->Fill(std::get<0>(*(ptpair_off.rbegin())));
	    }
	    Boost_PtTau_or->Fill(OffTau[1].Pt());
	  }
	  if((selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY)&&(selDiTau &&L1_DoubleIsoTauXXer)){
	    Boost_PtTau_and->Fill(OffTau[1].Pt());
	    PtTauPair_and->Fill(std::get<0>(*(ptpair_off_passditau.rbegin())));
	    N_DiTau_DiTauPair+=1;
	  }
	  if(selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY){
	    if(!(selDiTau && L1_DoubleIsoTauXXer)){
	      N_seedDiTauPt_noseedDiTau +=1;
	      N_noDiTau_DiTauPair+=1;
	      if(ptpair_off.size()>0)  PtTauPair_noDiTau_DiTauPt->Fill(std::get<0>(*(ptpair_off.rbegin())));	   
	      Boost_PtTau_noDiTau_DiTauPair->Fill(OffTau[1].Pt());
	    }
	    if(!(selDiTau)&&L1_DoubleIsoTauXXer)  PtTauPair_noDiTau_BoostbySel->Fill(std::get<0>(*(ptpair_off.rbegin())));
	    if((selDiTau)&&!(L1_DoubleIsoTauXXer))  PtTauPair_noDiTau_BoostbySeed->Fill(std::get<0>(*(ptpair_off.rbegin())));
	    if(!(selDiTau)&&!(L1_DoubleIsoTauXXer))  PtTauPair_noDiTau_BoostbyBoth->Fill(std::get<0>(*(ptpair_off.rbegin())));



	  }
		     
	  if((selDiTau &&L1_DoubleIsoTauXXer)&&!(selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY))  {
	    if(ptpair_off_passditau.size()>0) PtTauPair_DiTau_noDiTauPt->Fill(std::get<0>(*(ptpair_off_passditau.rbegin())));
	    Boost_PtTau_DiTau_noDiTauPair->Fill(OffTau[1].Pt());
	    N_DiTau_noDiTauPair+=1;

	    
	    
	    
	  }
	  ///////////////// OR
	  if ((selDiTauPt &&L1_DoubleIsoTauXXer_PtTauTauYY)||(selOR &&L1_DoubleIsoTauYYer)) {
	    DiTauORBoost = true;
	    N_seedOR+=1;
	  }
	  if((DiTauORBoost)||(selDiTau &&L1_DoubleIsoTauXXer)) {
	    
	    if((DiTauORBoost)&&(selDiTau &&L1_DoubleIsoTauXXer)){
	      N_DiTau_OR+=1;
	    }
	    if(DiTauORBoost){
	      if(!(selDiTau && L1_DoubleIsoTauXXer)){
		
		N_noDiTau_OR+=1;
	      }
	    }
	    
	    if((selDiTau &&L1_DoubleIsoTauXXer)&&!(DiTauORBoost))  {
	      N_DiTau_noOR+=1;
	    }    
	    
	    
	    
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
  std::cout.unsetf ( std::ios::floatfield ); 
  std::cout.precision(5);

  cout<<endl;
  cout<<"VBF"<<endl;
  cout<<"Only DiTau "<<N_selDiTau_noVBF<<endl;
  cout<<"DiTau and VFB "<<N_selDiTau_VBF<<endl;
  cout<<"Only VBF "<<N_noDiTau_VBF<<endl;
   cout<<"Gain "<<N_noDiTau_VBF/(N_selDiTau_noVBF+N_selDiTau_VBF)<<endl;  
  cout<<endl;
  cout<<"DiTauPair"<<endl;
  cout<<"Only DiTau "<<N_DiTau_noDiTauPair<<endl;
  cout<<"DiTau and DiTauPair "<<N_DiTau_DiTauPair<<endl;
  cout<<"Only DiTauPair"<<N_noDiTau_DiTauPair<<endl;
  cout<<"Gain "<<N_noDiTau_DiTauPair/(N_DiTau_noDiTauPair+N_DiTau_DiTauPair)<<endl;  
  cout<<endl;
  cout<<"N_seedDiTau\tN_DiTau_noDiTauPair\tN_DiTau_DiTauPair\tN_noDiTau_DiTauPair\tN_seedDiTauPt\tGain\tRatio"<<endl; 
  cout<<N_seedDiTau<<"\t"<<N_DiTau_noDiTauPair<<"\t"<<N_DiTau_DiTauPair<<"\t"<<N_noDiTau_DiTauPair<<"\t"<<N_seedDiTauPt<<"\t"<<N_noDiTau_DiTauPair/(N_DiTau_noDiTauPair+N_DiTau_DiTauPair)<<"\t"<<N_seedDiTauPt/N_seedDiTau<<endl; 
  cout<<endl;
  cout<<"DoubleIsoTau"<<DiTauOR_OnlinePTtaucut<<"er OR DoubleIsoTau"<<Boost_OnlinePTtaucut<<"er_PtTauTau"<<Boost_OnlinePTpaircut<<" wrt DoubleIsoTau"<<DiTau_OnlinePTtaucut<<"er"<<endl;
  cout<<"N_seedDiTau\tN_DiTau_noOR\tN_DiTau_OR\tN_noDiTau_OR\tN_seedOR\tGain\tRatio"<<endl;
  cout<<N_seedDiTau<<"\t"<<N_DiTau_noOR<<"\t"<<N_DiTau_OR<<"\t"<<N_noDiTau_OR<<"\t"<<N_seedOR<<"\t"<<N_noDiTau_OR/(N_DiTau_noOR+N_DiTau_OR)<<"\t"<<N_seedOR/N_seedDiTau<<endl; 
}




