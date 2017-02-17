// evaluates rate  - compile with c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
// and launch the executable 
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

using namespace std;




int main(int argc, char** argv){

  float nbStudiedRun = 96;
  float thisLumiRun = 7.6E32;
  float scaleToLumi = 1.8E34;
  float scale = 0.001*(nbStudiedRun*11245.6)*scaleToLumi/thisLumiRun;  

  

  cout << "Scale factor: " << scale << endl;
  //ZeroBias sample L1
  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias14Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events/";
  TFile *file = TFile::Open(Form("%sL1total.root", directory.Data()),"read");
  TTree * tInput = (TTree*) file->Get("L1Tree/L1Tree");
  //VBF sample
  //TString directory = "/data_CMS/cms/amendola/LLRHTauTauNtuples/HiggsTauTauOutput_VBFHToTauTau_-1Events_0Skipped_1487239220.05/";
  //TFile *file = TFile::Open(Form("%sVBFHTauTau_tau25tau25_jet30_tautau100.root", directory.Data()),"read");
  //TTree * tInput = (TTree*) file->Get("HTauTauTreeSel");

  bool isOffSel = false;
  TString fOutNameVBF;
  TString fOutNameXtrigger;
  TString fOutNameTaus;

  if(!isOffSel){
    fOutNameVBF = directory+"rateL1_VBF.root";
    fOutNameXtrigger = directory+"rateL1_xtrigger.root";
    fOutNameTaus = directory+"rateL1_taus.root";
  }
  
  
  
  Int_t lumi ;
  Int_t stage2_tauN ;
  std::vector<Float_t>* stage2_tauEt=0;
  std::vector<Float_t>* stage2_tauEta=0;
  std::vector<Float_t>* stage2_tauPhi=0;
  std::vector<int> *stage2_tauIso=0;  
  
  Int_t stage2_jetN ;
  std::vector<Float_t> *stage2_jetEt=0;
  std::vector<Float_t> *stage2_jetEta=0;
  std::vector<Float_t> *stage2_jetPhi=0;

  
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;
  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  
  tInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
  tInput ->SetBranchAddress("stage2_tauN", &stage2_tauN , &b_stage2_tauN );
  tInput ->SetBranchAddress("stage2_tauEta", &stage2_tauEta, &b_stage2_tauEta);
  tInput ->SetBranchAddress("stage2_tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
  tInput ->SetBranchAddress("stage2_tauEt", &stage2_tauEt, &b_stage2_tauEt);
  tInput ->SetBranchAddress("stage2_tauIso", &stage2_tauIso, &b_stage2_tauIso);
  tInput ->SetBranchAddress("stage2_jetN", &stage2_jetN , &b_stage2_jetN);
  tInput ->SetBranchAddress("stage2_jetEta", &stage2_jetEta, &b_stage2_jetEta);
  tInput ->SetBranchAddress("stage2_jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
  tInput ->SetBranchAddress("stage2_jetEt", &stage2_jetEt, &b_stage2_jetEt);
  
  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");
  TFile* fOutTaus = new TFile (Form("%s",fOutNameTaus.Data()), "recreate");
  TFile* fOutXtrigger = new TFile (Form("%s",fOutNameXtrigger.Data()), "recreate");
  
  fOutTaus->cd(); 
  //taus 
  TH1D* LeadTauPt_Pass = new TH1D ("LeadTauPt_Pass", "LeadTauPt_Pass", 50, 0, 100);
  TH1D* SingleTauPt_Pass = new TH1D ("SingleTauPt_Pass", "SingleTauPt_Pass", 50, 0, 100);
  TH1D* SubleadTauPt_Pass = new TH1D ("SubleadTauPt_Pass", "SubleadTauPt_Pass", 50, 0, 100);
  TH1D* SubleadTauPt_Pass_jet = new TH1D ("SubleadTauPt_Pass_jet", "SubleadTauPt_Pass", 50, 0, 100);
  TH1D* DiTauPt_Pass = new TH1D ("DiTauPt_Pass", "DiTauPt_Pass", 50, 0, 100);
  TH1D* RelativeRate_singleTau = new TH1D ("RelativeRate_singleTau", "Relative rate - single tau", 50, 0, 100);
  TH1D* RelativeRate_diTau = new TH1D ("RelativeRate_diTau", "Relative rate - double tau", 50, 0, 100);
  TH1D* RelativeRate_diTauSublead = new TH1D ("RelativeRate_diTauSublead", "Relative rate - double tau", 50, 0, 100);
  TH1D* Rate_diTauSublead = new TH1D ("Rate_diTauSublead", "rate - double tau", 50, 0, 100);
  TH1D* Ratio_diTauSublead = new TH1D ("Ratio_diTauSublead", "ratio - double tau", 50, 0, 100);
  TH1D* Rate_diTauSublead_jet = new TH1D ("Rate_diTauSublead_jet", "rate - double tau", 50, 0, 100);
  TH1D* Ratio_diTauSublead_jet = new TH1D ("Ratio_diTauSublead_jet", "ratio - double tau", 50, 0, 100);
  TH2D* DiTau2D_Pass = new TH2D ("DiTau2D_Pass", "double tau",  100, 0,400, 100, 25,125); 
  TH2D* Ratio_DiTau2D = new TH2D ("Ratio_DiTau2D", "Ratio - 2 taus", 100, 0,400, 100, 25,125); 
  TH2D* Rate_DiTau2D = new TH2D ("Rate_DiTau2D", "Rate - 2 taus", 100, 0,400, 100, 25,125); 
  TH1D* PtTauTau = new TH1D ("PtTauTau", "", 100, 0, 400);
  TH1D* MTauTau = new TH1D ("MTauTau", "", 100, 0, 400);
  
  fOutVBF->cd(); 
  //VBF
  TH2D* DiJet2D_Pass = new TH2D ("DiJet2D", "double jet const", 100, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D = new TH2D ("Ratio_DiJet2D", "Ratio - 2 jets, E_{T} > 30 GeV ", 100, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D = new TH2D ("Rate_DiJet2D", "Rate - 2 jets, E_{T} > 30 GeV ", 100, 0, 800, 100, 30, 130); 
  
  TH1D* ETA = new TH1D ("ETA", "", 1000, -4, 4);
  TH1D* ETA_strip1 = new TH1D ("ETA_strip400", "", 1000, -4, 4);
  TH1D* ETA_strip2 = new TH1D ("ETA_strip467", "", 1000, -4, 4);
  
  TH2D* DiJet2D_Pass_strips = new TH2D ("DiJet2D_strips", "check", 100, 0, 800, 100,30, 130);
  TH2D* DiJet2D_Pass_strips_nolead = new TH2D ("DiJet2D_strips_nolead", "check", 100, 0, 800, 100,30, 130);
  
  fOutXtrigger->cd(); 
  //Di-Tau+jet  
  TH2D* DiTauJet2D_Pass = new TH2D ("DiTauJet2D", "double tau + jet", 100, 0, 100, 100,0, 100);
  TH2D* Ratio_DiTauJet2D = new TH2D ("Ratio_DiTauJet2D", "Ratio - 2 taus, jet E_{T} > 30 GeV", 100, 0, 100, 100, 0, 100);
  TH2D* Rate_DiTauJet2D = new TH2D ("Rate_DiTauJet2D", "Rate - 2 taus, jet E_{T} > 30 GeV ", 100, 0, 100, 100, 0, 100); 

  TH2D* DiTauJet2D_Pass_jet = new TH2D ("DiTauJet2D_jet", "double tau + jet", 100, 0, 100, 100,30, 130);
  TH3D* DiTauJet3D_Pass_jet = new TH3D ("DiTauJet3D_jet", "double tau + jet", 100, 0, 100,100,0,100, 100,30, 130);
  TH2D* Ratio_DiTauJet2D_jet = new TH2D ("Ratio_DiTauJet2D_jet", "Ratio - 2 taus, jet E_{T} > 30 GeV", 100, 0, 100, 100, 30, 130);
  TH3D* Ratio_DiTauJet3D_jet = new TH3D ("Ratio_DiTauJet3D_jet", "Ratio - 2 taus 3D, jet E_{T} > 30 GeV", 100, 0, 100,100,0,100,100, 30, 130);
  TH2D* Rate_DiTauJet2D_jet = new TH2D ("Rate_DiTauJet2D_jet", "Rate - 2 taus, jet E_{T} > 30 GeV ", 100, 0, 100, 100, 30, 130); 
  TH3D* Rate_DiTauJet3D_jet = new TH3D ("Rate_DiTauJet3D_jet", "Rate - 2 taus 3D, jet E_{T} > 30 GeV ", 100, 0, 100, 100,0,100,100, 30, 130); 
  TH1D* IsoTau = new TH1D ("IsoTau", "", 20, 0, 10);


  ifstream PUFile("utils/PU_per_LS.txt");

  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);
      temp.ReplaceAll("5412,283171,","");
      int pos_coma = temp.First(",");
      TString LS_str(temp,pos_coma);
      //cout<<LS_str<<endl;
      TString Replacing = LS_str ;
      Replacing += ",";
      temp.ReplaceAll(Replacing.Data(),"");
      TString PU_str = temp;
      //cout<<PU_str<<endl;
      std::istringstream ss_LS(LS_str.Data());
      Int_t LS ;
      ss_LS >> LS;
      std::istringstream ss_PU(PU_str.Data());
      Float_t PU ;
      ss_PU >> PU;     
      PU_per_LS.insert(std::pair<Int_t,Float_t>(LS , PU ));

    }
  // analyze data    
  long int nEvents = tInput->GetEntries(); 
  int nEventsPass = 0;

  std::vector<double> ptTau_pass; 
  std::vector<double> ptJet_pass;
  std::vector<object> tau;
  std::vector<object> tauNoOverlap;
  std::vector<object> tau25noOverlap;
  std::vector<object> jet30;   
  std::vector<object> jet30noOverlap;   

  std::vector< std::tuple<double,int,int> > et_ditau_pass; //et of tau pair  
  std::vector< std::tuple<double,int,int> > m_ditau_pass; //m of tau pair  
  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF


  double singletau = 0;
  double ditau20 = 0;  
  double ditaujet = 0;  
  double ditau = 0;  
  double ditau20jet = 0;  

  bool L1_DoubleIsoTau25er_Jet50= false;
  bool L1_DoubleIsoTau32er= false;
  bool  L1_DoubleIsoTau25er_PtTauTau70 = false;

  
  double    acceptance_ditau =0;
  double     acceptance_ditaujet =0;
  double     acceptance_ditauPt=0;
  double     acceptance_ditau_ditaujet =0;
  double     acceptance_ditau_ditauPt =0;
  double      pureRate_ditaujet=0;
  double    pureRate_ditauPt =0;
  
  //nEvents = 1000000;
  // loop on all events
  for (long int iEv = 0; iEv < nEvents; iEv++){

    tInput->GetEntry(iEv);   
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;

    jet30.clear();
    jet30noOverlap.clear();
    tauNoOverlap.clear();
    tau25noOverlap.clear();
    mjj_pass.clear();
    et_ditau_pass.clear();
    m_ditau_pass.clear();
    tau.clear(); 
    
    if (!isOffSel){
      if(lumi<48 || lumi>221) continue;
      if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;
    }
    Float_t weight = PU_per_LS[48]/PU_per_LS[lumi];

    nEventsPass ++;
    L1_DoubleIsoTau25er_PtTauTau70 = false;
    L1_DoubleIsoTau25er_Jet50= false;
    L1_DoubleIsoTau32er= false;
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
     }

     std::sort (jet30.begin(),jet30.end());
     std::sort (tau.begin(),tau.end());

     //overap removal
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



     
     //DiTau+PtDiTau

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
       PtTauTau->Fill(std::get<0>(*(et_ditau_pass.rbegin())));
       MTauTau->Fill(std::get<0>(*(m_ditau_pass.rbegin())));
       DiTau2D_Pass -> Fill ( std::get<0>(*(et_ditau_pass.rbegin())),tau25noOverlap[1].Et(),weight);
     } else{
      DiTau2D_Pass -> Fill (-1,-1);
     }
       

	  //VBF
       
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
      
      DiJet2D_Pass -> Fill ( std::get<0>(*(mjj_pass.rbegin())),jet30[0].Et(),weight);


      
      
      /*      if(std::get<1>(*(mjj_pass.rbegin()))<=1 && std::get<2>(*(mjj_pass.rbegin()))<=1){
	      DiJet2D_Pass_strips->  Fill ( std::get<0>(*(mjj_pass.rbegin())),jet30[0].Et(),weight);
	      }
	      
	      //eta checks for strips
	      if(std::get<1>(*(mjj_pass.rbegin()))>0 && std::get<2>(*(mjj_pass.rbegin()))>0){
	      DiJet2D_Pass_strips_nolead->  Fill ( std::get<0>(*(mjj_pass.rbegin())),jet30[0].Et(),weight);
	      }
	      
	      if(std::get<0>(*(mjj_pass.rbegin()))>350 &&std::get<0>(*(mjj_pass.rbegin()))<400){
	      ETA->Fill(jet30[std::get<2>(*(mjj_pass.rbegin()))].Eta());
	      ETA->Fill(jet30[std::get<1>(*(mjj_pass.rbegin()))].Eta());
	      }
	      if(std::get<0>(*(mjj_pass.rbegin()))>400 &&std::get<0>(*(mjj_pass.rbegin()))<460){
	      ETA_strip1->Fill(jet30[std::get<2>(*(mjj_pass.rbegin()))].Eta());
	ETA_strip1->Fill(jet30[std::get<1>(*(mjj_pass.rbegin()))].Eta());
	}
	if(std::get<0>(*(mjj_pass.rbegin()))>460 &&std::get<0>(*(mjj_pass.rbegin()))<600){
	ETA_strip2->Fill(jet30[std::get<2>(*(mjj_pass.rbegin()))].Eta());
	ETA_strip2->Fill(jet30[std::get<1>(*(mjj_pass.rbegin()))].Eta());
	}*/
      
    } else {
      DiJet2D_Pass -> Fill (-1, -1);
    }

    
    
    if(!isOffSel){
      //taus and ditau+jet
      if (tauNoOverlap.size() >= 2 )
	{
	  LeadTauPt_Pass -> Fill(tauNoOverlap[0].Et(),weight); 
	  SingleTauPt_Pass -> Fill(tauNoOverlap[0].Et(),weight); 
	  SubleadTauPt_Pass -> Fill(tauNoOverlap[1].Et(),weight);
	  IsoTau -> Fill(tauNoOverlap[0].Iso());
	  DiTauPt_Pass -> Fill (tauNoOverlap[0].Et()+tauNoOverlap[1].Et(),weight);
	  if (jet30noOverlap.size()>0){
	    SubleadTauPt_Pass_jet -> Fill(tauNoOverlap[1].Et(),weight);
	    DiTauJet2D_Pass->Fill( tauNoOverlap[0].Et(),tauNoOverlap[1].Et(),weight );
	    DiTauJet2D_Pass_jet->Fill( tauNoOverlap[1].Et(),jet30noOverlap[0].Et(),weight);
	    //  DiTauJet3D_Pass_jet->Fill( tauNoOverlap[0].Et(),tauNoOverlap[1].Et(),jet30noOverlap[0].Et());
	  }else{
	  SubleadTauPt_Pass -> Fill(-1);
	  SubleadTauPt_Pass_jet -> Fill(-1);
	  DiTauJet2D_Pass->Fill( -1,-1 );
	  DiTauJet2D_Pass_jet->Fill( -1,-1);
	  //  DiTauJet3D_Pass_jet->Fill( -1,-1,-1 );
	}
      } 
    
    else if (tauNoOverlap.size() == 1){
      LeadTauPt_Pass -> Fill ( tauNoOverlap[0].Et(),weight ); 
      SingleTauPt_Pass -> Fill ( tauNoOverlap[0].Et(),weight ); 
      SubleadTauPt_Pass -> Fill ( -1 );
      SubleadTauPt_Pass_jet -> Fill ( -1 );
      DiTauPt_Pass -> Fill ( -1); 
      DiTauJet2D_Pass->Fill( -1,-1 );
      DiTauJet2D_Pass_jet->Fill( -1,-1 );
      //      DiTauJet3D_Pass_jet->Fill( -1,-1,-1 );
    } 
    
    else
      {
	LeadTauPt_Pass -> Fill ( -1 ); 
	SingleTauPt_Pass -> Fill ( -1 ); 
	SubleadTauPt_Pass -> Fill ( -1 );        
	SubleadTauPt_Pass_jet -> Fill ( -1 );
	DiTauPt_Pass -> Fill ( -1 ); 
	DiTauJet2D_Pass->Fill( -1,-1 );       
	DiTauJet2D_Pass_jet->Fill( -1,-1 );
	DiTauJet3D_Pass_jet->Fill( -1,-1,-1 );
      }       
      
    }



    //counters
    if(!isOffSel){

      if(tauNoOverlap.size()>=2){
	if (tauNoOverlap[1].Et()>32) L1_DoubleIsoTau32er = true;

	if(jet30noOverlap.size()>0){
	  if (tauNoOverlap[1].Et()>25 && jet30noOverlap[0].Et()>50) L1_DoubleIsoTau25er_Jet50 = true;
	}    

      }	
      if(tau25noOverlap.size()>=2){	
	if (tau25noOverlap[1].Et()>25){
	  if(std::get<0>(*(et_ditau_pass.rbegin()))>70) L1_DoubleIsoTau25er_PtTauTau70 = true;
	}
      }

      if (L1_DoubleIsoTau25er_Jet50 && !L1_DoubleIsoTau32er) pureRate_ditaujet +=1;
      if (L1_DoubleIsoTau25er_PtTauTau70 && !L1_DoubleIsoTau32er) pureRate_ditauPt +=1;

      
    }else{
      if (L1_DoubleIsoTau32er) acceptance_ditau +=1;
      if (L1_DoubleIsoTau25er_Jet50) acceptance_ditaujet +=1;
      if (L1_DoubleIsoTau25er_PtTauTau70) acceptance_ditauPt +=1;
      if (L1_DoubleIsoTau32er || L1_DoubleIsoTau25er_Jet50) acceptance_ditau_ditaujet +=1;
      if (L1_DoubleIsoTau32er || L1_DoubleIsoTau25er_PtTauTau70) acceptance_ditau_ditauPt +=1;
      
    } 
    
  }
    if(isOffSel){
      acceptance_ditau /= nEventsPass;
      acceptance_ditaujet /= nEventsPass;
      acceptance_ditauPt /= nEventsPass;
      acceptance_ditau_ditaujet /= nEventsPass;
      acceptance_ditau_ditauPt /= nEventsPass;
    }else{   

      pureRate_ditauPt *=scale/nEventsPass; 
    }



  // compute rate plots

  if(isOffSel){
    cout<<"Acceptance L1_DoubleIsoTau32er                                  "<<acceptance_ditau<<endl;
    cout<<"Acceptance L1_DoubleIsoTau25er_Jet50                            "<<acceptance_ditaujet<<endl;
    cout<<"Acceptance L1_DoubleIsoTau25er_PtTauTau70                       "<<acceptance_ditauPt<<endl;
    cout<<"Acceptance L1_DoubleIsoTau32er + L1_DoubleIsoTau25er_Jet50      "<<acceptance_ditau_ditaujet<<endl;
    cout<<"Acceptance L1_DoubleIsoTau32er + L1_DoubleIsoTau25er_PtTauTau70 "<<acceptance_ditau_ditauPt<<endl;
  }else{

    cout<<"events passed "<<pureRate_ditaujet<<endl;
    pureRate_ditaujet *=scale/nEventsPass; 
    cout<<"Pure Rate L1_DoubleIsoTau25er_Jet50        "<<pureRate_ditaujet<<" kHz"<<endl;
    cout<<"Pure Rate L1_DoubleIsoTau25er_PtTauTau70   "<<pureRate_ditauPt<<" kHz"<<endl;
  }
  cout << endl;
  cout << "Computing rates..." << endl; 
  //taus rates
  for (int i = 1; i <= 50; i++){

    int Tot_diTauSublead = SubleadTauPt_Pass -> Integral (1,51);    
    int Tot_singleTau = SingleTauPt_Pass -> Integral (1,51);
    int Tot_diTau = DiTauPt_Pass -> Integral (1,51);        

    double relRateSingleTau = 1.*(SingleTauPt_Pass->Integral(i, 51))/Tot_singleTau;
    double relRateDiTau = 1.*(DiTauPt_Pass->Integral(i, 51))/Tot_diTau;
    double binDiTauSublead = 1.*(SubleadTauPt_Pass->Integral(i, 51))/nEventsPass;
    double binDiTauSublead_jet = 1.*(SubleadTauPt_Pass_jet->Integral(i, 51))/nEventsPass;
    double relRateDiTauSublead =  1.*(SubleadTauPt_Pass->Integral(i, 51))/Tot_diTauSublead;

    RelativeRate_singleTau -> SetBinContent (i, relRateSingleTau);
    RelativeRate_diTau -> SetBinContent (i, relRateDiTau);        
    Ratio_diTauSublead -> SetBinContent (i, binDiTauSublead);        
    binDiTauSublead *=scale;
    Rate_diTauSublead -> SetBinContent (i, binDiTauSublead);        

    Ratio_diTauSublead_jet -> SetBinContent (i, binDiTauSublead_jet);        
    binDiTauSublead_jet *=scale;
    Rate_diTauSublead_jet -> SetBinContent (i, binDiTauSublead_jet);        


    RelativeRate_diTauSublead -> SetBinContent (i, relRateDiTauSublead);        
  }

  

  //VBF
  for (int i = 1; i <=DiJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Pass->GetNbinsY();j++ ){
      double binDiJet2D = 1.*(DiJet2D_Pass->Integral(i, DiJet2D_Pass->GetNbinsX()+1,j,DiJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      
    }
  }


  //ditau + ptditau
  for (int i = 1; i <=DiTau2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiTau2D_Pass->GetNbinsY();j++ ){
      double binDiTau2D = 1.*(DiTau2D_Pass->Integral(i, DiTau2D_Pass->GetNbinsX()+1,j,DiTau2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTau2D -> SetBinContent (i, j, binDiTau2D);        
      binDiTau2D *=scale;
      Rate_DiTau2D -> SetBinContent (i, j, binDiTau2D);        
      
    }
  }



  //diTau+jet
 for (int i = 1; i <=DiTauJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiTauJet2D_Pass->GetNbinsY();j++ ){
      double binDiTauJet2D = 1.*(DiTauJet2D_Pass->Integral(i, DiTauJet2D_Pass->GetNbinsX()+1,j,DiTauJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTauJet2D -> SetBinContent (i, j, binDiTauJet2D);        
      binDiTauJet2D *=scale;
      Rate_DiTauJet2D -> SetBinContent (i, j, binDiTauJet2D);        
      
    }
  }



 for (int i = 1; i <=DiTauJet2D_Pass_jet->GetNbinsX(); i++){
    for (int j = 1; j<= DiTauJet2D_Pass_jet->GetNbinsY();j++ ){
      double binDiTauJet2D_jet = 1.*(DiTauJet2D_Pass_jet->Integral(i, DiTauJet2D_Pass_jet->GetNbinsX()+1,j,DiTauJet2D_Pass_jet->GetNbinsY()+1))/nEventsPass;
      Ratio_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);        
      binDiTauJet2D_jet *=scale;
      Rate_DiTauJet2D_jet -> SetBinContent (i, j, binDiTauJet2D_jet);        
      
    }
  }


  /*   for (int i = 1; i <=DiTauJet3D_Pass_jet->GetNbinsX(); i++){
   cout<<"looping on i "<<i<<endl;
    for (int j = 1; j<= DiTauJet3D_Pass_jet->GetNbinsY();j++ ){
      for (int k = 1; k<= DiTauJet3D_Pass_jet->GetNbinsZ();k++ ){
	double binDiTauJet3D_jet = 1.*(DiTauJet3D_Pass_jet->Integral(i, DiTauJet3D_Pass_jet->GetNbinsX()+1,j,DiTauJet3D_Pass_jet->GetNbinsY()+1,k,DiTauJet3D_Pass_jet->GetNbinsZ()+1))/nEventsPass;
	Ratio_DiTauJet3D_jet -> SetBinContent (i, j, k, binDiTauJet3D_jet);        
	binDiTauJet3D_jet *=scale;
	Rate_DiTauJet3D_jet -> SetBinContent (i, j, k,binDiTauJet3D_jet);        
      }
   }

   }*/
 
 
 if(!isOffSel){
   fOutTaus -> Write();
   fOutVBF -> Write();
   fOutXtrigger -> Write();
 }
}

