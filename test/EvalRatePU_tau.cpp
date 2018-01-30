#include <sstream>
#include <fstream>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <regex>
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
#include "../interface/object.h"

using namespace std;

void appendFromFileList (TChain* chain, TString filename)
{
 
  std::ifstream infile(filename.Data());
  std::string line;
  while (std::getline(infile, line))
    {
      line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
      while (line.find(" ") != std::string::npos) line = line.erase(line.find(" "), 1); // remove white spaces
      while (line.find("\n") != std::string::npos) line = line.erase(line.find("\n"), 1); // remove new line characters
      while (line.find("\r") != std::string::npos) line = line.erase(line.find("\r"), 1); // remove carriage return characters
      if (!line.empty()) // skip empty lines
	chain->Add(line.c_str());
    }
  return;
}


int main(int argc, char** argv){

  int Mjj_cut = atoi(argv[1]);
  int jets_cut = atoi(argv[2]);
  int tau_cut = atoi(argv[3]);
  std::stringstream Mjj_ss;
  Mjj_ss << Mjj_cut;
  TString Mjj_str = Mjj_ss.str();
  std::stringstream jets_ss;
  jets_ss << jets_cut;
  string jets_str = jets_ss.str();
  std::stringstream tau_ss;
  tau_ss << tau_cut;
  string tau_str = tau_ss.str();
  TString directory = argv[4];
  TString fileList  = argv[5];
  TString PUperLumiFile = argv[6];
  ifstream PUFile(PUperLumiFile);
  int trigger_ID = atoi(argv[7]);
  bool emulated  =false;
  if(atoi(argv[8])==1) emulated = true; 
  int Run = atoi(argv[9]);
  int Fill = atoi(argv[10]);
  int maxevents = atoi(argv[11]);
  std::stringstream max_ss;
  max_ss << maxevents;
  TString max_str = max_ss.str();

  int njob = atoi(argv[12]);
  cout<<trigger_ID<<"trigger"<<endl; 

  cout <<Mjj_cut<<" "<<jets_cut<<" "<<tau_cut<<endl; 
  
  TChain * cInput;
  TChain * lInput;
  TChain * ugtInput; 

  if (emulated){
    cInput = new TChain ("l1UpgradeEmuTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
  }else{
    cInput= new TChain ("l1UpgradeTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
  }

  lInput = new TChain ("l1EventTree/L1EventTree");
  appendFromFileList(lInput, fileList);
  ugtInput = new TChain ("l1uGTTree/L1uGTTree");
  appendFromFileList(ugtInput, fileList);
  
  TString fOutNameVBF = "VBFtau_ratePU" ;

  fOutNameVBF += Form("_%d",njob);
  fOutNameVBF += ".root";



  lInput->SetMakeClass(1);
  cInput->SetMakeClass(1);
  ugtInput->SetMakeClass(1);  
      
  Int_t lumi ;
  UInt_t run ;
  UShort_t stage2_tauN ;
  std::vector<float> stage2_tauEt;
  std::vector<float> stage2_tauEta;
  std::vector<float> stage2_tauPhi;
  std::vector<short> stage2_tauIso;
  std::vector<short> stage2_tauBx;  
      
      
  UShort_t stage2_jetN ;
  std::vector<float> stage2_jetEt;
  std::vector<float> stage2_jetEta;
  std::vector<float> stage2_jetPhi;
  std::vector<short> stage2_jetBx;
  std::vector<bool>  trigger_algo;    
      
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_run ;

  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;
  TBranch *b_stage2_tauBx;

  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetBx;
  TBranch *b_stage2_jetPhi;
  TBranch *b_m_algoDecisionInitial;
  

  

  lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
  lInput ->SetBranchAddress("run", &run , &b_run );
  cInput ->SetBranchAddress("nTaus", &stage2_tauN , &b_stage2_tauN );
  cInput ->SetBranchAddress("tauEta", &stage2_tauEta, &b_stage2_tauEta);
  cInput ->SetBranchAddress("tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
  cInput ->SetBranchAddress("tauEt", &stage2_tauEt, &b_stage2_tauEt);
  cInput ->SetBranchAddress("tauIso", &stage2_tauIso, &b_stage2_tauIso);
  cInput ->SetBranchAddress("tauBx", &stage2_tauBx, &b_stage2_tauBx);
  cInput ->SetBranchAddress("nJets", &stage2_jetN , &b_stage2_jetN);
  cInput ->SetBranchAddress("jetEta", &stage2_jetEta, &b_stage2_jetEta);
  cInput ->SetBranchAddress("jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
  cInput ->SetBranchAddress("jetEt", &stage2_jetEt, &b_stage2_jetEt);
  cInput ->SetBranchAddress("jetBx", &stage2_jetBx, &b_stage2_jetBx);
  ugtInput ->SetBranchAddress("m_algoDecisionInitial", &trigger_algo, &b_m_algoDecisionInitial);

  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s/%s",directory.Data(),fOutNameVBF.Data()), "recreate");

  fOutVBF->cd(); 
  //VBF

  int nbins = 60;
  int low = 30;
  int high = 60;

  int nbinsLS = 1850;
  int lowLS = 50;
  int highLS = 1900;

  //PU
  TH1D* nEventsPass     = new TH1D ("nEventsPass", "Events per PU", nbins, low, high);
  TH1D* VBF_Pass        = new TH1D ("VBF_Pass", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  TH1D* VBF_Ratio       = new TH1D ("VBF_Ratio", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  TH1D* VBF_Pass_rej    = new TH1D ("VBF_Pass_rej", "rejVBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  TH1D* VBF_Ratio_rej   = new TH1D ("VBF_Ratio_rej", "rejVBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  TH1D* VBF_Pass_rej60    = new TH1D ("VBF_Pass_rej60", "rej60VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  TH1D* VBF_Ratio_rej60   = new TH1D ("VBF_Ratio_rej60", "rej60VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  //LS
  TH1D* nEventsPassLS   = new TH1D ("nEventsPassLS", "Events per LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_PassLS      = new TH1D ("VBF_PassLS", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_RatioLS     = new TH1D ("VBF_RatioLS", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Pass_rejLS  = new TH1D ("VBF_Pass_rejLS", "rejVBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Ratio_rejLS = new TH1D ("VBF_Ratio_rejLs", "rejVBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Pass_rej60LS  = new TH1D ("VBF_Pass_rej60LS", "rej60VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Ratio_rej60LS = new TH1D ("VBF_Ratio_rej60Ls", "rej60VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);

  //PU (trigger)
  TH1D* VBF_PassTrg        = new TH1D ("VBF_PassTrg", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str, nbins, low, high);
  //LS (trigger)
  TH1D* VBF_PassLSTrg      = new TH1D ("VBF_PassLSTrg", "VBF_"+Mjj_str+"_"+jets_str+"_"+tau_str+"_LS", nbinsLS, lowLS, highLS);


  
  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);

      
      regex reg(Form("%d,%d,",Fill,Run));
      temp = regex_replace(str, reg, "");
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
  long int nEvents = 0;
 

  std::vector<object> jet;
  std::vector<object> tau;   
  std::vector<object> jetrej;
  std::vector<object> taurej;
  std::vector<object> jetrej60;
  std::vector<object> taurej60;   
 


  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej60; //VBF

   

  

  for (Long64_t iEv =0 ;true; ++iEv){
   
    lumi = 0;			   
    run = 0;
    stage2_tauN=0;		   
    stage2_tauEt.clear(); 
    stage2_tauEta.clear();
    stage2_tauPhi.clear();
    stage2_tauIso.clear();
    stage2_tauBx.clear();    
		   
    stage2_jetN =0;		   
    stage2_jetEt.clear(); 
    stage2_jetEta.clear();
    stage2_jetBx.clear();
    stage2_jetPhi.clear();

    trigger_algo.clear();

    if(iEv == maxevents) break;
    int got = 0;

    got = lInput->GetEntry(iEv);

    if (got == 0) break;
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;

    cInput->GetEntry(iEv);
    ugtInput->GetEntry(iEv);
    jet.clear();
    tau.clear();
    taurej.clear();
    jetrej.clear();;
    taurej60.clear();
    jetrej60.clear();;
    mjj_pass.clear();
    mjj_rej.clear();
    mjj_rej60.clear();
    


     if(run!=Run) continue;
    
    
    if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;


 
					      
    


    nEventsPass->Fill(PU_per_LS[lumi]);
    nEventsPassLS->Fill(lumi);




   
    
    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      for (long int iL1 = 0; iL1 < stage2_tauEt.size(); iL1++){ //loop on taus
	// selections
	double tauPt  = stage2_tauEt.at(iL1);
	double tauEta  = fabs(stage2_tauEta.at(iL1));
	double tauIso  = fabs(stage2_tauIso.at(iL1));
	int tauBx  = stage2_tauBx.at(iL1);
	
	if(tauPt>tau_cut && tauBx == 0 && tauIso>0.) {
	  tau.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;

	}
	
      }
      std::sort (tau.begin(),tau.end());
      std::sort (taurej.begin(),taurej.end());
      
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));
      int jetBx  = stage2_jetBx.at(iL1);
      
      if(jetPt>jets_cut && jetBx == 0 && tau.size()>0) {
       
	  TLorentzVector tlv_jet;
	  tlv_jet.SetPtEtaPhiM(
			       stage2_jetEt.at(iL1),
			       stage2_jetEta.at(iL1),
			       stage2_jetPhi.at(iL1),
			       0.
			       );

	  if (tau[0].DeltaRTLV(tlv_jet) < 0.2) continue;
	  jet.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
      }

    }
    std::sort (jet.begin(),jet.end());

    //removing TT28
    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      for (long int iL1 = 0; iL1 < stage2_tauEt.size(); iL1++){ //loop on taus
	// selections
	double tauPt  = stage2_tauEt.at(iL1);
	double tauEta  = fabs(stage2_tauEta.at(iL1));
	double tauIso  = fabs(stage2_tauIso.at(iL1));
	int tauBx  = stage2_tauBx.at(iL1);
	
	if(tauPt>tau_cut && tauBx == 0 && tauIso>0.) {
	  if(tauEta>2.7 && tauEta<3.0){
	    if(tauPt>60.) taurej60.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;
	}else{
	  taurej.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;
	  taurej60.push_back(object(stage2_tauEt.at(iL1),stage2_tauEta.at(iL1),stage2_tauPhi.at(iL1),stage2_tauIso.at(iL1))) ;
	  }
	}
      }
      
      
      std::sort (taurej.begin(),taurej.end());
      std::sort (taurej60.begin(),taurej60.end());
      
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));
      int jetBx  = stage2_jetBx.at(iL1);
      
      if(jetPt>jets_cut && jetBx == 0 && taurej.size()>0) {
	
	TLorentzVector tlv_jet;
	tlv_jet.SetPtEtaPhiM(
			     stage2_jetEt.at(iL1),
			     stage2_jetEta.at(iL1),
			     stage2_jetPhi.at(iL1),
			     0.
			     );
	
	if (taurej[0].DeltaRTLV(tlv_jet) > 0.2){
	  if (!(jetEta>2.7 && jetEta<3.0)){
	    jetrej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	  }
	}
      }
       if(jetPt>jets_cut && jetBx == 0 && taurej60.size()>0) {
	
	TLorentzVector tlv_jet;
	tlv_jet.SetPtEtaPhiM(
			     stage2_jetEt.at(iL1),
			     stage2_jetEta.at(iL1),
			     stage2_jetPhi.at(iL1),
			     0.
			     );
	
	if (taurej60[0].DeltaRTLV(tlv_jet) > 0.2){
	  if(jetEta>2.7 && jetEta<3.0){
	    if (jetPt>60) jetrej60.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	  }else{
	    jetrej60.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	  }
	}
       }
    }
    std::sort (jetrej.begin(),jetrej.end());
    std::sort (jetrej60.begin(),jetrej60.end());
    
    
   

    //trigger decision

    if(trigger_ID!=-999){
      if (trigger_algo[trigger_ID]){
	VBF_PassTrg->Fill(PU_per_LS[lumi]);
	VBF_PassLSTrg->Fill(lumi);
      }else{
	VBF_PassTrg->Fill(-1);
	VBF_PassLSTrg->Fill(-1);
      }
    }
    
    //VBF
       
    if (jet.size() >= 2){

      for (int iJet = 0; iJet <jet.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet.size(); kJet++){      
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet[iJet].Et(),
			      jet[iJet].Eta(),
			      jet[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet[kJet].Et(),
			      jet[kJet].Eta(),
			      jet[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_pass.push_back(make_tuple(jetPair.M(),iJet,kJet));


	  }

	}
      }
      std::sort(mjj_pass.begin(),mjj_pass.end());
      if (std::get<0>(*(mjj_pass.rbegin()))>Mjj_cut){

	
	VBF_Pass->Fill(PU_per_LS[lumi]);
	VBF_PassLS->Fill(lumi);
	
      }else{
	
	VBF_Pass->Fill(-1);
	VBF_PassLS->Fill(-1);
      }
      
    } else {

      VBF_Pass->Fill(-1);
      VBF_PassLS->Fill(-1);
    }


    //VBF rejecting TT28
       
    if (jetrej.size() >= 2){

      for (int iJet = 0; iJet <jetrej.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jetrej.size(); kJet++){      
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetrej[iJet].Et(),
			      jetrej[iJet].Eta(),
			      jetrej[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetrej[kJet].Et(),
			      jetrej[kJet].Eta(),
			      jetrej[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_rej.push_back(make_tuple(jetPair.M(),iJet,kJet));


	  }

	}
      }
      std::sort(mjj_rej.begin(),mjj_rej.end());
      if (std::get<0>(*(mjj_rej.rbegin()))>Mjj_cut){

	
	VBF_Pass_rej->Fill(PU_per_LS[lumi]);
	VBF_Pass_rejLS->Fill(lumi);
	
      }else{
	
	VBF_Pass_rej->Fill(-1);
	VBF_Pass_rejLS->Fill(-1);
      }
      
    } else {

      VBF_Pass_rej->Fill(-1);
      VBF_Pass_rejLS->Fill(-1);
    }

    // VBF rejecting only TT28, pt<60 GeV
        if (jetrej60.size() >= 2){

      for (int iJet = 0; iJet <jetrej60.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jetrej60.size(); kJet++){      
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jetrej60[iJet].Et(),
			      jetrej60[iJet].Eta(),
			      jetrej60[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jetrej60[kJet].Et(),
			      jetrej60[kJet].Eta(),
			      jetrej60[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_rej60.push_back(make_tuple(jetPair.M(),iJet,kJet));


	  }

	}
      }
      std::sort(mjj_rej60.begin(),mjj_rej60.end());
      if (std::get<0>(*(mjj_rej60.rbegin()))>Mjj_cut){

	
	VBF_Pass_rej60->Fill(PU_per_LS[lumi]);
	VBF_Pass_rej60LS->Fill(lumi);
	
      }else{
	
	VBF_Pass_rej60->Fill(-1);
	VBF_Pass_rej60LS->Fill(-1);
      }
      
	} else {

      VBF_Pass_rej60->Fill(-1);
      VBF_Pass_rej60LS->Fill(-1);
    }

  }
  
  

  cout << "Making histograms..." << endl; 
  
  //VBF
  for (int i = 0; i <VBF_Pass->GetNbinsX(); i++){
    double bin = 0;
    if(nEventsPass->GetBinContent(i)>0)   bin = 1.*(VBF_Pass->GetBinContent(i))/nEventsPass->GetBinContent(i);
    VBF_Ratio -> SetBinContent (i, bin);        
    //rej
    bin = 1.*(VBF_Pass_rej->GetBinContent(i))/nEventsPass->GetBinContent(i);
    VBF_Ratio_rej -> SetBinContent (i, bin);
    //rej60
    bin = 1.*(VBF_Pass_rej60->GetBinContent(i))/nEventsPass->GetBinContent(i);
    VBF_Ratio_rej60 -> SetBinContent (i, bin);        
  }

  cout<< "Output saved in "<<directory<<"/"<<fOutNameVBF<<endl;

fOutVBF->cd();
fOutVBF -> Write();
}


