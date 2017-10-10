// evaluates rate  - compile with c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
// and launch the executable 
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
  int leadjet_cut = atoi(argv[2]);
  int subleadjet_cut = atoi(argv[3]);
  std::stringstream Mjj_ss;
  Mjj_ss << Mjj_cut;
  TString Mjj_str = Mjj_ss.str();
  std::stringstream leadjet_ss;
  leadjet_ss << leadjet_cut;
  string leadjet_str = leadjet_ss.str();
  std::stringstream subleadjet_ss;
  subleadjet_ss << subleadjet_cut;
  string subleadjet_str = subleadjet_ss.str();
  

  TString year = argv[4];
  cout<<year<<" year"<<endl;  

  //ZeroBias sample L1
  //TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/ZeroBiasRun2017A/";
  //2016
  //TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias26Apr2017HighPU_ibx0_BunchTrain0-5_2016H9Nov_-1Events/";
  //2017
  //

  //  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6245/";
  TString directory = argv[5];
  //TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6194/";
  //TString directory = "/eos/user/c/camendol/local/";
  //TFile *file = TFile::Open(Form("%sL1total.root", directory.Data()),"read");
  //2016
  //TString fileList = "fileLists/L1NtuplesL1Menu2017ZeroBiasRun2016H_BunchTrainsX.txt";
  //2017
  //  TString fileList = "fileLists/L1NtuplesL1Menu2017ZeroBiasFill6245.txt";
   //  TString fileList = "fileLists/ZeroBias2017_HighPU.list";

  //  TString fileList= "fileLists/ZeroBias2017_fill6245.list";	
  TString fileList  = argv[6];
  cout<<fileList<<" filelist"<<endl;
  bool emulated  =false;
  bool reweight = false;
    
  //  ifstream PUFile("utils/PU_per_LS.txt"); //2016
  //ifstream PUFile("utils/PU_per_LS_fill5824_2017.txt"); //2017
  //ifstream PUFile("utils/PU_per_LS_fill6194_2017.txt"); //2017

  //ifstream PUFile("utils/PU_per_LS_fill6245_2017.txt"); //2017
  TString PUperLumiFile = argv[7];
  ifstream PUFile(PUperLumiFile); //2017
      

  TChain * cInput;
  TChain * lInput; 
      
      
      
      
  if (emulated){
    cInput = new TChain ("l1UpgradeEmuTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);
  }else{
    cInput= new TChain ("l1UpgradeTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);
	
  }
      
  TString fOutNameVBF;
      
      
  if(emulated){
    fOutNameVBF = directory+"Emu_rateL1_VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_"+year+"_PU_LS.root";
	
  }else{
    fOutNameVBF = directory+"rateL1_VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_"+year+"_PU_testArgv.root";
	
  }
      
      
  lInput->SetMakeClass(1);
  cInput->SetMakeClass(1);  
      
  Int_t lumi ;
  UInt_t run ;
  UShort_t stage2_tauN ;
  std::vector<float> stage2_tauEt;
  std::vector<float> stage2_tauEta;
  std::vector<float> stage2_tauPhi;
  std::vector<short> stage2_tauIso;  
      
      
  UShort_t stage2_jetN ;
  std::vector<float> stage2_jetEt;
  std::vector<float> stage2_jetEta;
  std::vector<float> stage2_jetPhi;
      
      
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_run ;

  TBranch *b_stage2_tauN ;
  TBranch *b_stage2_tauEt;
  TBranch *b_stage2_tauEta;
  TBranch *b_stage2_tauPhi;
  TBranch *b_stage2_tauIso;

  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  
  

  

  lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
  lInput ->SetBranchAddress("run", &run , &b_run );
  cInput ->SetBranchAddress("nTaus", &stage2_tauN , &b_stage2_tauN );
  cInput ->SetBranchAddress("tauEta", &stage2_tauEta, &b_stage2_tauEta);
  cInput ->SetBranchAddress("tauPhi", &stage2_tauPhi, &b_stage2_tauPhi);
  cInput ->SetBranchAddress("tauEt", &stage2_tauEt, &b_stage2_tauEt);
  cInput ->SetBranchAddress("tauIso", &stage2_tauIso, &b_stage2_tauIso);
  cInput ->SetBranchAddress("nJets", &stage2_jetN , &b_stage2_jetN);
  cInput ->SetBranchAddress("jetEta", &stage2_jetEta, &b_stage2_jetEta);
  cInput ->SetBranchAddress("jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
  cInput ->SetBranchAddress("jetEt", &stage2_jetEt, &b_stage2_jetEt);

  
  ///////////////
  //// HISTO ////
  ///////////////
  
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");

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
  TH1D* VBF_Pass        = new TH1D ("VBF_Pass", "VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str, nbins, low, high);
  TH1D* VBF_Ratio       = new TH1D ("VBF_Ratio", "VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str, nbins, low, high);
  TH1D* VBF_Pass_rej    = new TH1D ("VBF_Pass_rej", "rejVBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str, nbins, low, high);
  TH1D* VBF_Ratio_rej   = new TH1D ("VBF_Ratio_rej", "rejVBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str, nbins, low, high);
  //LS
  TH1D* nEventsPassLS   = new TH1D ("nEventsPassLS", "Events per LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_PassLS      = new TH1D ("VBF_PassLS", "VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_RatioLS     = new TH1D ("VBF_RatioLS", "VBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Pass_rejLS  = new TH1D ("VBF_Pass_rejLS", "rejVBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_LS", nbinsLS, lowLS, highLS);
  TH1D* VBF_Ratio_rejLS = new TH1D ("VBF_Ratio_rejLs", "rejVBF_"+Mjj_str+"_"+leadjet_str+"_"+subleadjet_str+"_LS", nbinsLS, lowLS, highLS);


  
  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);

      //      temp.ReplaceAll("5412,283171,",""); //2016
      
      regex reg("6245,303948,");
      temp = regex_replace(str, reg, "");
      int pos_coma = temp.First(",");
      TString LS_str(temp,pos_coma);
      //cout<<LS_str<<endl;
      TString Replacing = LS_str ;
      Replacing += ",";
      temp.ReplaceAll(Replacing.Data(),"");
      TString PU_str = temp;
      //      cout<<PU_str<<endl;
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




 
 
  std::vector<object> jet30;   
  std::vector<object> jet30rej;   
 


  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_pass_sortPt; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_rej_sortPt; //VBF



   
  int Nfill = 0;  
  

  for (Long64_t iEv =0 ;true; ++iEv){
    lumi = 0;			   
    run = 0;
    stage2_tauN=0;		   
    stage2_tauEt.clear(); 
    stage2_tauEta.clear();
    stage2_tauPhi.clear();
    stage2_tauIso.clear();    
		   
		   
    stage2_jetN =0;		   
    stage2_jetEt.clear(); 
    stage2_jetEta.clear();
    stage2_jetPhi.clear();
    int got = 0;

    //  cInput->SetBranchStatus("*",0);
    
   
    got = lInput->GetEntry(iEv);
    
    if (got == 0) break;

   
  

    
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
   
    //    cInput->SetBranchStatus("*",1);
    cInput->GetEntry(iEv);
    jet30.clear();
    jet30rej.clear();;
    mjj_pass.clear();
    mjj_pass_sortPt.clear();
    mjj_rej.clear();
    mjj_rej_sortPt.clear();
    

    //      if (lumi < 890) continue;    
 

    if (lumi ==53) continue;
    if (lumi ==54) continue;
    if (lumi ==131) continue;
    if (lumi ==132) continue;
    if (lumi ==177) continue;
    if (lumi ==283) continue;
    if (lumi ==284) continue;
    if (lumi ==319) continue;
    if (lumi ==320) continue;
    if (lumi ==389) continue;
    if (lumi ==390) continue;
    if (lumi ==457) continue;
    if (lumi ==458) continue;
    if (lumi ==463) continue;
    if (lumi ==464) continue;
    if (lumi ==465) continue;
    if (lumi >559 && lumi<569) continue;
    if (lumi >1678 && lumi<1725) continue;
    
   
    //       if(lumi<56 || lumi>69) continue; //2016

    // if(run>303948) continue;
      ///if(lumi<135 || lumi>264) continue; //2017



    
    if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;



    
    Float_t weight = 0;
    if (reweight){
      weight = PU_per_LS[56]/PU_per_LS[lumi];
    }else{
      weight =  1.;
    }

    nEventsPass->Fill(PU_per_LS[lumi]);
    nEventsPassLS->Fill(lumi);



   

    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));

      if(jetPt>30.) {
	jet30.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	if(jetEta>2.7 && jetEta<3.0){
	 if(jetPt>60.) jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}else{
	  jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}
      }
      
    }
    
 
    std::sort (jet30.begin(),jet30.end());
    std::sort (jet30rej.begin(),jet30rej.end());
   
 
    //VBF
       
    if (jet30.size() >= 2){

      for (int iJet = 0; iJet <jet30.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30.size(); kJet++){      
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


	    if(jetPair.M()>=Mjj_cut){
	      // cout<<jet30[iJet].Et()<<" "<<jet30[kJet].Et()<<" "<<jetPair.M()<<endl;
	      mjj_pass_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30[iJet].Et(),jet30[kJet].Et()));
	    }
	  }

	}
      }
      std::sort(mjj_pass.begin(),mjj_pass.end());
      std::sort(mjj_pass_sortPt.begin(),mjj_pass_sortPt.end(),SortMjjByJetThreshold);
      if(mjj_pass_sortPt.size()>0) {

	if (std::get<0>(*(mjj_pass.rbegin()))>Mjj_cut && std::get<4>(*(mjj_pass_sortPt.rbegin()))>subleadjet_cut && jet30[0].Et() > leadjet_cut) {
	  //cout<<"pass "<<std::get<0>(*(mjj_pass.rbegin())) <<" "<<  std::get<4>(*(mjj_pass_sortPt.rbegin())) << " " << jet30[0].Et()<<endl;
	  VBF_Pass->Fill(PU_per_LS[lumi],weight);
	  VBF_PassLS->Fill(lumi,weight);

	}else{
	  //cout<<"fail"<<endl;
	  VBF_Pass->Fill(-1);
	  VBF_PassLS->Fill(-1);
	}
      }else{
	//	cout<<"fail"<<endl;
	VBF_Pass->Fill(-1);
	VBF_PassLS->Fill(-1);
      }
    } else {
      //	cout<<"fail"<<endl;
      VBF_Pass->Fill(-1);
      VBF_PassLS->Fill(-1);
    }
    //VBF rejecting TT28

    if (jet30rej.size() >= 2){
      for (int iJet = 0; iJet <jet30rej.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30rej.size(); kJet++){      
	  if (kJet!=iJet) {
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet30rej[iJet].Et(),
			      jet30rej[iJet].Eta(),
			      jet30rej[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet30rej[kJet].Et(),
			      jet30rej[kJet].Eta(),
			      jet30rej[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_rej.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>=Mjj_cut) mjj_rej_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30rej[iJet].Et(),jet30rej[kJet].Et()));
	  }
	}
      }
      std::sort(mjj_rej.begin(),mjj_rej.end());
      std::sort(mjj_rej_sortPt.begin(),mjj_rej_sortPt.end(),SortMjjByJetThreshold);

      if(mjj_rej_sortPt.size()>0) {
	if (std::get<0>(*(mjj_rej.rbegin()))>Mjj_cut && std::get<4>(*(mjj_pass_sortPt.rbegin()))>subleadjet_cut && jet30rej[0].Et() > leadjet_cut){
	  VBF_Pass_rej->Fill(PU_per_LS[lumi],weight);
	   VBF_Pass_rejLS->Fill(lumi,weight);
	}else{
	VBF_Pass_rej->Fill(-1);
	VBF_Pass_rejLS->Fill(-1);
	}
      }else{
	VBF_Pass_rej->Fill(-1);
	VBF_Pass_rejLS->Fill(-1);
      }
    } else {
      VBF_Pass_rej->Fill(-1);
      VBF_Pass_rejLS->Fill(-1);
    }  
  }
  
  //cout << endl;
  cout << "Making histograms..." << endl; 
  
  //VBF
  for (int i = 0; i <VBF_Pass->GetNbinsX(); i++){
    double bin = 0;
    if(nEventsPass->GetBinContent(i)>0)   bin = 1.*(VBF_Pass->GetBinContent(i))/nEventsPass->GetBinContent(i);
    VBF_Ratio -> SetBinContent (i, bin);        
    //rej
    bin = 1.*(VBF_Pass_rej->GetBinContent(i))/nEventsPass->GetBinContent(i);
    VBF_Ratio_rej -> SetBinContent (i, bin);        
  }



fOutVBF->cd();
fOutVBF -> Write();

}

