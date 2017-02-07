
// evaluates rate  - compile with c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
// and launch the executable 

#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

using namespace std;

class object{
public:
  object(float Et, float Eta, float Phi){Et_=Et; Eta_ = Eta;Phi_ = Phi;};
  float Et() const{
    return Et_;
  }
  float Phi() const {
    return Phi_;
  }
  float Eta() const {
    return Eta_;
  }
  
  class SortByDeltaR{
  public:    
    SortByDeltaR(float ref_Phi, float ref_Eta){ref_Phi_ = ref_Phi; ref_Eta_= ref_Eta;};
    bool operator() (const object& obj1,const object& obj2)
    {
      float dEta1 = obj1.Eta()-ref_Eta_;
      
      float dPhi1 = obj1.Phi()-ref_Phi_;
      if(dPhi1 > 3.14) dPhi1 -= 6.28;
      if(dPhi1 < -3.14) dPhi1 += 6.28;
      
      float dEta2 = obj2.Eta()-ref_Eta_;
      
      float dPhi2 = obj2.Phi()-ref_Phi_;
      if(dPhi2 > 3.14) dPhi2 -= 6.28;
      if(dPhi2 < -3.14) dPhi2 += 6.28;
      return sqrt(dEta1*dEta1 + dPhi1*dPhi1) < sqrt(dEta2*dEta2 + dPhi2*dPhi2);
    }
  private: 
    float ref_Phi_, ref_Eta_;  
    
  }; 
  
  bool inline operator <(const object& objref) const {
    return this->Et() > objref.Et();
  }
  
  
  friend std::ostream& operator << (std::ostream &out, const object &object);
private:
  float Et_, Eta_, Phi_;
};


std::ostream & operator << (std::ostream &out, const object &object)
{
  out<<object.Et();
  return out;
}


int main(int argc, char** argv){

  float nbStudiedRun = 96;
  float thisLumiRun = 7.6E32;
  float scaleToLumi = 1.8E34;
  float scale = 0.001*(nbStudiedRun*11245.6)*scaleToLumi/thisLumiRun;  

  

  cout << "Scale factor: " << scale << endl;
  // TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias25Gen2016HighPU_ibx0_100000Events/";
  //  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias6Feb2017HighPU_ibx0_-1Events/";
  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesOutput_ZeroBias7Feb2017HighPU_ibx0_BunchTrain1-5_2016H9Nov_-1Events/";

  //ZeroBias sample L1
  TFile *file = TFile::Open(Form("%sL1total.root", directory.Data()),"read");
  TTree * tInput = (TTree*) file->Get("L1Tree/L1Tree");
  
  TString fOutName = "rateL1_VBF";
  fOutName = directory+fOutName;
  fOutName += ".root";

    
  
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
  
 
  TH1D* LeadTauPt_Pass = new TH1D ("LeadTauPt_Pass", "LeadTauPt_Pass", 50, 0, 100);
  TH1D* SingleTauPt_Pass = new TH1D ("SingleTauPt_Pass", "SingleTauPt_Pass", 50, 0, 100);
  TH1D* SubleadTauPt_Pass = new TH1D ("SubleadTauPt_Pass", "SubleadTauPt_Pass", 50, 0, 100);
  TH1D* DiTauPt_Pass = new TH1D ("DiTauPt_Pass", "DiTauPt_Pass", 50, 0, 100);
  TH1D* DiTauSingleJetPt_Pass = new TH1D ("DiTauSingleJetPt_Pass", "DiTauSingleJetPt_Pass", 50, 0, 100);
  TH2D* DiJet2D_Pass = new TH2D ("DiJet2D", "double jet const", 100, 0, 800, 100,30, 130);

  TH1D* RelativeRate_singleTau = new TH1D ("RelativeRate_singleTau", "Relative rate - single tau", 50, 0, 100);
  TH1D* RelativeRate_diTau = new TH1D ("RelativeRate_diTau", "Relative rate - double tau", 50, 0, 100);
  TH1D* RelativeRate_diTauSublead = new TH1D ("RelativeRate_diTauSublead", "Relative rate - double tau", 50, 0, 100);
  TH1D* RelativeRate_diTauSingleJet = new TH1D ("RelativeRate_diTauSingleJet", "Relative rate - double tau plus jet", 50, 0, 100);
  TH2D* RelativeRate_SingleTauSingleJet2D = new TH2D ("RelativeRate_SingleTauSingleJet2D", "Relative rate - double tau ", 50, 0, 100, 50, 0, 100);
 
  //VBF
  TH2D* Ratio_DiJet2D = new TH2D ("Ratio_DiJet2D", "Ratio - 2 jets, E_{T} > 30 GeV ", 100, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D = new TH2D ("Rate_DiJet2D", "Rate - 2 jets, E_{T} > 30 GeV ", 100, 0, 800, 100, 30, 130); 

  
  
  // analyze data    
  long int nEvents = tInput->GetEntries(); 
  std::vector<double> ptTau_pass; 
  std::vector<double> ptJet_pass;

  std::vector<double> mjj_pass; //VBF
  double singletau = 0;
  double ditau = 0;  
  // loop on all events
  for (long int iEv = 0; iEv < nEvents; iEv++){
    
    tInput->GetEntry(iEv);
    
     if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
     std::vector<object> jet30;   
    //std::vector<object> tau;
    //vector <float> DeltaRmin_taujet; //needed to reject tau-jet overlap
    //vector <float> DeltaRmin_jetjet; //needed to reject jet-jet overlap    
     ptTau_pass.clear();
     ptJet_pass.clear();
     mjj_pass.clear();
     
     for (long int iL1 = 0; iL1 < stage2_tauN; iL1++){
       // selections
       double tauPt  = stage2_tauEt->at(iL1);
       ptTau_pass.push_back (tauPt);
     }
    
     for (long int iL1 = 0; iL1 < stage2_jetN; iL1++){
       // selections
       double jetPt  = (*stage2_jetEt)[iL1];
       ptJet_pass.push_back (jetPt);
       if(jetPt>30.)	jet30.push_back(object(stage2_jetEt->at(iL1),stage2_jetEta->at(iL1),stage2_jetPhi->at(iL1))) ;//VBF 
     }

     std::sort (ptTau_pass.begin(), ptTau_pass.end());
     std::sort (ptJet_pass.begin(), ptJet_pass.end());
    
     //VBF    
     std::sort (jet30.begin(),jet30.end());
     int iJet = 0;
     int kJet = 0;
    if (jet30.size() >= 2){
      for (iJet= 0; iJet <jet30.size(); iJet++){      
	for (kJet= 0; kJet <jet30.size(); kJet++){      
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
	    mjj_pass.push_back(jetPair.M());

	  }

	}

      }
      std::sort(mjj_pass.begin(),mjj_pass.end());
      DiJet2D_Pass -> Fill ( *(mjj_pass.rbegin()),jet30[0].Et());
    } else {
      DiJet2D_Pass -> Fill (-1, -1);
    }

    
    
    if (ptTau_pass.size() >= 2 )
      {
	ditau += 1;
	singletau += 1;
	LeadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
	SingleTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
	SubleadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin() + 1) );
	DiTauPt_Pass -> Fill (*(ptTau_pass.rbegin()) +*(ptTau_pass.rbegin() + 1) );
	if (ptJet_pass.size()>=1){
	  DiTauSingleJetPt_Pass->Fill( *(ptJet_pass.rbegin()) );
	}else{
	  DiTauSingleJetPt_Pass->Fill( -1 );
	}
      } 
    
    else if (ptTau_pass.size() == 1){
      singletau +=1;
      LeadTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
      SingleTauPt_Pass -> Fill ( *(ptTau_pass.rbegin()) ); 
      SubleadTauPt_Pass -> Fill ( -1 );
      DiTauPt_Pass -> Fill ( -1); 
      DiTauSingleJetPt_Pass->Fill( -1 );
    } 
    
    else
      {
	LeadTauPt_Pass -> Fill ( -1 ); 
	SingleTauPt_Pass -> Fill ( -1 ); 
	SubleadTauPt_Pass -> Fill ( -1 );        
	DiTauPt_Pass -> Fill ( -1 );        
	DiTauSingleJetPt_Pass->Fill( -1 );
      }       
  }
  
  // compute rate plots
  cout << "Computing rates..." << endl; 
  


  //taus
  for (int i = 1; i <= 50; i++){

    int Tot_diTauSublead = SubleadTauPt_Pass -> Integral (1,51);
    
    int Tot_singleTau = SingleTauPt_Pass -> Integral (1,51);
    int Tot_diTau = DiTauPt_Pass -> Integral (1,51);        
    int Tot_diTauSingleJet = DiTauSingleJetPt_Pass -> Integral (1,51);        

    double relRateSingleTau = 1.*(SingleTauPt_Pass->Integral(i, 51))/Tot_singleTau;
    double relRateDiTau = 1.*(DiTauPt_Pass->Integral(i, 51))/Tot_diTau;
    double relRateDiTauSingleJet = 1.*(DiTauSingleJetPt_Pass->Integral(i, 51))/Tot_diTauSingleJet;
    double relRateDiTauSublead =  1.*(SubleadTauPt_Pass->Integral(i, 51))/Tot_diTauSublead;

    
    RelativeRate_singleTau -> SetBinContent (i, relRateSingleTau);
    RelativeRate_diTau -> SetBinContent (i, relRateDiTau);        
    RelativeRate_diTauSublead -> SetBinContent (i, relRateDiTauSublead);        
    
    RelativeRate_diTauSingleJet -> SetBinContent (i, relRateDiTauSingleJet);        

  }

  cout<<"Fraction of events with at least 1 tau: "<<singletau/nEvents<<endl;
  cout<<"Fraction of events with at least 2 taus: "<<ditau/nEvents<<endl;   


  
  //VBF
  for (int i = 1; i <=DiJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Pass->GetNbinsY();j++ ){
      double binDiJet2D = 1.*(DiJet2D_Pass->Integral(i, DiJet2D_Pass->GetNbinsX()+1,j,DiJet2D_Pass->GetNbinsY()+1))/nEvents;
      Ratio_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      
    }
  }


  fOut -> Write();
}
