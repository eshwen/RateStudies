// evaluates rate  - compile with c++ -lm -o EvalRate_OnL1Legacy EvalRate_OnL1Legacy.cpp `root-config --glibs --cflags`
// and launch the executable 

#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

using namespace std;

void EvalRate(){
  //double scaleRunI = 1./1.713;
  double scaleRunI = 1.075/1.713;
  bool doScaleRunI = true;
  
  cout << "Do scale: " << doScaleRunI << endl;
  cout << "Scale factor: " << scaleRunI << endl;
  
  //    ZeroBias sample L1
  TFile *file = TFile::Open("test.root","read");
  TTree * tInput = (TTree*) file->Get("L1Tree/L1Tree");
  
    TString fOutName = "rateL1_test_";
    
    if (doScaleRunI) fOutName += "Scaled.root";
    else fOutName += "NotScaled.root";
    cout<<"check"<<endl;    


    //  Int_t stage2objNumber; //0 = taus, 1 = eg, 2 = jets, 3 = muons  

  
  Int_t stage2_tauN;
  std::vector<Float_t> stage2_tauEt;
  std::vector<Float_t> stage2_tauEta;
  std::vector<Float_t> stage2_tauPhi;

  
  Int_t stage2_jetN;
  std::vector<Float_t> stage2_jetEt;
  std::vector<Float_t> stage2_jetEta;
  std::vector<Float_t> stage2_jetPhi;

  // set branch and variables
      cout<<"declared"<<endl;    
  
  tInput ->SetBranchAddress ("stage2_tauN", &stage2_tauN);
  tInput ->SetBranchAddress ("stage2_tauEta", &stage2_tauEta);
  tInput ->SetBranchAddress ("stage2_tauPhi", &stage2_tauPhi);
  tInput ->SetBranchAddress ("stage2_tauPhi", &stage2_tauPhi);
  tInput ->SetBranchAddress ("stage2_jetN", &stage2_jetN);
  tInput ->SetBranchAddress ("stage2_jetEta", &stage2_jetEta);
  tInput ->SetBranchAddress ("stage2_jetPhi", &stage2_jetPhi);
  tInput ->SetBranchAddress ("stage2_jetPhi", &stage2_jetPhi);
       cout<<"addressed"<<endl;    

       tInput->Show();  

    ///////////////
    //// HISTO ////
    ///////////////

    TFile* fOut = new TFile (fOutName, "recreate");
   
    // the tree histograms are filed in identical way, structure is simply copied from the Stage 2 rate Evaluator
    TH1D* LeadPt_passIso = new TH1D ("LeadPt_passIso", "LeadPt_passIso", 5000, 0, 5000);
    TH1D* SubleadPt_passIso = new TH1D ("SubleadPt_passIso", "SubleadPt_passIso", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau = new TH1D ("RelativeRate_singleTau", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau = new TH1D ("RelativeRate_diTau", "Relative rate - double tau", 5000, 0, 5000);

    TH1D* LeadPt_NoIso = new TH1D ("LeadPt_NoIso", "LeadPt_NoIso", 5000, 0, 5000);
    TH1D* SubleadPt_NoIso = new TH1D ("SubleadPt_NoIso", "SubleadPt_NoIso", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau_NoIso = new TH1D ("RelativeRate_singleTau_NoIso", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau_NoIso = new TH1D ("RelativeRate_diTau_NoIso", "Relative rate - double tau", 5000, 0, 5000);

    TH1D* LeadPt_IsoAndShape = new TH1D ("LeadPt_IsoAndShape", "LeadPt_IsoAndShape", 5000, 0, 5000);
    TH1D* SubleadPt_IsoAndShape = new TH1D ("SubleadPt_IsoAndShape", "SubleadPt_IsoAndShape", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau_IsoAndShape = new TH1D ("RelativeRate_singleTau_IsoAndShape", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau_IsoAndShape = new TH1D ("RelativeRate_diTau_IsoAndShape", "Relative rate - double tau", 5000, 0, 5000);

        cout<<"histodeclared"<<endl;    
    // analyze data    
    long int nEvents = tInput->GetEntries(); 
    //long int nEvents = 100000; 
    
    std::vector<double> pt_pass; // just pass iso
    std::vector<double> pt_passAndShape; // also pass veto shape
    std::vector<double> pt_noiso;
      cout<<"vectorsdeclared"<<endl;    
    // loop on all events
    for (long int iEv = 0; iEv < nEvents; iEv++)
    {
        tInput->GetEntry(iEv);
  cout<<"loop get entry"<<endl;         
   if (iEv%100000 == 0) cout << iEv << " / " << nEvents << endl;
        
        // clear pt vector for this event and compute number of TT
        pt_pass.clear();
        pt_noiso.clear();
        pt_passAndShape.clear();
        
        
        // loop on all L1 taus --> save all taus passing Iso requirement + other selections
        //cout << iEv << " " << __stage2_tauN << " " << stage2_tauEta->size() << endl;
        for (int iL1 = 0; iL1 < stage2_tauN; iL1++)
        {
            // selections
            double abseta = TMath::Abs( (stage2_tauEta)[iL1] );
            double tauPt  = (stage2_tauPhi)[iL1];
            
            if (doScaleRunI) tauPt = scaleRunI*tauPt;

            if (abseta < 2.2)
            {
                pt_noiso.push_back (tauPt);
            
                // apply iso
                bool PassIso = true; // iso is applied by default in the collection
                if (PassIso)
                {
                    pt_pass.push_back (tauPt);
                    
                    
                    if (true) // there is no shape to test in Legacy trigger
					{
                        pt_passAndShape.push_back (tauPt);
                    }
                }
            }
        }

    
        // now that all taus in the event are analyzed, fill the histograms of lead and sublead pt (if missing, fil with a -1)
        std::sort (pt_pass.begin(), pt_pass.end());
        std::sort (pt_noiso.begin(), pt_noiso.end());
        std::sort (pt_passAndShape.begin(), pt_passAndShape.end());
        
        // no iso operations
        if (pt_noiso.size() >= 2 )
        {
            LeadPt_NoIso -> Fill ( *(pt_noiso.rbegin()) ); 
            SubleadPt_NoIso -> Fill ( *(pt_noiso.rbegin() + 1) );
        } 
        
        else if (pt_noiso.size() == 1) 
        {
            LeadPt_NoIso -> Fill ( *(pt_noiso.rbegin()) ); 
            SubleadPt_NoIso -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_NoIso -> Fill ( -1 ); 
            SubleadPt_NoIso -> Fill ( -1 );        
        }       


        
        // pass iso operations
        if (pt_pass.size() >= 2 )
        {
            LeadPt_passIso -> Fill ( *(pt_pass.rbegin()) ); 
            SubleadPt_passIso -> Fill ( *(pt_pass.rbegin() + 1) );
        } 
        
        else if (pt_pass.size() == 1) 
        {
            LeadPt_passIso -> Fill ( *(pt_pass.rbegin()) ); 
            SubleadPt_passIso -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_passIso -> Fill ( -1 ); 
            SubleadPt_passIso -> Fill ( -1 );        
        }      
        
        
        // also shape operations
        if (pt_passAndShape.size() >= 2 )
        {
            LeadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin()) ); 
            SubleadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin() + 1) );
        } 
        
        else if (pt_passAndShape.size() == 1) 
        {
            LeadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin()) ); 
            SubleadPt_IsoAndShape -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_IsoAndShape -> Fill ( -1 ); 
            SubleadPt_IsoAndShape -> Fill ( -1 );        
        }     
    }
    
    // compute rate plots
   cout << "Computing rates..." << endl; 
    
    for (int i = 1; i <= 5000; i++)
    {
        // no Iso
        int Tot_singleTau_noiso = LeadPt_NoIso -> Integral (0, 5001);
        int Tot_diTau_noiso = SubleadPt_NoIso -> Integral (0, 5001);
        
        double relRateSingle_noiso = 1.*(LeadPt_NoIso->Integral(i, 5001))/Tot_singleTau_noiso;
        double relRateDouble_noiso = 1.*(SubleadPt_NoIso->Integral(i, 5001))/Tot_diTau_noiso;
        
        RelativeRate_singleTau_NoIso -> SetBinContent (i, relRateSingle_noiso);
        RelativeRate_diTau_NoIso -> SetBinContent (i, relRateDouble_noiso);

        // with Iso
        int Tot_singleTau = LeadPt_passIso -> Integral (0, 5001);
        int Tot_diTau = SubleadPt_passIso -> Integral (0, 5001);
        
        double relRateSingle = 1.*(LeadPt_passIso->Integral(i, 5001))/Tot_singleTau;
        double relRateDouble = 1.*(SubleadPt_passIso->Integral(i, 5001))/Tot_diTau;
        
        RelativeRate_singleTau -> SetBinContent (i, relRateSingle);
        RelativeRate_diTau -> SetBinContent (i, relRateDouble);        

        // with Iso and shape
        int Tot_singleTau_shape = LeadPt_IsoAndShape -> Integral (0, 5001);
        int Tot_diTau_shape = SubleadPt_IsoAndShape -> Integral (0, 5001);
        
        double relRateSingle_shape = 1.*(LeadPt_IsoAndShape->Integral(i, 5001))/Tot_singleTau_shape;
        double relRateDouble_shape = 1.*(SubleadPt_IsoAndShape->Integral(i, 5001))/Tot_diTau_shape;
        
        RelativeRate_singleTau_IsoAndShape -> SetBinContent (i, relRateSingle_shape);
        RelativeRate_diTau_IsoAndShape -> SetBinContent (i, relRateDouble_shape);        

    }
    
    fOut -> Write();
}
