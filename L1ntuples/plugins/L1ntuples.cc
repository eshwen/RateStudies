// -*- C++ -*-
//
// Package:    RateStudies/L1ntuples
// Class:      L1ntuples
// 
/**\class L1ntuples L1ntuples.cc RateStudies/L1ntuples/plugins/L1ntuples.cc

 Description: Makes ntuples with L1 trigger information out of miniAOD files

*/
//
// Original Author:  Chiara Amendola
//         Created:  Thu, 19 Jan 2017 10:19:07 GMT
//
//


// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"


//from htautauntuplizer
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>

//#include <DataFormats/PatCandidates/interface/Muon.h>
//#include <DataFormats/PatCandidates/interface/MET.h>
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>

//#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
//#include <DataFormats/PatCandidates/interface/Electron.h>
//#include <DataFormats/METReco/interface/PFMET.h>
//#include <DataFormats/METReco/interface/PFMETCollection.h>
//#include <DataFormats/JetReco/interface/PFJet.h>
//#include <DataFormats/JetReco/interface/PFJetCollection.h>
//#include <DataFormats/PatCandidates/interface/PackedCandidate.h>




#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h>


#include "Geometry/Records/interface/HcalParametersRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/GeometryObjects/interface/HcalParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

//
// class declaration
//



class L1ntuples : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1ntuples(const edm::ParameterSet&);
      ~L1ntuples();
   private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void Initialize();
  int* FillStage2(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event);
  // ----------member data ---------------------------
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  
  TTree *myTree;
  
  edm::EDGetTokenT<BXVector<l1t::Tau> > theStage2TauTag;
  edm::EDGetTokenT<BXVector<l1t::Jet> > theStage2JetTag;
  
  Int_t *stage2objNumber; //0 = taus, 1 = eg, 2 = jets, 3 = muons  
  ULong64_t _indexevents;
  Int_t _runNumber;
  Int_t _lumi;
  
  Int_t _stage2_tauN;
  std::vector<Float_t> _stage2_tauEt;
  std::vector<Float_t> _stage2_tauEta;
  std::vector<Float_t> _stage2_tauPhi;
    std::vector<int> _stage2_tauIso;

  
  Int_t _stage2_jetN;
  std::vector<Float_t> _stage2_jetEt;
  std::vector<Float_t> _stage2_jetEta;
  std::vector<Float_t> _stage2_jetPhi;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

L1ntuples::L1ntuples(const edm::ParameterSet& pset) : 
  theStage2TauTag      (consumes<BXVector<l1t::Tau>>                     (pset.getParameter<edm::InputTag>("stage2TauCollection"))),
  theStage2JetTag      (consumes<BXVector<l1t::Jet>>                     (pset.getParameter<edm::InputTag>("stage2JetCollection")))
{
   //now do what ever initialization is needed

  Initialize();
}


L1ntuples::~L1ntuples()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1ntuples::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  Initialize();
   using namespace edm;

  edm::Handle<BXVector<l1t::Tau>>stage2TauHandle;
  edm::Handle<BXVector<l1t::Jet>>stage2JetHandle; 

  event.getByToken(theStage2TauTag,stage2TauHandle);
  event.getByToken(theStage2JetTag,stage2JetHandle);

  const BXVector<l1t::Tau>* stage2Tau = stage2TauHandle.product();
  const BXVector<l1t::Jet>* stage2Jet = stage2JetHandle.product(); 
  
  _indexevents = event.id().event();
  _runNumber = event.id().run();
  _lumi=event.luminosityBlock();

  stage2objNumber=FillStage2(stage2Tau, stage2Jet, event);
  _stage2_tauN = stage2objNumber[0];
  _stage2_jetN = stage2objNumber[2];
  myTree->Fill();  
   

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1ntuples::beginJob()
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("L1Tree","L1Tree");
  myTree->Branch("EventNumber",&_indexevents,"EventNumber/l");
  myTree->Branch("RunNumber",&_runNumber,"RunNumber/I");
  myTree->Branch("lumi",&_lumi,"lumi/I");
  
  myTree->Branch("stage2_tauN",&_stage2_tauN,"Stage2tausNumber/I");
  myTree->Branch("stage2_tauEt",&_stage2_tauEt);
  myTree->Branch("stage2_tauEta",&_stage2_tauEta);
  myTree->Branch("stage2_tauPhi",&_stage2_tauPhi);
  myTree->Branch("stage2_tauIso",&_stage2_tauIso);
  
  myTree->Branch("stage2_jetN",&_stage2_jetN,"Stage2jetsNumber/I");
  myTree->Branch("stage2_jetEt",&_stage2_jetEt);
  myTree->Branch("stage2_jetEta",&_stage2_jetEta);
  myTree->Branch("stage2_jetPhi",&_stage2_jetPhi);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1ntuples::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1ntuples::Initialize(){
  _indexevents=0;
  _runNumber=0;
  _lumi=0;
  _stage2_tauN = 0;
  _stage2_tauEt.clear();
  _stage2_tauEta.clear();
  _stage2_tauPhi.clear();
  _stage2_tauIso.clear();
  _stage2_jetN = 0;
  _stage2_jetEt.clear();
  _stage2_jetEta.clear();
  _stage2_jetPhi.clear();
  
}

int* L1ntuples::FillStage2(const BXVector<l1t::Tau>* taus, const BXVector<l1t::Jet>* jets, const edm::Event& event){
  static Int_t nObj[4];
  for(int ii =0; ii<4; ii++){
    nObj[ii]=0;
  }
  Initialize();
  //  for (int ibx = taus->getFirstBX(); ibx <= taus->getLastBX(); ++ibx)
  //  {
  int ibx = 0;
      for (BXVector<l1t::Tau>::const_iterator it=taus->begin(ibx); it!=taus->end(ibx); it++)
	{
          if (it->pt() > 0){
            nObj[0]++;
            _stage2_tauEt.push_back(it->et());
            _stage2_tauEta.push_back(it->eta());
            _stage2_tauPhi.push_back(it->phi());
	    _stage2_tauIso.push_back(it->hwIso());
	    if (ibx != 0) cout<<"argh"<<endl;
	   
          }
	}
      // }
      //for (int ibx = jets->getFirstBX(); ibx <= jets->getLastBX(); ++ibx)
      //  {
      for (BXVector<l1t::Jet>::const_iterator it=jets->begin(ibx); it!=jets->end(ibx); it++)
        {
          if (it->pt() > 0){
            nObj[2]++;
            _stage2_jetEt.push_back(it->et());
            _stage2_jetEta.push_back(it->eta());
            _stage2_jetPhi.push_back(it->phi());
          }
        }
      //  }
  return nObj;
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1ntuples);
