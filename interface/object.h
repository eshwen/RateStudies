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


using namespace std;

class object{
public:
  object(float Et, float Eta, float Phi, int Iso){Et_=Et; Eta_ = Eta;Phi_ = Phi;Iso_ = Iso;};
  float Et() const{
    return Et_;
  }
  float Phi() const {
    return Phi_;
  }
  float Eta() const {
    return Eta_;
  }

  int Iso() const {
    return Iso_;
  }
  
  float DeltaR(const object& objref){
      float dEta = objref.Eta()-this->Eta();
      
      float dPhi =  objref.Phi()-this->Phi();
      if(dPhi > 3.14) dPhi -= 6.28;
      if(dPhi < -3.14) dPhi += 6.28;
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      return dR;
    }
    
  bool inline operator <(const object& objref) const {
    return this->Et() > objref.Et();
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
  
  friend std::ostream& operator << (std::ostream &out, const object &object);
private:
  float Et_, Eta_, Phi_, Iso_;
};


std::ostream & operator << (std::ostream &out, const object &object)
{
  out<<object.Et();
  return out;
}
