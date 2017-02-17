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
  
  
  friend std::ostream& operator << (std::ostream &out, const object &object);
private:
  float Et_, Eta_, Phi_, Iso_;
};


std::ostream & operator << (std::ostream &out, const object &object)
{
  out<<object.Et();
  return out;
}
