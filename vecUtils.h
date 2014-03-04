
#include "Math/LorentzVector.h"

#include "Math/VectorUtil.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > LV;
namespace ROOT {
  namespace Math {
    namespace VectorUtil {
      float DeltaPhi     ( const LV& v1, const LV& v2) {return DeltaPhi     <LV,LV>(v1,v2);}
      float DeltaR       ( const LV& v1, const LV& v2) {return DeltaR       <LV,LV>(v1,v2);}
      float CosTheta     ( const LV& v1, const LV& v2) {return CosTheta     <LV,LV>(v1,v2);}
      float Angle        ( const LV& v1, const LV& v2) {return Angle        <LV,LV>(v1,v2);}
      float InvariantMass( const LV& v1, const LV& v2) {return InvariantMass<LV,LV>(v1,v2);}
    }
  }
}

#ifdef __CINT__
#pragma link C++ class ROOT::Math::PtEtaPhiE4D<float>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >+;
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >::operator+(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >::operator-(LV);
#pragma link C++ function ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >::Dot(LV);
#endif
