#ifndef g4PSIEvent_h
#define g4PSIEvent_h 1

//---------------------------------------------------------------------------
//!
//! ClassName:   g4PSIEvent
//!
//----------------------------------------------------------------------------
//

#include <TVector.h>

class g4PSIEvent {
public:
  void reset() {;};
  void add(TVector p) {
  };
private:
  Double_t p_;
  Int_t nhit_;
  std::vector<TVector> cos_;
};

#endif
