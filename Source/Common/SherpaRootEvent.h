////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaRootEvent.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPA_ROOT_EVENT_H
#define SHERPA_ROOT_EVENT_H

#include <Rtypes.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

class TTree;

////////////////////////////////////////////////////////////////////////////////////////////////////

struct SherpaRootEvent
{
    static const size_t max_nparticle   = 100;  // See Sherpa/AddOns/Root/Output_RootNtuple.H
    static const size_t max_nuwgt       = 18;   // See Sherpa/AddOns/Root/Output_RootNtuple.H

    Int_t       id                      = 0;

    Int_t       nparticle               = 0;

    Float_t     px[max_nparticle]       = {};   // [nparticle]
    Float_t     py[max_nparticle]       = {};   // [nparticle]
    Float_t     pz[max_nparticle]       = {};   // [nparticle]
    Float_t     E [max_nparticle]       = {};   // [nparticle]

    Double_t    alphas                  = 0;

    Int_t       kf[max_nparticle]       = {};   // [nparticle]

    Double_t    weight                  = 0;
    Double_t    weight2                 = 0;
    Double_t    me_wgt                  = 0;
    Double_t    me_wgt2                 = 0;

    Double_t    x1                      = 0;
    Double_t    x2                      = 0;
    Double_t    x1p                     = 0;
    Double_t    x2p                     = 0;

    Int_t       id1                     = 0;
    Int_t       id2                     = 0;

    Double_t    fac_scale               = 0;
    Double_t    ren_scale               = 0;

    Int_t       nuwgt                   = 0;
    
    Double_t    usr_wgts[max_nuwgt]     = {};   // [nuwgt]

    Char_t      alphasPower             = 0;
    Char_t      part[2]                 = {};
  
public:

    void SetInputTree(  TTree * pTree );
    void SetOutputTree( TTree * pTree );
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPA_ROOT_EVENT_H
