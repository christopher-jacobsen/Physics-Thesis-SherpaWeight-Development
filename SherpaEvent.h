////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaEvent.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMEX_SHERPAEVENT_H
#define SMEX_SHERPAEVENT_H

#include <RTypes.h>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

class TTree;

////////////////////////////////////////////////////////////////////////////////////////////////////

struct SherpaEvent
{
    Int_t                       id              = 0;

    Int_t                       nparticle       = 0;

    Float_t                     px[4]           = {};           // [nparticle]
    Float_t                     py[4]           = {};           // [nparticle]
    Float_t                     pz[4]           = {};           // [nparticle]
    Float_t                     E[4]            = {};           // [nparticle]
  //std::vector<Float_t> *      px              = nullptr;      // [nparticle]
  //std::vector<Float_t> *      py              = nullptr;      // [nparticle]
  //std::vector<Float_t> *      pz              = nullptr;      // [nparticle]
  //std::vector<Float_t> *      E               = nullptr;      // [nparticle]

    Double_t                    alphas          = 0;

    Int_t                       kf[4]           = {};           // [nparticle]
  //std::vector<Int_t> *        kf              = nullptr;      // [nparticle]

    Double_t                    weight          = 0;
    Double_t                    weight2         = 0;
    Double_t                    me_wgt          = 0;
    Double_t                    me_wgt2         = 0;

    Double_t                    x1              = 0;
    Double_t                    x2              = 0;
    Double_t                    x1p             = 0;
    Double_t                    x2p             = 0;

    Int_t                       id1             = 0;
    Int_t                       id2             = 0;

    Double_t                    fac_scale       = 0;
    Double_t                    ren_scale       = 0;

    Int_t                       nuwgt           = 0;

    Double_t                    usr_wgts[1]     = {};           // [nuwgt]
  //std::vector<Double_t> *     usr_wgts        = nullptr;      // [nuwgt]

    Char_t                      alphasPower     = 0;
    Char_t                      part[2]         = {};
  
public:

    void SetInputTree(  TTree * pTree );
    void SetOutputTree( TTree * pTree );
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SMEX_SHERPAEVENT_H
