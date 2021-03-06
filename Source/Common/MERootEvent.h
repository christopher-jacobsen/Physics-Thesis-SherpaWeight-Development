////////////////////////////////////////////////////////////////////////////////////////////////////
//  MERootEvent.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ME_ROOT_EVENT_H
#define ME_ROOT_EVENT_H

#include <Rtypes.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

class TTree;

////////////////////////////////////////////////////////////////////////////////////////////////////

struct MERootEvent
{
    Int_t       id   = 0;
    Double_t    me   = 0;

    
    void SetInputTree(  TTree * pTree );
    void SetOutputTree( TTree * pTree );
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // ME_ROOT_EVENT_H
