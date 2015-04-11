////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaRootEvent.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaRootEvent.h"

#include "common.h"
#include <sstream>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Root include files

#include <TTree.h>
#include <TBranch.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
static bool AttachBranchToVariable( TTree * pTree, const char * branchName, T & variable )
{
    TBranch * pBranch = pTree->GetBranch( branchName );
    if (!pBranch)
    {
        LogMsgError( "Branch %hs does not exist in tree.", FMT_HS(branchName) );
        return false;
    }
    
    pBranch = nullptr;
    
    pTree->SetBranchStatus(  branchName, 1 );                       // enable branch
    pTree->SetBranchAddress( branchName, &variable, &pBranch );     // connect variable
    pTree->AddBranchToCache( pBranch );                             // cache branch
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
static bool CreateBranchForVariable( TTree * pTree, const char * branchName, T * pVariable, const char * pLeafList = nullptr )
{
    if (!pLeafList)
        pTree->Branch( branchName, pVariable );                 // derive size from type T
    else
        pTree->Branch( branchName, pVariable, pLeafList );      // derive size from pLeafList content
    return true;
}

template<typename T>
static bool CreateBranchForVariable( TTree * pTree, const char * branchName, T & variable, const char * pLeafList = nullptr  )
{
    return CreateBranchForVariable(pTree, branchName, &variable, pLeafList );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEvent::SetInputTree( TTree * pTree )
{
    pTree->SetMakeClass(1);
    
    // disable all input branches as default
    pTree->SetBranchStatus( "*", 0 );
    
    AttachBranchToVariable( pTree, "id",          id          );
    AttachBranchToVariable( pTree, "nparticle",   nparticle   );
    AttachBranchToVariable( pTree, "px",          px          );
    AttachBranchToVariable( pTree, "py",          py          );
    AttachBranchToVariable( pTree, "pz",          pz          );
    AttachBranchToVariable( pTree, "E",           E           );
    AttachBranchToVariable( pTree, "alphas",      alphas      );
    AttachBranchToVariable( pTree, "kf",          kf          );
    AttachBranchToVariable( pTree, "weight",      weight      );
    AttachBranchToVariable( pTree, "weight2",     weight2     );
    AttachBranchToVariable( pTree, "me_wgt",      me_wgt      );
    AttachBranchToVariable( pTree, "me_wgt2",     me_wgt2     );
    AttachBranchToVariable( pTree, "x1",          x1          );
    AttachBranchToVariable( pTree, "x2",          x2          );
    AttachBranchToVariable( pTree, "x1p",         x1p         );
    AttachBranchToVariable( pTree, "x2p",         x2p         );
    AttachBranchToVariable( pTree, "id1",         id1         );
    AttachBranchToVariable( pTree, "id2",         id2         );
    AttachBranchToVariable( pTree, "fac_scale",   fac_scale   );
    AttachBranchToVariable( pTree, "ren_scale",   ren_scale   );
    AttachBranchToVariable( pTree, "nuwgt",       nuwgt       );
    AttachBranchToVariable( pTree, "usr_wgts",    usr_wgts    );
    AttachBranchToVariable( pTree, "alphasPower", alphasPower );
    AttachBranchToVariable( pTree, "part",        part        );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEvent::SetOutputTree( TTree * pTree )
{
    CreateBranchForVariable( pTree, "id",          id           );
    CreateBranchForVariable( pTree, "nparticle",   nparticle    );
    CreateBranchForVariable( pTree, "px",          px,          "px[nparticle]/F" );
    CreateBranchForVariable( pTree, "py",          py,          "py[nparticle]/F" );
    CreateBranchForVariable( pTree, "pz",          pz,          "pz[nparticle]/F" );
    CreateBranchForVariable( pTree, "E",           E,           "E[nparticle]/F"  );
    CreateBranchForVariable( pTree, "alphas",      alphas       );
    CreateBranchForVariable( pTree, "kf",          kf,          "kf[nparticle]/I" );
    CreateBranchForVariable( pTree, "weight",      weight       );
    CreateBranchForVariable( pTree, "weight2",     weight2      );
    CreateBranchForVariable( pTree, "me_wgt",      me_wgt       );
    CreateBranchForVariable( pTree, "me_wgt2",     me_wgt2      );
    CreateBranchForVariable( pTree, "x1",          x1           );
    CreateBranchForVariable( pTree, "x2",          x2           );
    CreateBranchForVariable( pTree, "x1p",         x1p          );
    CreateBranchForVariable( pTree, "x2p",         x2p          );
    CreateBranchForVariable( pTree, "id1",         id1          );
    CreateBranchForVariable( pTree, "id2",         id2          );
    CreateBranchForVariable( pTree, "fac_scale",   fac_scale    );
    CreateBranchForVariable( pTree, "ren_scale",   ren_scale    );
    CreateBranchForVariable( pTree, "nuwgt",       nuwgt        );
    CreateBranchForVariable( pTree, "usr_wgts",    usr_wgts,    "usr_wgts[nuwgt]/D" );
    CreateBranchForVariable( pTree, "alphasPower", alphasPower  );
    CreateBranchForVariable( pTree, "part",        part,        "part[2]/C" );
}
