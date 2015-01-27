////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaEvent.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaEvent.h"

#include "common.h"
#include <sstream>

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
static bool CreateBranchForVariable( TTree * pTree, const char * branchName, T * pVariable ) //, const char * pLeafList = nullptr )
{
    pTree->Branch( branchName, pVariable );
    return true;
}

template<typename T>
static bool CreateBranchForVariable( TTree * pTree, const char * branchName, T & variable )
{
    return CreateBranchForVariable(pTree, branchName, &variable );
}
/*
template<typename T, int n>
static bool CreateBranchForVariable( TTree * pTree, const char * branchName, T variable[n] )
{
    std::stringstream stream;
    stream << branchName << '[' << n << ']';
    return CreateBranchForVariable(pTree, branchName, variable, stream.str().c_str() );
}*/


////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaEvent::SetInputTree( TTree * pTree )
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
void SherpaEvent::SetOutputTree( TTree * pTree )
{
    CreateBranchForVariable( pTree, "id",          id          );
    CreateBranchForVariable( pTree, "nparticle",   nparticle   );
    CreateBranchForVariable( pTree, "px[4]",       px          );
    CreateBranchForVariable( pTree, "py[4]",       py          );
    CreateBranchForVariable( pTree, "pz[4]",       pz          );
    CreateBranchForVariable( pTree, "E[4]",        E           );
    CreateBranchForVariable( pTree, "alphas",      alphas      );
    CreateBranchForVariable( pTree, "kf[4]",       kf          );
    CreateBranchForVariable( pTree, "weight",      weight      );
    CreateBranchForVariable( pTree, "weight2",     weight2     );
    CreateBranchForVariable( pTree, "me_wgt",      me_wgt      );
    CreateBranchForVariable( pTree, "me_wgt2",     me_wgt2     );
    CreateBranchForVariable( pTree, "x1",          x1          );
    CreateBranchForVariable( pTree, "x2",          x2          );
    CreateBranchForVariable( pTree, "x1p",         x1p         );
    CreateBranchForVariable( pTree, "x2p",         x2p         );
    CreateBranchForVariable( pTree, "id1",         id1         );
    CreateBranchForVariable( pTree, "id2",         id2         );
    CreateBranchForVariable( pTree, "fac_scale",   fac_scale   );
    CreateBranchForVariable( pTree, "ren_scale",   ren_scale   );
    CreateBranchForVariable( pTree, "nuwgt",       nuwgt       );
    CreateBranchForVariable( pTree, "usr_wgts[1]", usr_wgts    );
    CreateBranchForVariable( pTree, "alphasPower", alphasPower );
    CreateBranchForVariable( pTree, "part[2]",     part        );
}
