////////////////////////////////////////////////////////////////////////////////////////////////////
//  MERootEvent.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "MERootEvent.h"

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Root include files

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"

#include <TTree.h>
#include <TBranch.h>

#pragma clang diagnostic pop

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

////////////////////////////////////////////////////////////////////////////////////////////////////
void MERootEvent::SetInputTree( TTree * pTree )
{
    pTree->SetMakeClass(1);
    
    // disable all input branches as default
    pTree->SetBranchStatus( "*", 0 );
    
    AttachBranchToVariable( pTree, "id", id );
    AttachBranchToVariable( pTree, "me", me );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void MERootEvent::SetOutputTree( TTree * pTree )
{
    CreateBranchForVariable( pTree, "id", id );
    CreateBranchForVariable( pTree, "me", me );
}
