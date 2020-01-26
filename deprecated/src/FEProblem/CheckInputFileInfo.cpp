#include "FEProblem/FEProblem.h"

void FEProblem::CheckInputFileInfo()
{
    // first check elmts settings
    //  here we mainly check the dofs is the same as or not as in the [dofs] block
    string dofname,blockname;
    int n;
    bool DofInvalid,BlockInvalid;
    vector<int> DofAssignFlag;
    DofAssignFlag.resize(dofHandler.GetDofsNumPerNode(),0);
    for(int i=1;i<=elmtSystem.GetElmtBlockNum();i++)
    {
        n=elmtSystem.GetIthElmtBlock(i)._ElmtDofNameList.size();
        blockname=elmtSystem.GetIthElmtBlock(i)._ElmtDomainBlockName;
        for(int j=0;j<n;++j)
        {
            DofInvalid=true;
            dofname=elmtSystem.GetIthElmtBlock(i)._ElmtDofNameList[j];
            if(dofHandler.GetDofIndexViaName(dofname)){
                DofAssignFlag[dofHandler.GetDofIndexViaName(dofname)-1]=1;
            }
            for(int k=1;k<=dofHandler.GetDofNameListLenght();k++)
            {
                if(dofHandler.GetIthDofNameFromList(k)==dofname)
                {
                    DofInvalid=false;
                    break;
                }
            }
            if(DofInvalid)
            {
                Msg_InputCheck_DofInvalidInElmtBlock(j+1,elmtSystem.GetIthElmtBlock(i)._ElmtBlockName);
                Msg_AsFem_Exit();
            }
        }
        BlockInvalid=true;
        for(auto it:mesh.GetBulkGroupList()){
            if(it==blockname){
                BlockInvalid=false;
                break;
            }
        }
        if(BlockInvalid)
        {
            Msg_ElmtBlockBlockNameInvalid(elmtSystem.GetIthElmtBlock(i)._ElmtBlockName);
            Msg_AsFem_Exit();
        }
    }
    // check whether all the dofs are applied to elmts
    for(unsigned int i=0;i<DofAssignFlag.size();++i){
        if(DofAssignFlag[i]==0){
            dofname=dofHandler.GetIthDofNameFromList(i+1);
            Msg_InputCheck_DofIsNotAssignedToElmt(dofname);
            Msg_AsFem_Exit();
            break;
        }
    }
    // cross-check between element blocks
    string iblockname,jblockname;
    string iblockdomain,jblockdomain;
    for(int i=1;i<=elmtSystem.GetElmtBlockNum();i++){
        iblockdomain=elmtSystem.GetIthElmtBlock(i)._ElmtDomainBlockName;
        iblockname=elmtSystem.GetIthElmtBlock(i)._ElmtBlockName;
        for(int j=1;j<=elmtSystem.GetElmtBlockNum();j++){
            jblockdomain=elmtSystem.GetIthElmtBlock(j)._ElmtDomainBlockName;
            jblockname=elmtSystem.GetIthElmtBlock(j)._ElmtBlockName;
            if(i!=j){
                if(iblockdomain==jblockdomain){
                    Msg_InputCheck_ElmtBlockDuplicateDomain(jblockname);
                    Msg_AsFem_Exit();
                }
            }
        }
    }
    //************************************************
    // Check the boundary condition block
    //************************************************
    for(int i=1;i<=bcSystem.GetBCBlockNum();i++)
    {
        
        DofInvalid=true;
        dofname=bcSystem.GetIthBCBlock(i)._DofName;
        for(int j=1;j<=dofHandler.GetDofNameListLenght();++j)
        {
            if(dofHandler.GetIthDofNameFromList(j)==dofname)
            {
                DofInvalid=false;
                break;
            }
        }
        if(DofInvalid)
        {
            Msg_InputCheck_DofInvalidInBCBlock(bcSystem.GetIthBCBlock(i)._BCBlockName);
            Msg_AsFem_Exit();
        }
        BlockInvalid=true;
        blockname=bcSystem.GetIthBCBlock(i)._BoundaryName;
        for(auto it:mesh.GetBounGroupList())
        {
            if(it==blockname)
            {
                BlockInvalid=false;
                break;
            }
        }
        if(BlockInvalid)
        {
            Msg_BCBlockBlockNameInvalid(bcSystem.GetIthBCBlock(i)._BCBlockName);
            Msg_AsFem_Exit();
        }
    }
    //************************************************
    // Check the initial condition block
    //************************************************
    for(int i=1;i<=icSystem.GetICBlockNum();i++)
    {
        DofInvalid=true;
        dofname=icSystem.GetIthICBlock(i)._DofName;
        for(int j=1;j<=dofHandler.GetDofNameListLenght();++j)
        {
            if(dofHandler.GetIthDofNameFromList(j)==dofname)
            {
                DofInvalid=false;
                break;
            }
        }
        if(DofInvalid)
        {
            Msg_InputCheck_DofInvalidInICBlock(icSystem.GetIthICBlock(i)._ICBlockName);
            Msg_AsFem_Exit();
        }

        BlockInvalid=true;
        blockname=icSystem.GetIthICBlock(i)._BlockName;
        for(auto it:mesh.GetBulkGroupList())
        {
            if(it==blockname)
            {
                BlockInvalid=false;
                break;
            }
        }
        if(BlockInvalid)
        {
            Msg_ICBlockBlockNameInvalid(icSystem.GetIthICBlock(i)._ICBlockName);
            Msg_AsFem_Exit();
        }
    }

    //************************************************
    // Check the material settings is compatible with the elmts or not
    //************************************************
    string mateblockname,elmtblockname,elmtmatename,elmttypename;
    bool IsMateMatch;
    for(int i=1;i<=elmtSystem.GetElmtBlockNum();i++)
    {
        IsMateMatch=false;
        elmtmatename=elmtSystem.GetIthElmtBlock(i)._ElmtMateName;
        elmtblockname=elmtSystem.GetIthElmtBlock(i)._ElmtBlockName;
        elmttypename=elmtSystem.GetIthElmtBlock(i)._ElmtTypeName;
        if(elmtmatename.length()<1) 
        {
            // no material is given, then we use the default material settings
            // if(elmtSystem.GetIthElmtBlock(i)._ElmtTypeID==ElmtType::POISSON){
            //     elmtSystem.SetIthElmtBlockMateID(i,MaterialType::CONSTPOISSONMATE);
            //     IsMateMatch=true;
            //     continue;
            // }
            // else if(elmtSystem.GetIthElmtBlock(i)._ElmtTypeID==ElmtType::DIFFUSION){
            //     elmtSystem.SetIthElmtBlockMateID(i,MaterialType::CONSTDIFFUSIONMATE);
            //     IsMateMatch=true;
            //     continue;
            // }
            // else{
            //     Msg_InputCheck_MateNotFoundForElmtBlock(elmtmatename,elmtblockname);
            //     Msg_AsFem_Exit();
            // }
            Msg_InputCheck_MateNotFoundForElmtBlock(elmtmatename,elmtblockname);
            Msg_AsFem_Exit();
        }
        // else if(elmtmatename=="nonlinear"&&elmttypename=="poisson"){
        //     elmtSystem.SetIthElmtBlockMateID(i,MaterialType::NONLINEARPOISSONMATE);
        //     IsMateMatch=true;
        //     continue;
        // }
        for(int j=1;j<=materialSystem.GetMateBlockNum();j++)
        {
            mateblockname=materialSystem.GetIthMateBlock(j)._MateBlockName;
            if(elmtmatename==mateblockname)
            {
                elmtSystem.SetIthElmtBlockMateID(i,materialSystem.GetIthMateBlock(j)._MateTypeID);
                IsMateMatch=true;
                break;
            }
        }
        if(!IsMateMatch)
        {
            Msg_InputCheck_MateNotFoundForElmtBlock(elmtmatename,elmtblockname);
            Msg_AsFem_Exit();
        }
    }
}