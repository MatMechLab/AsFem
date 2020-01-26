#ifndef ASFEM_ELMTBLOCK_H
#define ASFEM_ELMTBLOCK_H

#include <iostream>
#include <vector>

#include "ElmtType.h"
#include "MaterialSystem/MaterialType.h"

class ElmtBlock
{
public:
    //  [solid]
    //    type=solid2d
    //    dofs=u1 u2
    //    mate=mate1 [can be ignored]
    //    block=all  [can be ignored]
    //  [end]
    std::string _ElmtBlockName="";          //  [solid]
    std::string _ElmtTypeName="";           //    type=solid2d
    std::string _ElmtMateName="";           //    mate=mate1 [can be ignored]
    std::string _ElmtDomainBlockName="alldomain";       //    block=all  [can be ignored]
    std::vector<std::string> _ElmtDofNameList;
    std::vector<int> _ElmtDofIndexList;
    ElmtType _ElmtTypeID=ElmtType::NULLELMT;
    MaterialType _MateTypeID=MaterialType::NULLMATE;


    inline void Clean()
    {
        _ElmtBlockName="";
        _ElmtTypeName="";
        _ElmtMateName="";
        _ElmtDomainBlockName="alldomain";
        _ElmtDofNameList.clear();
        _ElmtDofIndexList.clear();
        _ElmtTypeID=ElmtType::NULLELMT;
        _MateTypeID=MaterialType::NULLMATE;
    }

    void PrintInfo() const
    {
        std::printf("*** Elmt block name= %-32s      ***\n",_ElmtBlockName.c_str());
        std::printf("***    type name = %-34s      ***\n",_ElmtTypeName.c_str());
        std::cout<<"***    dof name  = ";
        for(unsigned int i=0;i<_ElmtDofNameList.size();++i)
        {
            std::cout<<_ElmtDofNameList[i]<<" ";
        }
        std::cout<<std::endl;
        std::cout<<"***    dof index  = ";
        for(unsigned int i=0;i<_ElmtDofIndexList.size();++i)
        {
            std::cout<<_ElmtDofIndexList[i]<<" ";
        }
        std::cout<<std::endl;
        std::printf("***    material  = %-34s      ***\n",_ElmtMateName.c_str());
        std::printf("***    block     = %-34s      ***\n",_ElmtDomainBlockName.c_str());
    }
};

#endif 