#ifndef ASFEM_MATERIALBLOCK_H
#define ASFEM_MATERIALBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "MaterialType.h"

class MaterialBlock
{
public:
    std::string _MateBlockName;
    std::string _MateTypeName;
    std::vector<double> _MateParams;
    MaterialType _MateTypeID=MaterialType::CONSTPOISSONMATE;

    void Clean()
    {
        _MateBlockName="";
        _MateTypeName="";
        _MateParams.clear();
        _MateTypeID=MaterialType::NULLMATE;
    }

    void PrintMateBlockInfo() const
    {
        std::printf("*** Mateblockname= %-34s      ***\n",_MateBlockName.c_str());
        std::printf("***    type name = %-34s      ***\n",_MateTypeName.c_str());
        if(_MateParams.size()>0)
        {
            std::printf("***    params    =");
            for(unsigned int i=0;i<_MateParams.size();i++)
            {
                std::cout<<_MateParams[i]<<" ";
            }
            std::cout<<std::endl;
        }
    }
};

#endif // ASFEM_MATERIALBLOCK_H