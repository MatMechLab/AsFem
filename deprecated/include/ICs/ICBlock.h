#ifndef ASFEM_ICBLOCK_H
#define ASFEM_ICBLOCK_H


#include <iostream>
#include <iomanip>
#include <string>

#include "ICType.h"

class ICBlock
{
public:
    std::string _ICBlockName;    // [ic1]
    std::string _DofName;        //  dof=u1
    std::string _BlockName;   //  block=left
    std::string _ICElmtName;     //  type=dirichlet [user1]
    std::vector<double> _value;
    ICType _ICTypeID;

    inline void Clear()
    {
        _ICBlockName="";
        _DofName="";
        _BlockName="";
        _ICElmtName=""; 
        _value.clear();
        _ICTypeID=ICType::CONST;
    }

    void PrintICBlockInfo() const
    {
        std::printf("*** IC block name= %-34s      ***\n",_ICBlockName.c_str());
        std::printf("***    type name = %-34s      ***\n",_ICElmtName.c_str());
        std::printf("***    dof name  = %-34s      ***\n",_DofName.c_str());
        if(_value.size()>0)
        {
            std::printf("***    value     =");
            for(unsigned int i=0;i<_value.size();i++)
            {
                std::cout<<_value[i]<<" ";
            }
            std::cout<<std::endl;
        }
        std::printf("***    block     = %-34s      ***\n",_BlockName.c_str());
    }
};


#endif 