#ifndef ASFEM_BCBLOCK_H
#define ASFEM_BCBLOCK_H

#include <iostream>
#include <iomanip>
#include <string>



class BCBlock
{
public:
    std::string _BCBlockName;    // [bc1]
    std::string _DofName;        //  dof=u1
    std::string _BoundaryName;   //  boundary=left
    std::string _BCElmtName;     //  type=dirichlet [user1]
    double _value=0.0;
    bool _IsTimeDependent=false;

    inline void Clear()
    {
        _BCBlockName="";
        _DofName="";
        _BoundaryName="";
        _BCElmtName=""; 
        _value=0.0;
        _IsTimeDependent=false;
    }

    void PrintBCBlockInfo() const{
        std::printf("*** BC block name= %-34s      ***\n",_BCBlockName.c_str());
        std::printf("***    type name = %-34s      ***\n",_BCElmtName.c_str());
        std::printf("***    dof name  = %-34s      ***\n",_DofName.c_str());
        std::printf("***    value     = %-16.8e                        ***\n",_value);
        std::printf("***    boundary  = %-34s      ***\n",_BoundaryName.c_str());
        if(_IsTimeDependent){
            std::printf("***    time dependent=true                                 ***\n");
        }
    }
};

#endif //ASFEM_BCBLOCK_H