//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Utils/StringUtils.h"

bool IsUniqueStrVec(vector<string> &strvec){
    if(strvec.size()<=1){
        return true;
    }
    else if(strvec.size()==2){
        if(strvec[0]==strvec[1]){
            return false;
        }
        else{
            return true;
        }
    }
    else if(strvec.size()==3){
        if(strvec[0]!=strvec[1]&&strvec[0]!=strvec[2]&&strvec[1]!=strvec[2]){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        bool IsDuplicate=false;
        for(unsigned int i=0;i<strvec.size();++i){
            for(unsigned int j=0;j<strvec.size();++j){
                if(i!=j){
                    if(strvec[i]==strvec[j]){
                        IsDuplicate=true;
                        break;
                    }
                }
            }
        }
        if(IsDuplicate){
            return false;
        }
        else{
            return true;
        }
    }
    return true;
}