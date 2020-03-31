//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "DofHandler/DofHandler.h"

void DofHandler::AddNameToDofNameList(vector<string> dofnames){
    _DofNameList.clear();
    _DofNameToIDMap.clear();
    _DofIDToNameMap.clear();
    PetscInt k=0;
    for(unsigned int i=0;i<dofnames.size();i++){
        _DofNameList.push_back(dofnames[i]);
        k+=1;
        _DofIDToNameMap.push_back(make_pair(k,dofnames[i]));
        _DofNameToIDMap.push_back(make_pair(dofnames[i],k));
    }
    _nMaxDofsPerNode=int(dofnames.size());
}