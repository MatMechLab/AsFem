//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "DofHandler/DofHandler.h"

DofHandler::DofHandler(){
    _nMaxDim=0;_nMinDim=0;
    _nMaxDofsPerNode=0;
    _nMaxDofsPerElmt=0;// for bulk elmt
    _nDofsPerSurfaceElmt=0;_nDofsPerLineElmt=0;
    _nDofs=0;_nActiveDofs=0;
    _nBulkElmts=0;_nElmts=0;_nNodes=0;
    _DofNameList.clear();
    _DofIndexList.clear();
    _DofNameToIDMap.clear();
    _DofIDToNameMap.clear();
}
//*********************************************