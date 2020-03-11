//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalHistToGlobal(const int &elmtid,const int &gpInd,const vector<double> &gpHist,Vec &Hist){
    int i,iInd;
    for(i=1;i<=_nHist;i++){
        iInd=(elmtid-1)*(_nHist*_nGPoints)+(gpInd-1)*_nHist+i-1;
        VecSetValue(Hist,iInd,gpHist[i-1],INSERT_VALUES);
    }
}