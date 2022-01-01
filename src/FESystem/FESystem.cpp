//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.30
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

FESystem::FESystem(){
    _BulkVolumes=0.0;
    
    _elConn.clear();_elDofs.clear();
    _elDofsActiveFlag.clear();
    _elU.clear();_elV.clear();
    _gpU.clear();_gpV.clear();
    _gpHist.clear();_gpHistOld.clear();_gpProj.clear();
    _gpGradU.clear();_gpGradV.clear();
    _MaterialValues.clear();
    _nHist=0;_nProj=0;
    _MaxKMatrixValue=-1.0e3;_KMatrixFactor=0.1;

    _localK.Clean();_localR.Clean();
}