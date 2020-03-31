//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/FE.h"


FE::FE(){
    _IsInit=false;
    _nDim=0;_nDimMin=0;_nOrder=2;_nBCOrder=2;
    _QPointType="gauss";
}