//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_PPSSYSTEM_H
#define ASFEM_PPSSYSTEM_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <algorithm>

#include "petsc.h"


class PPSSystem{
public:
    PPSSystem();

    void RunPostProcess();

private:
    void VolumeIntegrationPostProcess();
    void SurfaceIntegrationPostProcess();

    void NodeValuePostProcess();
    void SurfaceValuePostProcess();
    void VolumeValuePostProcess();

private:
    int _nPPS;
    int _nMaxDim,_nMinDim;
};

#endif // ASFEM_PPSSYSTEM_H