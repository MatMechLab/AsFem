//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "TimeStepping/TimeStepping.h"


TimeStepping::TimeStepping(){
    _FinalTime=0.1;
    _dtmax=0.1;_dtmin=1.0e-13;
    _dt0=1.0e-5;
    _interval=1;
    _IsAdaptive=false;
    _CutFactor=0.85;_GrowthFactor=1.1;
    _nOpts=4;
}