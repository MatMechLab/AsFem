//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.08.22
//+++ Purpose: The projection libs for different projection method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/ProjectionSystem.h"

void ProjectionSystem::runProjectionLibs(const bool &Flag,
                                         const FECell &t_FECell,
                                         const int &NodesNum,
                                         const vector<int> &ElConn,
                                         const double &DetJac,
                                         const ShapeFun &Shp,
                                         const MaterialsContainer &Mate,
                                         ProjectionData &Data){
    switch (m_ProjType)
    {
    case ProjectionType::DEFAULT:
    case ProjectionType::LEASTSQUARE:{
        if(Flag){
            LeastSquareProjection::localProjectionAction(NodesNum,ElConn,DetJac,Shp,Mate,Data);
        }
        else{
            LeastSquareProjection::globalProjectionAction(t_FECell,Data);
        }
        break;
    }
    default:
        MessagePrinter::printErrorTxt("Unsupported projection method, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}