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

void ProjectionSystem::runProjectionLibs(const bool &flag,
                                         const FECell &t_fecell,
                                         const int &nodesnum,
                                         const vector<int> &t_elconn,
                                         const double &detjac,
                                         const ShapeFun &t_shp,
                                         const MaterialsContainer &t_mate,
                                         ProjectionData &t_data){
    switch (m_proj_type)
    {
    case ProjectionType::DEFAULT:
    case ProjectionType::LEASTSQUARE:{
        if(flag){
            LeastSquareProjection::localProjectionAction(nodesnum,t_elconn,detjac,t_shp,t_mate,t_data);
        }
        else{
            LeastSquareProjection::globalProjectionAction(t_fecell,t_data);
        }
        break;
    }
    default:
        MessagePrinter::printErrorTxt("Unsupported projection method, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}