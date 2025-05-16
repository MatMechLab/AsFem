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
//+++ Date   : 2022.07.22
//+++ Purpose: Defines the abstract class for general projection,
//+++          the projection is used for extrapolate the guass point
//+++          quantities to the nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "Utils/MessagePrinter.h"

#include "ProjectionSystem/ProjectionData.h"
#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "FE/ShapeFun.h"
#include "MateSystem/MateSystem.h"
#include "MateSystem/MaterialsContainer.h"
#include "FE/FE.h"
#include "SolutionSystem/SolutionSystem.h"
#include "FEProblem/FEControlInfo.h"

using std::string;
using std::vector;


/**
 * The abstract class for general projection
 */
class ProjectionBase{
protected:
    /**
     * initialize the projection system
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof handler class
     */
    virtual void initMyProjection(const FECell &t_FECell,const DofHandler &t_DofHandler)=0;

    /**
     * execute my own projection method based on the child class
     * @param t_FECell the FECell class
     * @param t_DofHandler the DofHandler class
     * @param t_ElmtSystem the ElmtSystem class
     * @param t_MateSystem the MateSystem class
     * @param t_FE the FE class
     * @param t_SolnSystem the SolutionSystem
     * @param t_FECtrlInfo the FECtrlInfo structure
     * @param Data the projection data
     */
    virtual void executeMyProjection(const FECell &t_FECell,
                                     const DofHandler &t_DofHandler,
                                     const ElmtSystem &t_ElmtSystem,
                                     MateSystem &t_MateSystem,
                                     FE &t_FE,
                                     SolutionSystem &t_SolnSystem,
                                     const FEControlInfo &t_FECtrlInfo,
                                     ProjectionData &Data)=0;

};