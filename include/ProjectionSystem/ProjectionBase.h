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
#include "MateSystem/MaterialsContainer.h"

using std::string;
using std::vector;


/**
 * The abstract class for general projection
 */
class ProjectionBase{
protected:
    /**
     * the global projection action
     * @param t_fecell the fe cell class
     * @param t_data the projection data structure
     */
    virtual void globalProjectionAction(const FECell &t_fecell,ProjectionData &t_data)=0;
    /**
     * the local projection action
     * @param nodesnum the nodes number of current element
     * @param t_elconn the local element connectivity
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param m_mate the material container from material system
     * @param t_data the projection data structure
     */
    virtual void localProjectionAction(const int &nodesnum,
                                       const vector<int> &t_elconn,
                                       const double &detjac,
                                       const ShapeFun &t_shp,
                                       const MaterialsContainer &m_mate,
                                       ProjectionData &t_data)=0;

};