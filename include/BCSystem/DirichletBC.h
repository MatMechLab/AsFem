//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.10.06
//+++ Purpose: implement the dirichlet boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "BCSystem/DirichletBCBase.h"

/**
 * This class implement the dirichelt boundary condition calculation
 */
class DirichletBC:public DirichletBCBase{
public:
    DirichletBC();

    /**
     * This function calculate the 'displacement' value for dirichlet boundary condition,
     * and reset the K and RHS with the penalty method
     * @param calctype the FEM calculation type, it could be jacobian and residual
     * @param bcvalue the boundary condition value to be applied
     * @param params the parameters taken from the input file
     * @param elmtinfo the local element information
     * @param dofids the ids of each dof
     * @param nodecoords the coordinates of current node, its node , not the gauss point !!!
     * @param K the K matrix of the system
     * @param RHS the residual of the system
     * @param U the solution vector of the system
     */ 
    void ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const vector<int> &dofids, const Vector3d &nodecoords, Mat &K, Vec &RHS, Vec &U) override;

private:
    /**
     * This function calculate the 'displacement' value for dirichlet boundary condition.
     * @param dofids the id of dofs
     * @param bcvalue the boundary condition value to be applied
     * @param params the parameters taken from the input file
     * @param elmtinfo the local element information
     * @param nodecoords the coordinates of current node, its node , not the gauss point !!!
     * @param localU the local solution vector
     */ 
    void ComputeU(const vector<int> &dofids,const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const Vector3d &nodecoords, VectorXd &localU) override;

};

