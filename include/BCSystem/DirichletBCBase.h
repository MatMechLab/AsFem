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
//+++ Date   : 2021.10.06
//+++ Purpose: define the abstract class for dirichlet type boundary
//+++          condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Utils/MessagePrinter.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"
#include "Utils/MathFuns.h"

#include "ElmtSystem/LocalElmtData.h"
#include "FESystem/FECalcType.h"

#include "petsc.h"

using namespace std;

/**
 * This class defines the abstract class for dirichlet type boundary boundary condition
 */
class DirichletBCBase{
public:
    /**
     * execute the boundary condition value for different (dirichlet type) boundary conditions.
     * \f$u=u_{g}\f$ is the final output
     * @param calctype the calculation type for either residual or jacbobian
     * @param bcvalue the boundary value defined in the input file
     * @param params the boundary condition related parameters
     * @param elmtinfo the basic information for current element
     * @param dofid the id of the applied dof(start from 0, the global one)
     * @param nodecoords the coordinate of current node
     * @param K the system jacbobian matrix(to be applied with, i.e., dirichlet bc)
     * @param RHS the system residual
     * @param U the system solution
     */
    virtual void ComputeBCValue(const FECalcType &calctype,const double &bcvalue,const vector<double> &params,const LocalElmtInfo &elmtinfo,const vector<int> &dofids,const Vector3d &nodecoords,Mat &K,Vec &RHS,Vec &U)=0;

    /**
     * calculate the 'displacement' value of current dofs
     * \f$\int_{\partial\Omega}R_{u}dS\f$ is the final outpu
     * @param bcvalue the boundary value defined in the input file
     * @param params the boundary condition related parameters
     * @param elmtinfo the basic information for current element
     * @param nodecoords the coordinate of current node
     * @param localU the solution vector of current node
     */
    virtual void ComputeU(const vector<int> &dofids,const double &bcvalue,const vector<double> &params,const LocalElmtInfo &elmtinfo,const Vector3d &nodecoords,VectorXd &localU)=0; 

protected:
    VectorXd localU;

};


