//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.10
//+++ Purpose: Define the initial condition system in AsFem
//+++          Here we can apply:
//+++               1) constant value ic
//+++               2) random value ic
//+++               3) other type or user defined ic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


//******************************************
//*** for AsFem own header
//******************************************
#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

#include "MathUtils/Vector.h"

#include "ICSystem/ICBlock.h"
#include "ICSystem/ICType.h"

/**
 * For built-in and user-defined ICs
 */
#include "ICSystem/ConstantIC.h"
#include "ICSystem/RandomIC.h"
#include "ICSystem/CircleIC.h"


/**
 * This class implement the initial condition management for the FEM simulation
 */
class ICSystem:public ConstantIC,
               public RandomIC,
               public CircleIC{
public:
    /**
     * constructor
     */
    ICSystem();
    /**
     * initialize the ICSystem
     * @param the maximum nodal dofs
     */
    void init(const int &nodal_dofs);

    /**
     * add single icblock to the list
     * @param t_icblock the single ic block read from input file
     */
    void addICBlock2List(const ICBlock &t_icblock);

    /**
     * apply the initial conditions to the solution
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofHandler class
     * @param U0 the initial solution
     */
    void applyInitialConditions(const Mesh &t_mesh,const DofHandler &t_dofhandler,Vector &U0);

    /**
     * print out the initial condition info
     */
    void printICSystemInfo()const;
    /**
     * release the allocated memory
     */
    void releaseMemory();

    /**
     * get the number of ic blocks read from input file
     */
    inline int getICBlocksNum()const{return m_icblocks_num;}

private:
    //*******************************************************
    //*** for different initial conditions
    //*******************************************************
    /**
     * compute the initial condition value for each dof (of each node)
     * @param t_type the initial condition type
     * @param t_params parameters read from json input file
     * @param icvalue the initial condition value
     * @param dim the dimensions of current mesh
     * @param dofs the dofs number of current IC
     * @param nodecoords the current node's coordinates
     * @param localU the initial solution values of current IC
     */
    void runICLibs(const ICType &t_type,
                   const nlohmann::json &t_params,
                   const double &icvalue,
                   const int &dim,
                   const int &dofs,
                   const Vector3d &nodecoords,
                   VectorXd &localU);

private:
    vector<ICBlock> m_icblock_list;/**< the ic block list read from input file */
    int m_icblocks_num;/**< number of total ic blocks */

    VectorXd m_localU;/**< the local initial solution vector */

private:
    PetscMPIInt m_rank;/**< for the local process id */
    PetscMPIInt m_size;/**< for the local process id */

    PetscRandom m_rnd;/**< random seed */
    double m_minval,m_maxval,m_value;/**< for the max, min, and intermediate value */

};