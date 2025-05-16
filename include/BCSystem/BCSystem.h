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
//+++ Date   : 2020.07.10
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) dirichlet bc, i.e. displacement, temperature ...
//+++               2) neuman bc, i.e. flux, force
//+++               3) robin bc as well as user-defined-bc (ubc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * for AsFem's own header files
 */
#include "BCSystem/BCType.h"
#include "BCSystem/BCBlock.h"

#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "FESystem/FECalcType.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/Vector.h"

#include "MathUtils/MatrixXd.h"
#include "MathUtils/SparseMatrix.h"

#include "MathUtils/Rank2Tensor.h"
#include "MathUtils/Rank4Tensor.h"

/**
 * for built-in and user-defined boundary conditions
 */
#include "BCSystem/DirichletBC.h"
#include "BCSystem/User1DirichletBC.h"
#include "BCSystem/User2DirichletBC.h"
#include "BCSystem/User3DirichletBC.h"
#include "BCSystem/User4DirichletBC.h"
#include "BCSystem/User5DirichletBC.h"
#include "BCSystem/Poisson2DBenchmarkBC.h"
#include "BCSystem/RotatedDirichletBC.h"

#include "BCSystem/NeumannBC.h"
#include "BCSystem/PressureBC.h"
#include "BCSystem/TractionBC.h"


/**
 * This class implements the boundary condition management and calculation.
 */
class BCSystem:public DirichletBC,
               public User1DirichletBC,
               public User2DirichletBC,
               public User3DirichletBC,
               public User4DirichletBC,
               public User5DirichletBC,
               public Poisson2DBenchmarkBC,
               public RotatedDirichletBC,
               // for integrated bc
               public NeumannBC,
               public PressureBC,
               public TractionBC{
public:
    /**
     * constructor
     */
    BCSystem();
    /**
     * add bcblock to the list, the bc block must be unique
     * @param t_BCBlock the bc block to be added
     */
    void addBCBlock2List(const BCBlock &t_BCBlock);
    /**
     * set the dirichlet penalty coefficient
     */
    inline void setDirichletPenalty(const double &val){m_DirichletPenalty=val;}

    /**
     * initialize the bc system
     * @param t_dofs max dofs per bc element
     */
    void init(const int &t_dofs);

    //*******************************************************
    //*** general gettings
    //*******************************************************
    /**
     * get the total number of bc blocks defined in the input file
     */
    inline int getBCBlocksNum()const{return m_BCBlocksNum;}
    /**
     * get the i-th bc block
     * @param i the integer index for bc block
     */
    inline BCBlock getIthBCBlock(const int &i)const{
        if(i<1||i>m_BCBlocksNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bcblocks");
            MessagePrinter::exitAsFem();
        }
        return m_BCBlockList[i-1];
    }
    /**
     * get the penalty coefficient for dirichlet bc
     */
    inline double getDirichletPenalty()const{return m_DirichletPenalty;}

    //*******************************************************
    //*** for different boundary conditions
    //*******************************************************
    /**
     * apply the boundary conditions
     * @param CalcType the calculation type
     * @param t current time
     * @param Ctan the time derivative vector
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dofhandler class
     * @param t_FE the fe space class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyBoundaryConditions(const FECalcType &CalcType,
                                 const double &t,const double (&Ctan)[3],
                                 const FECell &t_FECell,
                                 const DofHandler &t_DofHandler,
                                 FE &t_FE,
                                 Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS);

    /**
     * apply the preset boundary condition, this will only modify the U vector, K and R will not be touched.
     * @param CalcType the calculation type
     * @param t current time
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dofhandler class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyPresetBoundaryConditions(const FECalcType &CalcType,
                                 const double &t,
                                 const FECell &t_FECell,
                                 const DofHandler &t_DofHandler,
                                 Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS);

    /**
     * release the memory
     */
    void releaseMemory();

    /**
     * print out the boundary system info
     */
    void printBCSystemInfo()const;

private:
    /**
     * apply different boundary conditions
     * @param CalcType the calculation type
     * @param t_BCType the boundary condition type
     * @param BCValue the boundary value
     * @param Ctan the time derivative vector
     * @param Params the json content for bc block
     * @param Normal the normal vector of current bc element (on qpoint)
     * @param ElmtInfo the local element info 
     * @param ElmtSoln the local element solution
     * @param Shp the shape function class
     * @param LocalK the local K matrix
     * @param LocalR the local residual vector
     */
    void runBCLibs(const FECalcType &CalcType,
                   const BCType &t_BCType,
                   const double &BCValue,
                   const double (&Ctan)[3],
                   const nlohmann::json &Params,
                   const Vector3d &Normal,
                   const LocalElmtInfo &ElmtInfo,
                   const LocalElmtSolution &ElmtSoln,
                   const LocalShapeFun &Shp,
                   MatrixXd &LocalK,
                   VectorXd &LocalR);

    /**
     * for dirichlet boundary condition
     * @param CalcType the calculation type
     * @param BCValue the boundary value
     * @param t_BCType the boundary condition type
     * @param Params the json content of parameters for bc block
     * @param DofIDs the local dofs id of current bc block
     * @param BCNameList the boundary name list
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dofhandler class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyDirichletBC(const FECalcType &CalcType,
                          const double &BCValue,
                          const BCType &t_BCType,
                          const nlohmann::json &Params,
                          const vector<int> &DofIDs,
                          const vector<string> &BCNameList,
                          const FECell &t_FECell,const DofHandler &t_DofHandler,
                          Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                          Vector &V,
                          SparseMatrix &AMATRIX,
                          Vector &RHS);

    /**
     * for integrated boundary condition
     * @param CalcType the calculation type
     * @param BCValue the boundary value
     * @param t_BCType the boundary condition type
     * @param Ctan the time integration array
     * @param Parameters the json parameter content for bc block
     * @param DofIDs the local dofs id of current bc block
     * @param BCNameList the boundary name list
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dofhandler class
     * @param t_FE the FE space class
     * @param U the solution vector
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyIntegratedBC(const FECalcType &CalcType,
                           const double &BCValue,
                           const BCType &t_BCType,
                           const double (&Ctan)[3],
                           const nlohmann::json &Parameters,
                           const vector<int> &DofIDs,
                           const vector<string> &BCNameList,
                           const FECell &t_FECell,const DofHandler &t_DofHandler,
                           FE &t_FE,
                           Vector &U,Vector &Uold,Vector &Uolder,
                           Vector &V,
                           SparseMatrix &AMATRIX,
                           Vector &RHS);
    
private:
    /**
     * assemble local residual to the global RHS
     * @param Dofs the number of dofs for current bc element
     * @param DofIDs the local dof ids
     * @param I the node index, global one
     * @param JxW the jacobian*weight value
     * @param t_DofHandler the dofhandler class
     * @param LocalR the local residual vector
     * @param RHS the system residual vector
     */
    void assembleLocalResidual2Global(const int &Dofs,const vector<int> &DofIDs,
                                      const int &I,const double &JxW,
                                      const DofHandler &t_DofHandler,
                                      const VectorXd &LocalR,Vector &RHS);
    
    /**
     * assemble local residual to the global RHS
     * @param Dofs the number of dofs for current bc element
     * @param DofIDs the local dof ids
     * @param I the node i-index, global one
     * @param J the node j-index, global one
     * @param JxW the jacobian*weight value
     * @param t_DofHandler the dofhandler class
     * @param LocalK the local K matrix
     * @param AMATRIX the system K matrix
     */
    void assembleLocalJacobian2Global(const int &Dofs,const vector<int> &DofIDs,
                                      const int &I,const int &J,
                                      const double &JxW,
                                      const DofHandler &t_DofHandler,
                                      const MatrixXd &LocalK,
                                      SparseMatrix &AMATRIX);



private:
    vector<BCBlock> m_BCBlockList;/**< vector for bc blocks */
    int m_BCBlocksNum; /**< number of bc blocks */
    double m_DirichletPenalty;/**< the penalty coefficient for dirichlet bc */

private:
    PetscMPIInt m_Rank;/**< for the rank id of current cpu */
    PetscMPIInt m_Size;/**< for the size of total cpus */

private:
    VectorXd m_LocalR;/**< for the local 'element's vector, used for one single element/model */
    MatrixXd m_LocalK;/**< for the local 'element's matrix, used for one single element/model */

    int m_BCElmtNodesNum;/**< the nodes number of the bc element */
    Nodes m_Nodes;/**< for the nodal coordinates of current bc element (current configuration) */
    Nodes m_Nodes0;/**< for the nodal coordinates of current bc element (reference configuration) */

    LocalElmtInfo m_LocalElmtInfo;/**< for the local element information */
    LocalElmtSolution m_LocalElmtSoln;/**< for the local element solution */
    LocalShapeFun m_LocalShp;/**< for the local shape function */

    Vector3d m_Normal;/**< the normal vector */
    Rank2Tensor m_XS;/**< for the derivative in norm vector calculation */

};