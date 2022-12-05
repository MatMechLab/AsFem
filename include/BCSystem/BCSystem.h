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

#include "Mesh/Mesh.h"
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
     * @param t_bcblock the bc block to be added
     */
    void addBCBlock2List(const BCBlock &t_bcblock);
    /**
     * set the dirichlet penalty coefficient
     */
    inline void setDirichletPenalty(const double &val){m_dirichlet_penalty=val;}

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
    inline int getBCBlocksNum()const{return m_bcblocks_num;}
    /**
     * get the i-th bc block
     * @param i the integer index for bc block
     */
    inline BCBlock getIthBCBlock(const int &i)const{
        if(i<1||i>m_bcblocks_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bcblocks");
            MessagePrinter::exitAsFem();
        }
        return m_bclock_list[i-1];
    }
    /**
     * get the penalty coefficient for dirichlet bc
     */
    inline double getDirichletPenalty()const{return m_dirichlet_penalty;}

    //*******************************************************
    //*** for different boundary conditions
    //*******************************************************
    /**
     * apply the boundary conditions
     * @param t_calctype the calculation type
     * @param t current time
     * @param ctan the time derivative vector
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_fe the fe space class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyBoundaryConditions(const FECalcType &t_calctype,
                                 const double &t,const double (&ctan)[3],
                                 const Mesh &t_mesh,
                                 const DofHandler &t_dofhandler,
                                 FE &t_fe,
                                 Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS);

    /**
     * apply the preset boundary condition, this will only modify the U vector, K and R will not be touched.
     * @param t_calctype the calculation type
     * @param t current time
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyPresetBoundaryConditions(const FECalcType &t_calctype,
                                 const double &t,
                                 const Mesh &t_mesh,
                                 const DofHandler &t_dofhandler,
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
     * @param t_calctype the calculation type
     * @param t_bctype the boundary condition type
     * @param t_bcvalue the boundary value
     * @param ctan the time derivative vector
     * @param t_json the json content for bc block
     * @param t_normal the normal vector of current bc element (on qpoint)
     * @param t_elmtinfo the local element info 
     * @param t_elmtsoln the local element solution
     * @param t_shp the shape function class
     * @param localK the local K matrix
     * @param localR the local residual vector
     */
    void runBCLibs(const FECalcType &t_calctype,
                   const BCType &t_bctype,
                   const double &t_bcvalue,
                   const double (&ctan)[3],
                   const nlohmann::json &t_json,
                   const Vector3d &t_normal,
                   const LocalElmtInfo &t_elmtinfo,
                   const LocalElmtSolution &t_elmtsoln,
                   const LocalShapeFun &t_shp,
                   MatrixXd &localK,
                   VectorXd &localR);

    /**
     * for dirichlet boundary condition
     * @param t_calctype the calculation type
     * @param t_bcvalue the boundary value
     * @param t_bctype the boundary condition type
     * @param t_json the json content for bc block
     * @param t_dofids the local dofs id of current bc block
     * @param t_bcnamelist the boundary name list
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param U the solution vector
     * @param Ucopy the solution vector's copy
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyDirichletBC(const FECalcType &t_calctype,
                          const double &t_bcvalue,
                          const BCType &t_bctype,
                          const nlohmann::json &t_json,
                          const vector<int> &t_dofids,
                          const vector<string> &t_bcnamelist,
                          const Mesh &t_mesh,const DofHandler &t_dofhandler,
                          Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                          Vector &V,
                          SparseMatrix &AMATRIX,
                          Vector &RHS);

    /**
     * for integrated boundary condition
     * @param t_calctype the calculation type
     * @param t_bcvalue the boundary value
     * @param t_bctype the boundary condition type
     * @param ctan the time integration array
     * @param t_parameters the json parameter content for bc block
     * @param t_dofids the local dofs id of current bc block
     * @param t_bcnamelist the boundary name list
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_fe the FE space class
     * @param U the solution vector
     * @param Uold the solution vector of previous step
     * @param Uolder the solution vector of pre-previous step
     * @param V the velocity vector
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void applyIntegratedBC(const FECalcType &t_calctype,
                           const double &t_bcvalue,
                           const BCType &t_bctype,
                           const double (&ctan)[3],
                           const nlohmann::json &t_parameters,
                           const vector<int> &t_dofids,
                           const vector<string> &t_bcnamelist,
                           const Mesh &t_mesh,const DofHandler &t_dofhandler,
                           FE &t_fe,
                           Vector &U,Vector &Uold,Vector &Uolder,
                           Vector &V,
                           SparseMatrix &AMATRIX,
                           Vector &RHS);
    
private:
    /**
     * assemble local residual to the global RHS
     * @param t_dofs the number of dofs for current bc element
     * @param t_dofids the local dof ids
     * @param I the node index, global one
     * @param jxw the jacobian*weight value
     * @param t_dofhandler the dofhandler class
     * @param t_localR the local residual vector
     * @param RHS the system residual vector
     */
    void assembleLocalResidual2Global(const int &t_dofs,const vector<int> &t_dofids,
                                      const int &I,const double &jxw,
                                      const DofHandler &t_dofhandler,
                                      const VectorXd &t_localR,Vector &RHS);
    
    /**
     * assemble local residual to the global RHS
     * @param t_dofs the number of dofs for current bc element
     * @param t_dofids the local dof ids
     * @param I the node i-index, global one
     * @param J the node j-index, global one
     * @param jxw the jacobian*weight value
     * @param t_dofhandler the dofhandler class
     * @param t_localK the local K matrix
     * @param AMATRIX the system K matrix
     */
    void assembleLocalJacobian2Global(const int &t_dofs,const vector<int> &t_dofids,
                                      const int &I,const int &J,
                                      const double &jxw,
                                      const DofHandler &t_dofhandler,
                                      const MatrixXd &t_localK,
                                      SparseMatrix &AMATRIX);



private:
    vector<BCBlock> m_bclock_list;/**< vector for bc blocks */
    int m_bcblocks_num; /**< number of bc blocks */
    double m_dirichlet_penalty;/**< the penalty coefficient for dirichlet bc */

private:
    PetscMPIInt m_rank;/**< for the rank id of current cpu */
    PetscMPIInt m_size;/**< for the size of total cpus */

private:
    VectorXd m_localR;/**< for the local 'element's vector, used for one single element/model */
    MatrixXd m_localK;/**< for the local 'element's matrix, used for one single element/model */

    int m_bcelmt_nodesnum;/**< the nodes number of the bc element */
    Nodes m_nodes;/**< for the nodal coordinates of current bc element (current configuration) */
    Nodes m_nodes0;/**< for the nodal coordinates of current bc element (reference configuration) */

    LocalElmtInfo m_local_elmtinfo;/**< for the local element information */
    LocalElmtSolution m_local_elmtsoln;/**< for the local element solution */
    LocalShapeFun m_local_shp;/**< for the local shape function */

    Vector3d m_normal;/**< the normal vector */
    Rank2Tensor m_xs;/**< for the derivative in norm vector calculation */

};