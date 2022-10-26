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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "petsc.h"

/**
 * For AsFem's own header files
 */
#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "EquationSystem/EquationSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/Vector.h"

#include "MathUtils/MatrixXd.h"
#include "MathUtils/SparseMatrix.h"

#include "FESystem/FECalcType.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"

#include "ElmtSystem/LocalElmtData.h"


/**
 * This class implements the general FEM calculation actions, 
 * i.e. local K and R calculation and assemble..., for AsFem.
 */
class BulkFESystem{
public:
    /**
     * constructor
     */
    BulkFESystem();

    //***************************************************************
    //*** general settings
    //***************************************************************
    /**
     * init the bulk FE system
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof handler class
     */
    void init(const Mesh &t_mesh, const DofHandler &t_dofhandler);
    
    /**
     * reset the maximum coefficient of K matrix
     */
    inline void resetMaxKMatrixCoeff(){m_max_k_coeff=-1.0e16;}

    /**
     * generate the system residual and jacobian based on different elements(PDEs/ODEs)
     * @param t_calctype the calculation action type, i.e., form residual or form jacobian
     * @param t double value for the current time
     * @param dt double value for the current time increment
     * @param ctan the 3-d vector for the time derivative coefficients
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofHandler class
     * @param t_fe the fe class for the shape function and gauss points
     * @param t_elmtsystem the element system class
     * @param t_matesystem the material system class
     * @param t_solutionsystem the solution system class for 'U', 'V', and 'A'
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void formBulkFE(const FECalcType &t_calctype,const double &t,const double &dt,const double (&ctan)[3],
                    Mesh &t_mesh,const DofHandler &t_dofhandler,FE &t_fe,
                    ElmtSystem &t_elmtsystem,MateSystem &t_matesystem,
                    SolutionSystem &t_solutionsystem,
                    SparseMatrix &AMATRIX,Vector &RHS);
    /**
     * release the allocated memory
     */
    void releaseMemory();

    //***************************************************************
    //*** general gettings
    //***************************************************************
    /**
     * get the maximum coefficient of current K matrix
     */
    inline double getMaxCoefOfKMatrix()const{return m_max_k_coeff;}

private:
    //***************************************************************
    //*** assemble functions
    //***************************************************************
    /**
     * assemble local residual to global residual
     * @param t_dofs the dofs number of current sub element
     * @param t_dofsid the dofs id of current sub element, the local one, not the global ids!
     * @param t_globalnodeid the global node index
     * @param t_dofhandler the dofHandler class
     * @param jxw the JxW for integration
     * @param t_subR the residual of current sub element
     * @param RHS the global residual
     */
    void assembleLocalResidual2GlobalR(const int &t_dofs,const vector<int> &t_dofsid,
                                       const int &t_globalnodeid,
                                       const DofHandler &t_dofhandler,
                                       const double &jxw,
                                       const VectorXd &t_subR,
                                       Vector &RHS);
    /**
     * assemble local residual to global residual
     * @param t_dofs the dofs number of current sub element
     * @param t_dofsid the dofs id of current sub element, the local one
     * @param t_globalnodeidI the node index (I-index), global one
     * @param t_globalnodeidJ the node index (J-index), global one
     * @param jxw JxW for integration
     * @param t_dofhandler the dofHandler class
     * @param t_subK the jacobian of current sub element
     * @param AMATRIX the global K matrix
     */
    void assembleLocalJacobian2GlobalK(const int &t_dofs,const vector<int> &t_dofsid,
                                       const int &t_globalnodeidI,const int &t_globalnodeidJ,
                                       const double &jxw,
                                       const DofHandler &t_dofhandler,
                                       const MatrixXd &t_subK,
                                       SparseMatrix &AMATRIX);


private:
    int m_max_nodal_dofs;/**< the maximum dofs number of each node */
    int m_max_elmt_dofs;/**< the maximum dofs number of each bulk element */
    int m_bulkelmt_nodesnum;/**< the nodes number of the bulk element */
    VectorXd m_localR;/**< for the local Residual vector, this is used for the whole element */
    MatrixXd m_localK;/**< for the local Jacobian matrix, this is used for the whole element */
    VectorXd m_subR;/**< for the local 'element's vector, used for one single element/model */
    MatrixXd m_subK;/**< for the local 'element's matrix, used for one single element/model */

    vector<int> m_elmtconn;/**< for local element's connectivity */
    vector<int> m_elmtdofsid;/**< for local elemental nodes' gloabl ids, start from 0 */
    int         m_subelmt_dofs;/**< for the dofs number of each sub element */
    vector<int> m_subelmtdofsid;/**< for local sub-elemental nodes' gloabl ids, start from 0 */

    double m_max_k_coeff;/**< the max(absolute) value of current K matrix */

    vector<double> m_elmtU;/**< 'displacemen' of current element */
    vector<double> m_elmtUold;/**< previous 'displacemen' of current element */
    vector<double> m_elmtUolder;/**< pre-previous 'displacemen' of current element */
    vector<double> m_elmtV;/**< 'velocity' of current element */
    vector<double> m_elmtA;/**< 'acceleration' of current element */

    Nodes m_nodes;/**< for the nodal coordinates of current bulk element (current configuration) */
    Nodes m_nodes0;/**< for the nodal coordinates of current bulk element (reference configuration) */

    LocalElmtInfo m_local_elmtinfo;/**< for the local element information */
    LocalElmtSolution m_local_elmtsoln;/**< for the local element solution */
    LocalShapeFun m_local_shp;/**< for the local shape function */

private:
    PetscMPIInt m_rank;/**< for the rank id of current cpu */
    PetscMPIInt m_size;/**< for the size of total cpus */

    
};