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
#include "FECell/Nodes.h"
#include "FECell/FECell.h"
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
     * @param t_fecell the fe cell class
     * @param t_dofhandler the dof handler class
     */
    void init(const FECell &t_fecell, const DofHandler &t_dofhandler);
    
    /**
     * reset the maximum coefficient of K matrix
     */
    inline void resetMaxKMatrixCoeff(){m_MaxKmatCoeff=-1.0e16;}

    /**
     * generate the system residual and jacobian based on different elements(PDEs/ODEs)
     * @param t_CalcType the calculation action type, i.e., form residual or form jacobian
     * @param T double value for the current time
     * @param Dt double value for the current time increment
     * @param Ctan the 3-d vector for the time derivative coefficients
     * @param t_FEcell the fe cell class
     * @param t_DofHandler the dofHandler class
     * @param t_FE the fe class for the shape function and gauss points
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_SolnSystem the solution system class for 'U', 'V', and 'A'
     * @param AMATRIX the system K matrix
     * @param RHS the system residual vector
     */
    void formBulkFE(const FECalcType &t_CalcType,
                    const double &T,
                    const double &Dt,
                    const double (&Ctan)[3],
                    FECell &t_FEcell,
                    const DofHandler &t_DofHandler,
                    FE &t_FE,
                    ElmtSystem &t_ElmtSystem,
                    MateSystem &t_MateSystem,
                    SolutionSystem &t_SolnSystem,
                    SparseMatrix &AMATRIX,
                    Vector &RHS);
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
    inline double getMaxCoefOfKMatrix()const{return m_MaxKmatCoeff;}

private:
    //***************************************************************
    //*** assemble functions
    //***************************************************************
    /**
     * assemble local residual to global residual
     * @param Dofs the dofs number of current sub element
     * @param DofIDs the dofs id of current sub element, the local one, not the global ids!
     * @param GlobalNodeID the global node index
     * @param t_DofHandler the dofHandler class
     * @param JxW the JxW for integration
     * @param SubR the residual of current sub element
     * @param RHS the global residual
     */
    void assembleLocalResidual2GlobalR(const int &Dofs,
                                       const vector<int> &DofIDs,
                                       const int &GlobalNodeID,
                                       const DofHandler &t_DofHandler,
                                       const double &JxW,
                                       const VectorXd &SubR,
                                       Vector &RHS);
    /**
     * assemble local residual to global residual
     * @param Dofs the dofs number of current sub element
     * @param DofIDs the dofs id of current sub element, the local one
     * @param GlobalNodeIDI the node index (I-index), global one
     * @param GlobalNodeIDJ the node index (J-index), global one
     * @param jxw JxW for integration
     * @param t_dofhandler the dofHandler class
     * @param SubK the jacobian of current sub element
     * @param AMATRIX the global K matrix
     */
    void assembleLocalJacobian2GlobalK(const int &Dofs,
                                       const vector<int> &DofIDs,
                                       const int &GlobalNodeIDI,
                                       const int &GlobalNodeIDJ,
                                       const double &JxW,
                                       const DofHandler &t_DofHandler,
                                       const MatrixXd &SubK,
                                       SparseMatrix &AMATRIX);


private:
    int m_MaxNodalDofs;/**< the maximum dofs number of each node */
    int m_MaxElmtDofs;/**< the maximum dofs number of each bulk element */
    int m_BulkElmtNodesNum;/**< the nodes number of the bulk element */
    VectorXd m_LocalR;/**< for the local Residual vector, this is used for the whole element */
    MatrixXd m_LocalK;/**< for the local Jacobian matrix, this is used for the whole element */
    VectorXd m_SubR;/**< for the local 'element's vector, used for one single element/model */
    MatrixXd m_SubK;/**< for the local 'element's matrix, used for one single element/model */

    vector<int> m_ElmtConn;/**< for local element's connectivity */
    vector<int> m_ElmtDofIDs;/**< for local elemental nodes' gloabl ids, start from 0 */
    int         m_SubElmtDofs;/**< for the dofs number of each sub element */
    vector<int> m_SubElmtDofIDs;/**< for local sub-elemental nodes' gloabl ids, start from 0 */

    double m_MaxKmatCoeff;/**< the max(absolute) value of current K matrix */

    vector<double> m_ElmtU;/**< 'displacemen' of current element */
    vector<double> m_ElmtUold;/**< previous 'displacemen' of current element */
    vector<double> m_ElmtUolder;/**< pre-previous 'displacemen' of current element */
    vector<double> m_ElmtV;/**< 'velocity' of current element */
    vector<double> m_ElmtA;/**< 'acceleration' of current element */

    Nodes m_Nodes;/**< for the nodal coordinates of current bulk element (current configuration) */
    Nodes m_Nodes0;/**< for the nodal coordinates of current bulk element (reference configuration) */

    LocalElmtInfo m_LocalElmtInfo;/**< for the local element information */
    LocalElmtSolution m_LocalElmtSoln;/**< for the local element solution */
    LocalShapeFun m_LocalShp;/**< for the local shape function */

private:

    
};