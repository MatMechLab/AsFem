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
//+++ Date   : 2025.01.24
//+++ Purpose: Implement full least square projection from integration
//+++          points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "ProjectionSystem/ProjectionBase.h"

/**
 * This class implements the least square projection from integration points to nodal points
 */
class FullLeastSquareProjection:public ProjectionBase{
public:
    FullLeastSquareProjection();

protected:
    /**
     * initialize the projection system
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof handler class
     */
    void initMyProjection(const FECell &t_FECell,const DofHandler &t_DofHandler) final;
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
    void executeMyProjection(const FECell &t_FECell,
                             const DofHandler &t_DofHandler,
                             const ElmtSystem &t_ElmtSystem,
                             MateSystem &t_MateSystem,
                             FE &t_FE,
                             SolutionSystem &t_SolnSystem,
                             const FEControlInfo &t_FECtrlInfo,
                             ProjectionData &Data) final;
private:
    /**
     * the global projection action
     * @param t_FECell the fe cell class
     * @param t_Data the projection data structure
     */
    void globalProjectionAction(const FECell &t_FECell,ProjectionData &Data);
    /**
     * the local projection action
     * @param NodesNum the nodes number of current element
     * @param ElmtID the element id
     * @param QPointsNum the number of qpoints
     * @param ALeft the extrapolate matrix to be solved
     * @param ARight the extrapolate matrix for the right hand side
     * @param ElConn the local element connectivity
     * @param Volume the volume of the current element
     * @param SolnSystem the solution system class
     * @param Data the projection data structure
     */
    void localProjectionAction(const int &NodesNum,
                               const int &ElmtID,
                               const int &QPointsNum,
                               Eigen::MatrixXd &ALeft,
                               Eigen::MatrixXd &ARight,
                               const vector<int> &ElConn,
                               const double &Volume,
                               const SolutionSystem &SolnSystem,
                               ProjectionData &Data);


    /**
     * the local projection action for scalar material
     * @param NodesNum the nodes number of current element
     * @param ElmtID the element id
     * @param QPointsNum the number of qpoints
     * @param ALeft the extrapolate matrix to be solved
     * @param ARight the extrapolate matrix for the right hand side
     * @param ElConn the elemental connectivity, start from 1
     * @param Volume the jacobian determinte of current qpoint
     * @param SolnSystem the solution system class
     * @param ProjNum the number of scalar material to be projected
     * @param ScalarMateNameList the scalar material name list
     * @param ScalarProjMateVecList the list of scalar projection vector
     */
    void projectLocalScalarMate2Global(const int &NodesNum,
                                       const int &ElmtID,
                                       const int &QPointsNum,
                                       Eigen::MatrixXd &ALeft,
                                       Eigen::MatrixXd &ARight,
                                       const vector<int> &ElConn,
                                       const double &Volume,
                                       const SolutionSystem &SolnSystem,
                                       const int &ProjNum,
                                       const vector<string> &ScalarMateNameList,
                                       vector<Vector> &ScalarProjMateVecList);

    /**
     * the local projection action for vector material
     * @param NodesNum the nodes number of current element
     * @param ElmtID the element id
     * @param QPointsNum the number of qpoints
     * @param ALeft the extrapolate matrix to be solved
     * @param ARight the extrapolate matrix for the right hand side
     * @param ElConn the elemental connectivity, start from 1
     * @param Volume the jacobian determinte of current qpoint
     * @param SolnSystem the solution system class
     * @param ProjNum the number of vector material to be projected
     * @param VectorMateNameList the vector material name list
     * @param VectorProjMateVecList the vector projection vector
     */
    void projectLocalVectorMate2Global(const int &NodesNum,
                                       const int &ElmtID,
                                       const int &QPointsNum,
                                       Eigen::MatrixXd &ALeft,
                                       Eigen::MatrixXd &ARight,
                                       const vector<int> &ElConn,
                                       const double &Volume,
                                       const SolutionSystem &SolnSystem,
                                       const int &ProjNum,
                                       const vector<string> &VectorMateNameList,
                                       vector<Vector> &VectorProjMateVecList);
    /**
     * the local projection action for rank-2 material
     * @param NodesNum the nodes number of current element
     * @param ElmtID the element id
     * @param QPointsNum the number of qpoints
     * @param ALeft the extrapolate matrix to be solved
     * @param ARight the extrapolate matrix for the right hand side
     * @param ElConn the elemental connectivity, start from 1
     * @param Volume the jacobian determinte of current qpoint
     * @param SolnSystem the solution system class
     * @param ProjNum the number of vector material to be projected
     * @param Rank2MateNameList the rank-2 material name list
     * @param Rank2ProjMateVecList the rank-2 projection vector
     */
    void projectLocalRank2Mate2Global(const int &NodesNum,
                                      const int &ElmtID,
                                      const int &QPointsNum,
                                      Eigen::MatrixXd &ALeft,
                                      Eigen::MatrixXd &ARight,
                                      const vector<int> &ElConn,
                                      const double &Volume,
                                      const SolutionSystem &SolnSystem,
                                      const int &ProjNum,
                                      const vector<string> &Rank2MateNameList,
                                      vector<Vector> &Rank2ProjMateVecList);

    /**
     * the local projection action for rank-4 material
     * @param NodesNum the nodes number of current element
     * @param ElmtID the element id
     * @param QPointsNum the number of qpoints
     * @param ALeft the extrapolate matrix to be solved
     * @param ARight the extrapolate matrix for the right hand side
     * @param ElConn the elemental connectivity, start from 1
     * @param Volume the jacobian determinte of current qpoint
     * @param SolnSystem the solution system class
     * @param ProjNum the number of vector material to be projected
     * @param Rank4MateNameList the rank-4 material name list
     * @param Rank4ProjMateVecList the rank-4 projection vector
     */
    void projectLocalRank4Mate2Global(const int &NodesNum,
                                      const int &ElmtID,
                                      const int &QPointsNum,
                                      Eigen::MatrixXd &ALeft,
                                      Eigen::MatrixXd &ARight,
                                      const vector<int> &ElConn,
                                      const double &Volume,
                                      const SolutionSystem &SolnSystem,
                                      const int &ProjNum,
                                      const vector<string> &Rank4MateNameList,
                                      vector<Vector> &Rank4ProjMateVecList);


private:
    int m_BulkElmtNodesNum;/**< the nodes number of the bulk element */
    vector<int> m_ElmtConn;/**< for local element's connectivity */
    int         m_SubElmtDofs;/**< for the dofs number of each sub element */
    vector<int> m_SubElmtDofIDs;/**< for local sub-elemental nodes' gloabl ids, start from 0 */

    Nodes m_Nodes;/**< for the nodal coordinates of current bulk element (current configuration) */
    Nodes m_Nodes0;/**< for the nodal coordinates of current bulk element (reference configuration) */

    LocalElmtInfo m_LocalElmtInfo;/**< for the local element information */
    LocalElmtSolution m_LocalElmtSoln;/**< for the local element solution */

    Eigen::MatrixXd A,AT,Al,Ar,W;/**< the extrapolate matrix */
    Eigen::VectorXd QpSolnVec,NodalSolnVec,RHSVec;/**< the right hand side vector and the solution vector */

    bool IsFirstSetup;


};