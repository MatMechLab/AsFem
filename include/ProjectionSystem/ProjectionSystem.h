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
//+++ Purpose: Implement the general projection system, which can do
//+++          the extropolation from guass points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector.h"

#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "FE/FE.h"
#include "SolutionSystem/SolutionSystem.h"
#include "FEProblem/FEControlInfo.h"

#include "ProjectionSystem/ProjectionType.h"
#include "ProjectionSystem/ProjectionData.h"
#include "ElmtSystem/LocalElmtData.h"

#include "ProjectionSystem/LeastSquareProjection.h"
#include "ProjectionSystem/FullLeastSquareProjection.h"

/**
 * This class implements the general projections method for extrapolating the gauss point's quantities
 */
class ProjectionSystem:public LeastSquareProjection,
                       public FullLeastSquareProjection{
public:
    /**
     * constructor
     */
    ProjectionSystem();

    /**
     * initialize the projection system
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof handler class
     */
    void init(const FECell &t_FECell,const DofHandler &t_DofHandler);

    //********************************************************
    //*** general settings
    //********************************************************
    /**
     * set the projection status
     * @param Flag boolean flag for projection status
     */
    void setProjectionStatus(const bool &Flag){m_IsProjection=Flag;}
    /**
     * set the number of scalar materials to be projected
     * @param num integer
     */
    void setScalarMaterialNum(const int &num){m_Data.m_ScalarProjMateNum=num;}
    /**
     * set the number of vector materials to be projected
     * @param num integer
     */
    void setVectorMaterialNum(const int &num){m_Data.m_VectorProjMateNum=num;}
    /**
     * set the number of rank-2 tensor materials to be projected
     * @param num integer
     */
    void setRank2TensorMaterialNum(const int &num){m_Data.m_Rank2ProjMateNum=num;}
    /**
     * set the number of rank-4 tensor materials to be projected
     * @param num integer
     */
    void setRank4TensorMaterialNum(const int &num){m_Data.m_Rank4ProjMateNum=num;}
    /**
     * set the projection type from input file
     * @param projtype the projection type
     */
    void setProjectionType(const ProjectionType &projtype){m_ProjType=projtype;}
    /**
     * add scalar material name to list
     * @param matename string name for scalar material
     */
    void addScalarMateName2List(const string &matename);
    /**
     * add vector material name to list
     * @param matename string name for scalar material
     */
    void addVectorMateName2List(const string &matename);
    /**
     * add rank-2 tensor material name to list
     * @param matename string name for scalar material
     */
    void addRank2MateName2List(const string &matename);
    /**
     * add rank-4 material name to list
     * @param matename string name for scalar material
     */
    void addRank4MateName2List(const string &matename);

    //********************************************************
    //*** general gettings
    //********************************************************
    /**
     * get the active status of projection
     */
    inline bool isProjectionActive()const{return m_IsProjection;}
    /**
     * get the nodes number
     */
    inline int getNodesNum()const{return m_NodesNum;}

    /**
     * get the number of scalar materials
     */
    inline int getScalarMaterialNum()const{return m_Data.m_ScalarProjMateNum;}
    /**
     * get the i-th scalar material name
     * @param i i index for scalar material name list
     */
    inline string getIthScalarMateName(const int &i)const{
        if(i<1||i>m_Data.m_ScalarProjMateNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_Data.m_ScalarProjMateNum)+") for scalar material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_Data.m_ScalarProjMateNameList[i-1];
    }

    /**
     * get the number of vector materials
     */
    inline int getVectorMaterialNum()const{return m_Data.m_VectorProjMateNum;}
    /**
     * get the i-th vector material name
     * @param i i index for vector material name list
     */
    inline string getIthVectorMateName(const int &i)const{
        if(i<1||i>m_Data.m_VectorProjMateNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_Data.m_VectorProjMateNum)+") for vector material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_Data.m_VectorProjMateNamelist[i-1];
    }

    /**
     * get the number of rank-2 tensor materials
     */
    inline int getRank2MaterialNum()const{return m_Data.m_Rank2ProjMateNum;}
    /**
     * get the i-th rank-2 material name
     * @param i i index for rank-2 material name list
     */
    inline string getIthRank2MateName(const int &i)const{
        if(i<1||i>m_Data.m_Rank2ProjMateNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_Data.m_Rank2ProjMateNum)+") for rank-2 material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_Data.m_Rank2ProjMateNameList[i-1];
    }

    /**
     * get the number of rank-4 materials
     */
    inline int getRank4MaterialNum()const{return m_Data.m_Rank4ProjMateNum;}
    /**
     * get the i-th rank-4 material name
     * @param i i index for rank-4 material name list
     */
    inline string getIthRank4MateName(const int &i)const{
        if(i<1||i>m_Data.m_Rank4ProjMateNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_Data.m_Rank4ProjMateNum)+") for rank-2 material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_Data.m_Rank4ProjMateNameList[i-1];
    }

    /**
     * get the reference of projection data
     */
    inline ProjectionData& getProjectionDataRef(){return m_Data;}
    /**
     * get the copy of projection data
     */
    inline ProjectionData getProjectionDataCopy()const{return m_Data;}

    //*********************************************
    //*** for general gettings
    //*********************************************
    /**
     * get the i-th node's scalar material value via its scalar material name
     * @param nodeid the global node id, starts from 0
     * @param matename the scalar material name
     */
    double getIthNodeScalarMateViaMateName(const int &nodeid,const string &matename);

    /**
     * get the i-th node's vector material value via its vector material name
     * @param nodeid the global node id, starts from 0
     * @param matename the vector material name
     */
    Vector3d getIthNodeVectorMateViaMateName(const int &nodeid,const string &matename);

    /**
     * get the i-th node's rank-2 tensor material value via its rank-2 tensor material name
     * @param nodeid the global node id, starts from 0
     * @param matename the rank-2 material name
     */
    Rank2Tensor getIthNodeRank2MateViaMateName(const int &nodeid,const string &matename);

    /**
     * get the i-th node's rank-4 tensor material value via its rank-4 tensor material name
     * @param nodeid the global node id, starts from 0
     * @param matename the rank-4 material name
     */
    Rank4Tensor getIthNodeRank4MateViaMateName(const int &nodeid,const string &matename);

    /**
     * execute the projection process
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dofhandler class
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_FE the fe space class
     * @param t_SolnSystem the solution system class
     * @param t_FECtrlInfo the fe control info
     */
    void executeProjection(const FECell &t_FECell,
                           const DofHandler &t_DofHandler,
                           const ElmtSystem &t_ElmtSystem,
                           MateSystem &t_MateSystem,
                           FE &t_FE,
                           SolutionSystem &t_SolnSystem,
                           const FEControlInfo &t_FECtrlInfo);

    /**
     * make the ghost copy for all the projection data vector
    */
    void makeGhostCopyOfProjectionData();
    /**
     * destroy the ghost copy for all the projection data vector
    */
    void destroyGhostCopyOfProjectionData();

    /**
     * release the allocated memory
     */
    void releaseMemory();

    /**
     * print out the projection information
     */
    void printProjectionInfo()const;
    

private:
    bool m_IsAllocated;/**< boolean flag for memory allocation status */
    bool m_IsProjection;/**< boolean flag for the projection status */
    int m_NodesNum;/** total nodes of bulk mesh */

    ProjectionType m_ProjType;/**< the type of projection method */
    ProjectionData m_Data;/**< the projection data */

    int m_BulkElmtNodesNum;/**< the nodes number of the bulk element */
    vector<int> m_ElmtConn;/**< for local element's connectivity */
    int         m_SubElmtDofs;/**< for the dofs number of each sub element */
    vector<int> m_SubElmtDofIDs;/**< for local sub-elemental nodes' gloabl ids, start from 0 */

    Nodes m_Nodes;/**< for the nodal coordinates of current bulk element (current configuration) */
    Nodes m_Nodes0;/**< for the nodal coordinates of current bulk element (reference configuration) */

    LocalElmtInfo m_LocalElmtInfo;/**< for the local element information */
    LocalElmtSolution m_LocalElmtSoln;/**< for the local element solution */

    PetscMPIInt m_Rank;/**< for the rank id of current cpu */
    PetscMPIInt m_Size;/**< for the size of total cpus */

};