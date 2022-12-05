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
//+++ Date   : 2022.07.22
//+++ Purpose: Implement the general projection system, which can do
//+++          the extropolation from guass points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector.h"

#include "Mesh/Mesh.h"
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

/**
 * This class implements the general projections method for extrapolating the gauss point's quantities
 */
class ProjectionSystem:public LeastSquareProjection{
public:
    /**
     * constructor
     */
    ProjectionSystem();

    /**
     * initialize the projection system
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof handler class
     */
    void init(const Mesh &t_mesh,const DofHandler &t_dofhandler);

    //********************************************************
    //*** general settings
    //********************************************************
    /**
     * set the projection status
     * @param flag boolean flag for projection status
     */
    void setProjectionStatus(const bool &flag){m_isprojection=flag;}
    /**
     * set the number of scalar materials to be projected
     * @param num integer
     */
    void setScalarMaterialNum(const int &num){m_data.m_scalarmate_num=num;}
    /**
     * set the number of vector materials to be projected
     * @param num integer
     */
    void setVectorMaterialNum(const int &num){m_data.m_vectormate_num=num;}
    /**
     * set the number of rank-2 tensor materials to be projected
     * @param num integer
     */
    void setRank2TensorMaterialNum(const int &num){m_data.m_rank2mate_num=num;}
    /**
     * set the number of rank-4 tensor materials to be projected
     * @param num integer
     */
    void setRank4TensorMaterialNum(const int &num){m_data.m_rank4mate_num=num;}
    /**
     * set the projection type from input file
     * @param projtype the projection type
     */
    void setProjectionType(const ProjectionType &projtype){m_proj_type=projtype;}
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
    inline bool isProjectionActive()const{return m_isprojection;}
    /**
     * get the nodes number
     */
    inline int getNodesNum()const{return m_nodesnum;}

    /**
     * get the number of scalar materials
     */
    inline int getScalarMaterialNum()const{return m_data.m_scalarmate_num;}
    /**
     * get the i-th scalar material name
     * @param i i index for scalar material name list
     */
    inline string getIthScalarMateName(const int &i)const{
        if(i<1||i>m_data.m_scalarmate_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_data.m_scalarmate_num)+") for scalar material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_data.m_scalarmate_namelist[i-1];
    }

    /**
     * get the number of vector materials
     */
    inline int getVectorMaterialNum()const{return m_data.m_vectormate_num;}
    /**
     * get the i-th vector material name
     * @param i i index for vector material name list
     */
    inline string getIthVectorMateName(const int &i)const{
        if(i<1||i>m_data.m_vectormate_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_data.m_vectormate_num)+") for vector material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_data.m_vectormate_namelist[i-1];
    }

    /**
     * get the number of rank-2 tensor materials
     */
    inline int getRank2MaterialNum()const{return m_data.m_rank2mate_num;}
    /**
     * get the i-th rank-2 material name
     * @param i i index for rank-2 material name list
     */
    inline string getIthRank2MateName(const int &i)const{
        if(i<1||i>m_data.m_rank2mate_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_data.m_rank2mate_num)+") for rank-2 material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_data.m_rank2mate_namelist[i-1];
    }

    /**
     * get the number of rank-4 materials
     */
    inline int getRank4MaterialNum()const{return m_data.m_rank4mate_num;}
    /**
     * get the i-th rank-4 material name
     * @param i i index for rank-4 material name list
     */
    inline string getIthRank4MateName(const int &i)const{
        if(i<1||i>m_data.m_rank4mate_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range ("+to_string(m_data.m_rank4mate_num)+") for rank-2 material,"
                                        +"please check either your code or your input file");
            MessagePrinter::exitAsFem();
        }
        return m_data.m_rank4mate_namelist[i-1];
    }

    /**
     * get the reference of projection data
     */
    inline ProjectionData& getProjectionDataRef(){return m_data;}
    /**
     * get the copy of projection data
     */
    inline ProjectionData getProjectionDataCopy()const{return m_data;}

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
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_elmtsystem the element system class
     * @param t_matesystem the material system class
     * @param t_fe the fe space class
     * @param t_solution the solution system class
     * @param t_fectrlinfo the fe control info 
     */
    void executeProjection(const Mesh &t_mesh,
                           const DofHandler &t_dofhandler,
                           const ElmtSystem &t_elmtsystem,
                           MateSystem &t_matesystem,
                           FE &t_fe,
                           SolutionSystem &t_solution,
                           const FEControlInfo &t_fectrlinfo);

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
    /**
     * for different projection methods
     * @param flag true for local, false for global
     * @param t_mesh the mesh class
     * @param nodesnum nodes number of current element
     * @param t_elconn the local element's connectivity
     * @param detjac the jacobian determinte
     * @param t_shp the shape function
     * @param t_mate the material container
     * @param t_data the projection data
     */
    void runProjectionLibs(const bool &flag,
                           const Mesh &t_mesh,
                           const int &nodesnum,
                           const vector<int> &t_elconn,
                           const double &detjac,
                           const ShapeFun &t_shp,
                           const MaterialsContainer &t_mate,
                           ProjectionData &t_data);
    

private:
    bool m_isallocated;/**< boolean flag for memory allocation status */
    bool m_isprojection;/**< boolean flag for the projection status */
    int m_nodesnum;/** total nodes of bulk mesh */

    ProjectionType m_proj_type;/**< the type of projection method */
    ProjectionData m_data;/**< the projection data */

    int m_bulkelmt_nodesnum;/**< the nodes number of the bulk element */
    vector<int> m_elmtconn;/**< for local element's connectivity */
    int         m_subelmt_dofs;/**< for the dofs number of each sub element */
    vector<int> m_subelmtdofsid;/**< for local sub-elemental nodes' gloabl ids, start from 0 */

    Nodes m_nodes;/**< for the nodal coordinates of current bulk element (current configuration) */
    Nodes m_nodes0;/**< for the nodal coordinates of current bulk element (reference configuration) */

    LocalElmtInfo m_local_elmtinfo;/**< for the local element information */
    LocalElmtSolution m_local_elmtsoln;/**< for the local element solution */

    PetscMPIInt m_rank;/**< for the rank id of current cpu */
    PetscMPIInt m_size;/**< for the size of total cpus */

};