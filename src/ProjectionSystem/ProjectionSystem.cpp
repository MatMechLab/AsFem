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

#include "ProjectionSystem/ProjectionSystem.h"

ProjectionSystem::ProjectionSystem(){
    m_IsAllocated=false;/**< boolean flag for memory allocation status */
    m_IsProjection=false;
    m_NodesNum=0;/** total nodes of bulk mesh */

    m_Data.m_scalarmate_num=0;/**< number of scalar material to be projected */
    m_Data.m_vectormate_num=0;/**< number of vector material to be projected */
    m_Data.m_rank2mate_num=0;/**< number of rank2tensor material to be projected */
    m_Data.m_rank4mate_num=0;/**< number of rank4tensor material to be projected */

    m_Data.m_scalarmate_namelist.clear();/**< vector for the name of scalar materials to be projected */
    m_Data.m_vectormate_namelist.clear();/**< vector for the name of vector materials to be projected */
    m_Data.m_rank2mate_namelist.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_Data.m_rank4mate_namelist.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_ProjType=ProjectionType::DEFAULT;/**< the type of projection method */
    m_Data.m_proj_scalarmate_vec.releaseMemory();/**< the vector stores the projected scalar materials */
    m_Data.m_proj_vectormate_vec.releaseMemory();/**< the vector stores the projected vector materials */
    m_Data.m_proj_rank2mate_vec.releaseMemory();/**< the vector stores the projected rank-2 materials */
    m_Data.m_proj_rank4mate_vec.releaseMemory();/**< the vector stores the projected rank-4 materials */
}
//******************************************************************
//*** init
//******************************************************************
void ProjectionSystem::init(const FECell &t_fecell,const DofHandler &t_dofhandler){
    m_NodesNum=t_fecell.getFECellNodesNum();
    if(getScalarMaterialNum()){
        m_Data.m_proj_scalarmate_vec.resize(m_NodesNum*(1+1*m_Data.m_scalarmate_num),0.0);
    }
    if(getVectorMaterialNum()){
        m_Data.m_proj_vectormate_vec.resize(m_NodesNum*(1+3*m_Data.m_vectormate_num),0.0);
    }
    if(getRank2MaterialNum()){
        m_Data.m_proj_rank2mate_vec.resize(m_NodesNum*(1+9*m_Data.m_rank2mate_num),0.0);
    }
    if(getRank4MaterialNum()){
        m_Data.m_proj_rank4mate_vec.resize(m_NodesNum*(1+36*m_Data.m_rank4mate_num),0.0);
    }



    int m_max_nodal_dofs=t_dofhandler.getMaxDofsPerNode();
    m_BulkElmtNodesNum=t_fecell.getFECellNodesNumPerBulkElmt();

    m_ElmtConn.resize(m_BulkElmtNodesNum,0);
    m_SubElmtDofIDs.resize(m_max_nodal_dofs+1,0);

    m_LocalElmtSoln.m_QpU.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpUold.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpUolder.resize(m_max_nodal_dofs+1,0.0);

    m_LocalElmtSoln.m_QpV.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpA.resize(m_max_nodal_dofs+1,0.0);

    m_LocalElmtSoln.m_QpGradU.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUold.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUolder.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_QpGradV.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_Qpgradu.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_Qpgradv.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtInfo.m_QpCoords=0.0;
    m_LocalElmtInfo.m_QpCoords0=0.0;

    m_Nodes.resize(m_BulkElmtNodesNum);
    m_Nodes0.resize(m_BulkElmtNodesNum);
}

void ProjectionSystem::makeGhostCopyOfProjectionData(){
    m_Data.m_proj_scalarmate_vec.makeGhostCopy();
    m_Data.m_proj_vectormate_vec.makeGhostCopy();
    m_Data.m_proj_rank2mate_vec.makeGhostCopy();
    m_Data.m_proj_rank4mate_vec.makeGhostCopy();
}
void ProjectionSystem::destroyGhostCopyOfProjectionData(){
    m_Data.m_proj_scalarmate_vec.destroyGhostCopy();
    m_Data.m_proj_vectormate_vec.destroyGhostCopy();
    m_Data.m_proj_rank2mate_vec.destroyGhostCopy();
    m_Data.m_proj_rank4mate_vec.destroyGhostCopy();
}
//******************************************************************
//*** release memory
//******************************************************************
void ProjectionSystem::releaseMemory(){
    m_IsAllocated=false;/**< boolean flag for memory allocation status */
    m_IsProjection=false;

    m_NodesNum=0;/** total nodes of bulk mesh */

    m_Data.m_scalarmate_num=0;/**< number of scalar material to be projected */
    m_Data.m_vectormate_num=0;/**< number of vector material to be projected */
    m_Data.m_rank2mate_num=0;/**< number of rank2tensor material to be projected */
    m_Data.m_rank4mate_num=0;/**< number of rank4tensor material to be projected */

    m_Data.m_scalarmate_namelist.clear();/**< vector for the name of scalar materials to be projected */
    m_Data.m_vectormate_namelist.clear();/**< vector for the name of vector materials to be projected */
    m_Data.m_rank2mate_namelist.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_Data.m_rank4mate_namelist.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_ProjType=ProjectionType::DEFAULT;/**< the type of projection method */
    m_Data.m_scalarmate_num=0;
    m_Data.m_vectormate_num=0;
    m_Data.m_rank2mate_num=0;
    m_Data.m_rank4mate_num=0;
    m_Data.m_proj_scalarmate_vec.releaseMemory();/**< the vector stores the projected scalar materials */
    m_Data.m_proj_vectormate_vec.releaseMemory();/**< the vector stores the projected vector materials */
    m_Data.m_proj_rank2mate_vec.releaseMemory();/**< the vector stores the projected rank-2 materials */
    m_Data.m_proj_rank4mate_vec.releaseMemory();/**< the vector stores the projected rank-4 materials */

    m_LocalElmtSoln.m_QpU.clear();
    m_LocalElmtSoln.m_QpUold.clear();
    m_LocalElmtSoln.m_QpUolder.clear();
    m_LocalElmtSoln.m_QpV.clear();
    m_LocalElmtSoln.m_QpA.clear();

    m_LocalElmtSoln.m_QpGradU.clear();
    m_LocalElmtSoln.m_QpGradUold.clear();
    m_LocalElmtSoln.m_QpGradUolder.clear();

    m_LocalElmtSoln.m_QpGradV.clear();

    m_LocalElmtSoln.m_Qpgradu.clear();
    m_LocalElmtSoln.m_Qpgradv.clear();

    m_LocalElmtInfo.m_QpCoords0=0.0;
    m_LocalElmtInfo.m_QpCoords =0.0;

    m_Nodes.clear();
    m_Nodes0.clear();
    
}
//******************************************************************
void ProjectionSystem::addScalarMateName2List(const string &matename){
    if(m_Data.m_scalarmate_namelist.size()<1){
        m_Data.m_scalarmate_namelist.push_back(matename);
        m_Data.m_scalarmate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_scalarmate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_scalarmate_namelist.push_back(matename);
            m_Data.m_scalarmate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate scalar material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addVectorMateName2List(const string &matename){
    if(m_Data.m_vectormate_namelist.size()<1){
        m_Data.m_vectormate_namelist.push_back(matename);
        m_Data.m_vectormate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_vectormate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_vectormate_namelist.push_back(matename);
            m_Data.m_vectormate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate vector material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addRank2MateName2List(const string &matename){
    if(m_Data.m_rank2mate_namelist.size()<1){
        m_Data.m_rank2mate_namelist.push_back(matename);
        m_Data.m_rank2mate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_rank2mate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_rank2mate_namelist.push_back(matename);
            m_Data.m_rank2mate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate rank-2 tensor material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//*******************************************
void ProjectionSystem::addRank4MateName2List(const string &matename){
    if(m_Data.m_rank4mate_namelist.size()<1){
        m_Data.m_rank4mate_namelist.push_back(matename);
        m_Data.m_rank4mate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_rank4mate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_rank4mate_namelist.push_back(matename);
            m_Data.m_rank4mate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate rank-4 tensor material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//**********************************************************
void ProjectionSystem::printProjectionInfo()const{
    string str;
    char buff[69];
    MessagePrinter::printNormalTxt("Projection system information summary");

    snprintf(buff,69,"scalar mate=%2d, vector mate=%2d, rank-2 mate=%2d, rank-4 mate=%2d",
                     getScalarMaterialNum(),getVectorMaterialNum(),
                     getRank2MaterialNum(),getRank4MaterialNum());
    str=buff;
    MessagePrinter::printNormalTxt(str);

    str.clear();
    if(getScalarMaterialNum()<1){
        MessagePrinter::printNormalTxt("  scalar material = (empty, no projection)");
    }
    else{
        for(const auto &it:m_Data.m_scalarmate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  scalar material ="+str);
    }

    if(getVectorMaterialNum()<1){
        MessagePrinter::printNormalTxt("  vector material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_vectormate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  vector material ="+str);
    }

    if(getRank2MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-2 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_rank2mate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-2 tensor material ="+str);
    }

    if(getRank4MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-4 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_rank4mate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-4 tensor material ="+str);
    }
    
    MessagePrinter::printStars();
    
}