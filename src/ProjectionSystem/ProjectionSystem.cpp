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

#include "ProjectionSystem/ProjectionSystem.h"

ProjectionSystem::ProjectionSystem(){
    m_isallocated=false;/**< boolean flag for memory allocation status */
    m_isprojection=false;
    m_nodesnum=0;/** total nodes of bulk mesh */

    m_data.m_scalarmate_num=0;/**< number of scalar material to be projected */
    m_data.m_vectormate_num=0;/**< number of vector material to be projected */
    m_data.m_rank2mate_num=0;/**< number of rank2tensor material to be projected */
    m_data.m_rank4mate_num=0;/**< number of rank4tensor material to be projected */

    m_data.m_scalarmate_namelist.clear();/**< vector for the name of scalar materials to be projected */
    m_data.m_vectormate_namelist.clear();/**< vector for the name of vector materials to be projected */
    m_data.m_rank2mate_namelist.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_data.m_rank4mate_namelist.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_proj_type=ProjectionType::DEFAULT;/**< the type of projection method */
    m_data.m_proj_scalarmate_vec.releaseMemory();/**< the vector stores the projected scalar materials */
    m_data.m_proj_vectormate_vec.releaseMemory();/**< the vector stores the projected vector materials */
    m_data.m_proj_rank2mate_vec.releaseMemory();/**< the vector stores the projected rank-2 materials */
    m_data.m_proj_rank4mate_vec.releaseMemory();/**< the vector stores the projected rank-4 materials */
}
//******************************************************************
//*** init
//******************************************************************
void ProjectionSystem::init(const Mesh &t_mesh,const DofHandler &t_dofhandler){
    m_nodesnum=t_mesh.getBulkMeshNodesNum();
    if(getScalarMaterialNum()){
        m_data.m_proj_scalarmate_vec.resize(m_nodesnum*(1+1*m_data.m_scalarmate_num),0.0);
    }
    if(getVectorMaterialNum()){
        m_data.m_proj_vectormate_vec.resize(m_nodesnum*(1+3*m_data.m_vectormate_num),0.0);
    }
    if(getRank2MaterialNum()){
        m_data.m_proj_rank2mate_vec.resize(m_nodesnum*(1+9*m_data.m_rank2mate_num),0.0);
    }
    if(getRank4MaterialNum()){
        m_data.m_proj_rank4mate_vec.resize(m_nodesnum*(1+36*m_data.m_rank4mate_num),0.0);
    }

    int m_max_nodal_dofs=t_dofhandler.getMaxDofsPerNode();
    m_bulkelmt_nodesnum=t_mesh.getBulkMeshNodesNumPerBulkElmt();

    m_elmtconn.resize(m_bulkelmt_nodesnum,0);
    m_subelmtdofsid.resize(m_max_nodal_dofs+1,0);

    m_local_elmtsoln.m_gpU.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpUold.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpUolder.resize(m_max_nodal_dofs+1,0.0);

    m_local_elmtsoln.m_gpV.resize(m_max_nodal_dofs+1,0.0);
    m_local_elmtsoln.m_gpA.resize(m_max_nodal_dofs+1,0.0);

    m_local_elmtsoln.m_gpGradU.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpGradUold.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpGradUolder.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtsoln.m_gpGradV.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtsoln.m_gpgradu.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_local_elmtsoln.m_gpgradv.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_local_elmtinfo.m_gpCoords=0.0;
    m_local_elmtinfo.m_gpCoords0=0.0;

    m_nodes.resize(m_bulkelmt_nodesnum);
    m_nodes0.resize(m_bulkelmt_nodesnum);
}

void ProjectionSystem::makeGhostCopyOfProjectionData(){
    m_data.m_proj_scalarmate_vec.makeGhostCopy();
    m_data.m_proj_vectormate_vec.makeGhostCopy();
    m_data.m_proj_rank2mate_vec.makeGhostCopy();
    m_data.m_proj_rank4mate_vec.makeGhostCopy();
}
void ProjectionSystem::destroyGhostCopyOfProjectionData(){
    m_data.m_proj_scalarmate_vec.destroyGhostCopy();
    m_data.m_proj_vectormate_vec.destroyGhostCopy();
    m_data.m_proj_rank2mate_vec.destroyGhostCopy();
    m_data.m_proj_rank4mate_vec.destroyGhostCopy();
}
//******************************************************************
//*** release memory
//******************************************************************
void ProjectionSystem::releaseMemory(){
    m_isallocated=false;/**< boolean flag for memory allocation status */
    m_isprojection=false;

    m_nodesnum=0;/** total nodes of bulk mesh */

    m_data.m_scalarmate_num=0;/**< number of scalar material to be projected */
    m_data.m_vectormate_num=0;/**< number of vector material to be projected */
    m_data.m_rank2mate_num=0;/**< number of rank2tensor material to be projected */
    m_data.m_rank4mate_num=0;/**< number of rank4tensor material to be projected */

    m_data.m_scalarmate_namelist.clear();/**< vector for the name of scalar materials to be projected */
    m_data.m_vectormate_namelist.clear();/**< vector for the name of vector materials to be projected */
    m_data.m_rank2mate_namelist.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_data.m_rank4mate_namelist.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_proj_type=ProjectionType::DEFAULT;/**< the type of projection method */
    m_data.m_scalarmate_num=0;
    m_data.m_vectormate_num=0;
    m_data.m_rank2mate_num=0;
    m_data.m_rank4mate_num=0;
    m_data.m_proj_scalarmate_vec.releaseMemory();/**< the vector stores the projected scalar materials */
    m_data.m_proj_vectormate_vec.releaseMemory();/**< the vector stores the projected vector materials */
    m_data.m_proj_rank2mate_vec.releaseMemory();/**< the vector stores the projected rank-2 materials */
    m_data.m_proj_rank4mate_vec.releaseMemory();/**< the vector stores the projected rank-4 materials */

    m_local_elmtsoln.m_gpU.clear();
    m_local_elmtsoln.m_gpUold.clear();
    m_local_elmtsoln.m_gpUolder.clear();
    m_local_elmtsoln.m_gpV.clear();
    m_local_elmtsoln.m_gpA.clear();

    m_local_elmtsoln.m_gpGradU.clear();
    m_local_elmtsoln.m_gpGradUold.clear();
    m_local_elmtsoln.m_gpGradUolder.clear();

    m_local_elmtsoln.m_gpGradV.clear();

    m_local_elmtsoln.m_gpgradu.clear();
    m_local_elmtsoln.m_gpgradv.clear();

    m_local_elmtinfo.m_gpCoords0=0.0;
    m_local_elmtinfo.m_gpCoords =0.0;

    m_nodes.clear();
    m_nodes0.clear();
    
}
//******************************************************************
void ProjectionSystem::addScalarMateName2List(const string &matename){
    if(m_data.m_scalarmate_namelist.size()<1){
        m_data.m_scalarmate_namelist.push_back(matename);
        m_data.m_scalarmate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_data.m_scalarmate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_data.m_scalarmate_namelist.push_back(matename);
            m_data.m_scalarmate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate scalar material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addVectorMateName2List(const string &matename){
    if(m_data.m_vectormate_namelist.size()<1){
        m_data.m_vectormate_namelist.push_back(matename);
        m_data.m_vectormate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_data.m_vectormate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_data.m_vectormate_namelist.push_back(matename);
            m_data.m_vectormate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate vector material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addRank2MateName2List(const string &matename){
    if(m_data.m_rank2mate_namelist.size()<1){
        m_data.m_rank2mate_namelist.push_back(matename);
        m_data.m_rank2mate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_data.m_rank2mate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_data.m_rank2mate_namelist.push_back(matename);
            m_data.m_rank2mate_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate rank-2 tensor material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//*******************************************
void ProjectionSystem::addRank4MateName2List(const string &matename){
    if(m_data.m_rank4mate_namelist.size()<1){
        m_data.m_rank4mate_namelist.push_back(matename);
        m_data.m_rank4mate_num=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_data.m_rank4mate_namelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_data.m_rank4mate_namelist.push_back(matename);
            m_data.m_rank4mate_num+=1;
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
        for(const auto &it:m_data.m_scalarmate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  scalar material ="+str);
    }

    if(getVectorMaterialNum()<1){
        MessagePrinter::printNormalTxt("  vector material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_data.m_vectormate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  vector material ="+str);
    }

    if(getRank2MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-2 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_data.m_rank2mate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-2 tensor material ="+str);
    }

    if(getRank4MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-4 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_data.m_rank4mate_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-4 tensor material ="+str);
    }
    
    MessagePrinter::printStars();
    
}