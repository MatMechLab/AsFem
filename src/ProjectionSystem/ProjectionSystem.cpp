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

    m_Data.m_ScalarProjMateNum=0;/**< number of scalar material to be projected */
    m_Data.m_VectorProjMateNum=0;/**< number of vector material to be projected */
    m_Data.m_Rank2ProjMateNum=0;/**< number of rank2tensor material to be projected */
    m_Data.m_Rank4ProjMateNum=0;/**< number of rank4tensor material to be projected */

    m_Data.m_ScalarProjMateNameList.clear();/**< vector for the name of scalar materials to be projected */
    m_Data.m_VectorProjMateNamelist.clear();/**< vector for the name of vector materials to be projected */
    m_Data.m_Rank2ProjMateNameList.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_Data.m_Rank4ProjMateNameList.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_ProjType=ProjectionType::DEFAULT;/**< the type of projection method */
    for (auto &it:m_Data.m_ScalarProjMateVecList) it.releaseMemory();
    for (auto &it:m_Data.m_VectorProjMateVecList) it.releaseMemory();
    for (auto &it:m_Data.m_Rank2ProjMateVecList) it.releaseMemory();
    for (auto &it:m_Data.m_Rank4ProjMateVecList) it.releaseMemory();
}
//******************************************************************
//*** init
//******************************************************************
void ProjectionSystem::init(const FECell &t_fecell,const DofHandler &t_dofhandler){
    if (m_ProjType==ProjectionType::DEFAULT||
        m_ProjType==ProjectionType::LEASTSQUARE) {
        LeastSquareProjection::initMyProjection(t_fecell,t_dofhandler);
    }
    else if (m_ProjType==ProjectionType::FULLLEASTSQUARE) {
        FullLeastSquareProjection::initMyProjection(t_fecell,t_dofhandler);
    }

    m_NodesNum=t_fecell.getFECellNodesNum();
    if(getScalarMaterialNum()){
        m_Data.m_ScalarProjMateVecList.resize(m_Data.m_ScalarProjMateNum);
        for (auto &it:m_Data.m_ScalarProjMateVecList) it.resize(m_NodesNum*(1+1),0.0);
    }
    if(getVectorMaterialNum()){
        m_Data.m_VectorProjMateVecList.resize(m_Data.m_VectorProjMateNum,0.0);
        for (auto &it:m_Data.m_VectorProjMateVecList) it.resize(m_NodesNum*(1+3),0.0);
    }
    if(getRank2MaterialNum()){
        m_Data.m_Rank2ProjMateVecList.resize(m_Data.m_Rank2ProjMateNum);
        for (auto &it:m_Data.m_Rank2ProjMateVecList) it.resize(m_NodesNum*(1+9),0.0);
    }
    if(getRank4MaterialNum()){
        m_Data.m_Rank4ProjMateVecList.resize(m_Data.m_Rank4ProjMateNum);
        for (auto &it:m_Data.m_Rank4ProjMateVecList) it.resize(m_NodesNum*(1+36),0.0);
    }
}

void ProjectionSystem::makeGhostCopyOfProjectionData(){
    for (auto &it:m_Data.m_ScalarProjMateVecList) it.makeGhostCopy();
    for (auto &it:m_Data.m_VectorProjMateVecList) it.makeGhostCopy();
    for (auto &it:m_Data.m_Rank2ProjMateVecList) it.makeGhostCopy();
    for (auto &it:m_Data.m_Rank4ProjMateVecList) it.makeGhostCopy();
}
void ProjectionSystem::destroyGhostCopyOfProjectionData(){
    for (auto &it:m_Data.m_ScalarProjMateVecList) it.destroyGhostCopy();
    for (auto &it:m_Data.m_VectorProjMateVecList) it.destroyGhostCopy();
    for (auto &it:m_Data.m_Rank2ProjMateVecList) it.destroyGhostCopy();
    for (auto &it:m_Data.m_Rank4ProjMateVecList) it.destroyGhostCopy();
}
//******************************************************************
//*** release memory
//******************************************************************
void ProjectionSystem::releaseMemory(){
    m_IsAllocated=false;/**< boolean flag for memory allocation status */
    m_IsProjection=false;

    m_NodesNum=0;/** total nodes of bulk mesh */

    m_Data.m_ScalarProjMateNum=0;/**< number of scalar material to be projected */
    m_Data.m_VectorProjMateNum=0;/**< number of vector material to be projected */
    m_Data.m_Rank2ProjMateNum=0;/**< number of rank2tensor material to be projected */
    m_Data.m_Rank4ProjMateNum=0;/**< number of rank4tensor material to be projected */

    m_Data.m_ScalarProjMateNameList.clear();/**< vector for the name of scalar materials to be projected */
    m_Data.m_VectorProjMateNamelist.clear();/**< vector for the name of vector materials to be projected */
    m_Data.m_Rank2ProjMateNameList.clear();/**< vector for the name of rank-2 tensor materials to be projected */
    m_Data.m_Rank4ProjMateNameList.clear();/**< vector for the name of rank-4 tensor materials to be projected */

    m_ProjType=ProjectionType::DEFAULT;/**< the type of projection method */
    for (auto &it:m_Data.m_ScalarProjMateVecList) it.releaseMemory();
    for (auto &it:m_Data.m_VectorProjMateVecList) it.releaseMemory();
    for (auto &it:m_Data.m_Rank2ProjMateVecList) it.releaseMemory();
    for (auto it:m_Data.m_Rank4ProjMateVecList) it.releaseMemory();

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

    if (m_ProjType==ProjectionType::FULLLEASTSQUARE) {

    }
    
}
//******************************************************************
void ProjectionSystem::addScalarMateName2List(const string &matename){
    if(m_Data.m_ScalarProjMateNameList.size()<1){
        m_Data.m_ScalarProjMateNameList.push_back(matename);
        m_Data.m_ScalarProjMateNum=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_ScalarProjMateNameList){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_ScalarProjMateNameList.push_back(matename);
            m_Data.m_ScalarProjMateNum+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate scalar material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addVectorMateName2List(const string &matename){
    if(m_Data.m_VectorProjMateNamelist.size()<1){
        m_Data.m_VectorProjMateNamelist.push_back(matename);
        m_Data.m_VectorProjMateNum=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_VectorProjMateNamelist){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_VectorProjMateNamelist.push_back(matename);
            m_Data.m_VectorProjMateNum+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate vector material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//************************************
void ProjectionSystem::addRank2MateName2List(const string &matename){
    if(m_Data.m_Rank2ProjMateNameList.size()<1){
        m_Data.m_Rank2ProjMateNameList.push_back(matename);
        m_Data.m_Rank2ProjMateNum=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_Rank2ProjMateNameList){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_Rank2ProjMateNameList.push_back(matename);
            m_Data.m_Rank2ProjMateNum+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate rank-2 tensor material name("+matename+"), can\'t add it to the projection list");
            MessagePrinter::exitAsFem();
        }
    }
}
//*******************************************
void ProjectionSystem::addRank4MateName2List(const string &matename){
    if(m_Data.m_Rank4ProjMateNameList.size()<1){
        m_Data.m_Rank4ProjMateNameList.push_back(matename);
        m_Data.m_Rank4ProjMateNum=1;
    }
    else{
        bool NoInList=true;
        for(const auto &it:m_Data.m_Rank4ProjMateNameList){
            if(matename==it){
                NoInList=false;break;
            }
        }
        if(NoInList){
            m_Data.m_Rank4ProjMateNameList.push_back(matename);
            m_Data.m_Rank4ProjMateNum+=1;
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

    if (m_ProjType==ProjectionType::DEFAULT||
        m_ProjType==ProjectionType::LEASTSQUARE) {
        MessagePrinter::printNormalTxt("  projection method= simplified least square (default)");
        }
    else if (m_ProjType==ProjectionType::FULLLEASTSQUARE) {
        MessagePrinter::printNormalTxt("  projection method= full least square");
    }

    snprintf(buff,69,"  scalar mate=%2d, vector mate=%2d, rank-2 mate=%2d, rank-4 mate=%2d",
                     getScalarMaterialNum(),getVectorMaterialNum(),
                     getRank2MaterialNum(),getRank4MaterialNum());
    str=buff;
    MessagePrinter::printNormalTxt(str);

    str.clear();
    if(getScalarMaterialNum()<1){
        MessagePrinter::printNormalTxt("  scalar material = (empty, no projection)");
    }
    else{
        for(const auto &it:m_Data.m_ScalarProjMateNameList) str+=it+" ";
        MessagePrinter::printNormalTxt("  scalar material ="+str);
    }

    if(getVectorMaterialNum()<1){
        MessagePrinter::printNormalTxt("  vector material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_VectorProjMateNamelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  vector material ="+str);
    }

    if(getRank2MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-2 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_Rank2ProjMateNameList) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-2 tensor material ="+str);
    }

    if(getRank4MaterialNum()<1){
        MessagePrinter::printNormalTxt("  rank-4 tensor material = (empty, no projection)");
    }
    else{
        str.clear();
        for(const auto &it:m_Data.m_Rank4ProjMateNameList) str+=it+" ";
        MessagePrinter::printNormalTxt("  rank-4 tensor material ="+str);
    }
    
    MessagePrinter::printStars();
    
}