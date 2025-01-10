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
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply dirichlet boundary condition via the 
//+++          penalty method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyDirichletBC(const FECalcType &CalcType,
                                const double &BCValue,
                                const BCType &t_BCType,
                                const nlohmann::json &Params,
                                const vector<int> &DofIDs,
                                const vector<string> &BCNameList,
                                const FECell &t_FECell,
                                const DofHandler &t_DofHandler,
                                Vector &U,
                                Vector &Ucopy,
                                Vector &Uold,
                                Vector &Uolder,
                                Vector &V,
                                SparseMatrix &AMATRIX,
                                Vector &RHS){
    int i,j,k,iInd;
    vector<int> globaldofids;
    globaldofids.resize(DofIDs.size(),0);

    m_LocalElmtInfo.m_DofsNum=static_cast<int>(DofIDs.size());

    Ucopy.makeGhostCopy();
    Uold.makeGhostCopy();
    Uolder.makeGhostCopy();
    V.makeGhostCopy();

    for(const auto &name:BCNameList){
        for (const auto &cell:t_FECell.getLocalMeshCellVectorCopyViaPhyName(name)){
            m_LocalElmtInfo.m_Dim=cell.Dim;
            m_LocalElmtInfo.m_DofsNum=static_cast<int>(DofIDs.size());
            for(i=1;i<=cell.NodesNumPerElmt;i++){
                j=cell.ElmtConn[i-1];
                m_LocalElmtInfo.m_QpCoords0(1)=cell.ElmtNodeCoords0(i,1);
                m_LocalElmtInfo.m_QpCoords0(2)=cell.ElmtNodeCoords0(i,2);
                m_LocalElmtInfo.m_QpCoords0(3)=cell.ElmtNodeCoords0(i,3);
                for(k=1;k<=m_LocalElmtInfo.m_DofsNum;k++){
                    iInd=t_DofHandler.getIthNodeJthDofID(j,DofIDs[k-1]);
                    globaldofids[k-1]=iInd;

                    m_LocalElmtSoln.m_QpU[k]=Ucopy.getIthValueFromGhost(iInd);
                    m_LocalElmtSoln.m_QpUold[k]=Uold.getIthValueFromGhost(iInd);
                    m_LocalElmtSoln.m_QpUolder[k]=Uolder.getIthValueFromGhost(iInd);
                    m_LocalElmtSoln.m_QpV[k]=V.getIthValueFromGhost(iInd);
                }

                switch (t_BCType)
                {
                case BCType::DIRICHLETBC:
                    DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::ROTATEDDIRICHLETBC:
                    RotatedDirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER1DIRICHLETBC:
                    User1DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER2DIRICHLETBC:
                    User2DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER3DIRICHLETBC:
                    User3DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER4DIRICHLETBC:
                    User4DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER5DIRICHLETBC:
                    User5DirichletBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::POISSON2DBENCHMARKBC:
                    Poisson2DBenchmarkBC::computeBCValue(CalcType,m_DirichletPenalty,BCValue,Params,
                                                m_LocalElmtInfo,m_LocalElmtSoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                default:
                    MessagePrinter::printErrorTxt("unsupported dirichlet type boundary condition in ApplyDirichletBC.cpp, plese check your input file or your code");
                    MessagePrinter::exitAsFem();
                    break;
                }
            }
        }

    }

    Ucopy.destroyGhostCopy();
    Uold.destroyGhostCopy();
    Uolder.destroyGhostCopy();
    V.destroyGhostCopy();

    // do the assemble
    U.assemble();
    if(CalcType==FECalcType::COMPUTERESIDUAL) RHS.assemble();
    if(CalcType==FECalcType::COMPUTEJACOBIAN) AMATRIX.assemble();


}