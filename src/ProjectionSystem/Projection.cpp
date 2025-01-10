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
//+++ Date   : 2022.08.22
//+++ Purpose: execute the projection from integration points to
//+++          nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/ProjectionSystem.h"

void ProjectionSystem::executeProjection(const FECell &t_FECell,
                                         const DofHandler &t_DofHandler,
                                         const ElmtSystem &t_ElmtSystem,
                                         MateSystem &t_MateSystem,
                                         FE &t_FE,
                                         SolutionSystem &t_SolnSystem,
                                         const FEControlInfo &t_FECtrlInfo){
    // if no projection is required, then return back
    if(getScalarMaterialNum()<1 &&
       getVectorMaterialNum()<1 &&
       getRank2MaterialNum()<1 &&
       getRank4MaterialNum()<1){
        return;
    }
    
    // before we do the projection, set all the vector to zeros
    m_Data.m_proj_scalarmate_vec.setToZero();
    m_Data.m_proj_vectormate_vec.setToZero();
    m_Data.m_proj_rank2mate_vec.setToZero();
    m_Data.m_proj_rank4mate_vec.setToZero();

    // for the current and the previous steps' solution array
    t_SolnSystem.m_Ucurrent.makeGhostCopy();//
    t_SolnSystem.m_Uold.makeGhostCopy();
    t_SolnSystem.m_Uolder.makeGhostCopy();

    // for the velocity and acceleration
    t_SolnSystem.m_V.makeGhostCopy();
    t_SolnSystem.m_A.makeGhostCopy();

    int QpointsNum;
    double xi,eta,zeta,w,J,JxW;

    int SubElmtBlockID;
    int GlobalDofID,GlobalNodeID;

    m_LocalElmtInfo.m_Dim=t_FECell.getFECellMaxDim();
    m_LocalElmtInfo.m_NodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    vector<SingleMeshCell> MyLocalCellVec;

    MyLocalCellVec=t_FECell.getLocalBulkFECellVecCopy();
    m_BulkElmtNodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    QpointsNum=t_FE.m_BulkQpoints.getQPointsNum();
    m_LocalElmtInfo.m_QpointsNum=QpointsNum;

    for (int e=1;e<=t_FECell.getLocalFECellBulkElmtsNum();e++) {
        m_LocalElmtInfo.m_Dim=MyLocalCellVec[e-1].Dim;
        m_LocalElmtInfo.m_NodesNum=MyLocalCellVec[e-1].NodesNumPerElmt;
        m_LocalElmtInfo.m_Dt=t_FECtrlInfo.Dt;
        m_LocalElmtInfo.m_T=t_FECtrlInfo.T;

        m_Nodes=MyLocalCellVec[e-1].ElmtNodeCoords;
        m_Nodes0=MyLocalCellVec[e-1].ElmtNodeCoords0;
        m_ElmtConn=MyLocalCellVec[e-1].ElmtConn;
        for (int qp=1;qp<=QpointsNum;qp++) {
            w=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,0);
            xi=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,1);
            if (m_LocalElmtInfo.m_Dim==1) {
                eta=0.0;zeta=0.0;
            }
            else if (m_LocalElmtInfo.m_Dim==2) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=0.0;
            }
            else if (m_LocalElmtInfo.m_Dim==3) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,3);
            }
            t_FE.m_BulkShp.calc(xi,eta,zeta,m_Nodes0,true);
            J=t_FE.m_BulkShp.getJacDet();
            JxW=J*w;

            m_LocalElmtInfo.m_QpointID=qp;
            m_LocalElmtInfo.m_QpCoords0=0.0;
            for (int i=1;i<=m_LocalElmtInfo.m_NodesNum;i++) {
                m_LocalElmtInfo.m_QpCoords0(1)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,1);
                m_LocalElmtInfo.m_QpCoords0(2)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,2);
                m_LocalElmtInfo.m_QpCoords0(3)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,3);
            }

            t_MateSystem.m_MaterialContainerOld.getScalarMaterialsRef()=t_SolnSystem.m_QpointsScalarMaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getVectorMaterialsRef()=t_SolnSystem.m_QpointsVectorMaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getRank2MaterialsRef()=t_SolnSystem.m_QpointsRank2MaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getRank4MaterialsRef()=t_SolnSystem.m_QpointsRank4MaterialsOld_Local[(e-1)*QpointsNum+qp-1];

            for (int SubElmt=1;SubElmt<=t_ElmtSystem.getLocalIthBulkElmtSubElmtsNum(e);SubElmt++) {
                SubElmtBlockID=t_ElmtSystem.getLocalIthBulkElmtJthSubElmtID(e,SubElmt);
                m_SubElmtDofs=static_cast<int>(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_DofIDs.size());
                m_LocalElmtInfo.m_DofsNum=m_SubElmtDofs;

                //*********************************************************************
                //*** calculate physical quantities on each integration point
                //*********************************************************************
                // here, the U/V/A vector's index should start from 1, to make it consistent with
                // 1st dof, 2nd dof, 3rd one , and so on. So, please let the index starts from 1 !!!
                for (int i=1;i<=m_LocalElmtInfo.m_DofsNum;i++){
                    m_SubElmtDofIDs[i-1]=t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_DofIDs[i-1];

                    m_LocalElmtSoln.m_QpU[i]=0.0;
                    m_LocalElmtSoln.m_QpUold[i]=0.0;
                    m_LocalElmtSoln.m_QpUolder[i]=0.0;

                    m_LocalElmtSoln.m_QpV[i]=0.0;
                    m_LocalElmtSoln.m_QpA[i]=0.0;

                    m_LocalElmtSoln.m_QpGradU[i]=0.0;
                    m_LocalElmtSoln.m_QpGradUold[i]=0.0;
                    m_LocalElmtSoln.m_QpGradUolder[i]=0.0;

                    m_LocalElmtSoln.m_QpGradV[i]=0.0;
                    for (int j=1;j<=m_BulkElmtNodesNum;j++) {
                        GlobalNodeID=MyLocalCellVec[e-1].ElmtConn[j-1];
                        GlobalDofID=t_DofHandler.getIthNodeJthDofID(GlobalNodeID,m_SubElmtDofIDs[i-1]);

                        m_LocalElmtSoln.m_QpUolder[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpUold[i]  +=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpU[i]     +=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpV[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpA[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_A.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpGradUolder[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUolder[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUolder[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        //
                        m_LocalElmtSoln.m_QpGradUold[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUold[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUold[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        //
                        m_LocalElmtSoln.m_QpGradU[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpGradV[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                    }
                }// end-of-sub-element-dofs-loop


                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                t_MateSystem.runBulkMateLibs(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_MateType,
                                             t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_JsonParams,
                                             m_LocalElmtInfo,
                                             m_LocalElmtSoln);

            }//end-of-sub-element-loop
            // true---> for the local projection
            runProjectionLibs(true,t_FECell,m_BulkElmtNodesNum,m_ElmtConn,JxW,t_FE.m_BulkShp,t_MateSystem.m_MaterialContainer,m_Data);
        }// end-of-qpoints-loop
    }// end-of-element-loop

    t_SolnSystem.m_Ucurrent.destroyGhostCopy();//
    t_SolnSystem.m_Uold.destroyGhostCopy();
    t_SolnSystem.m_Uolder.destroyGhostCopy();

    // for the velocity and acceleration
    t_SolnSystem.m_V.destroyGhostCopy();
    t_SolnSystem.m_A.destroyGhostCopy();

    //******************************************************
    //*** after the local projection(in the element loop)
    //*** one needs to execute the global projection behavior
    //******************************************************
    // false--> for global projection action
    runProjectionLibs(false,t_FECell,m_BulkElmtNodesNum,m_ElmtConn,JxW,t_FE.m_BulkShp,t_MateSystem.m_MaterialContainer,m_Data);

}