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
//+++ Date   : 2020.11.30
//+++ Purpose: Loop over all the bulk elements, calculate the system
//+++          residual, and K matrix
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"
#include "NonlinearSolver/SNESSolver.h"

void BulkFESystem::formBulkFE(const FECalcType &t_CalcType,
                              const double &T,
                              const double &Dt,
                              const double (&Ctan)[3],
                              FECell &t_FECell,
                              const DofHandler &t_DofHandler,
                              FE &t_FE,
                              ElmtSystem &t_ElmtSystem,
                              MateSystem &t_MateSystem,
                              SolutionSystem &t_SolnSystem,
                              SparseMatrix &AMATRIX,
                              Vector &RHS){
    if(T||Dt||Ctan[0]||t_FECell.getFECellBulkElmtsNum()||t_DofHandler.getBulkElmtsNum()||t_FE.getMaxDim()||
       t_ElmtSystem.getBulkElmtBlocksNum()){}
    if(t_CalcType==FECalcType::COMPUTERESIDUAL){
        RHS.setToZero();
    }
    else if(t_CalcType==FECalcType::COMPUTEJACOBIAN){
        AMATRIX.setToZero();
    }
    else if (t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        RHS.setToZero();
        AMATRIX.setToZero();
    }

    // for the current and the previous steps' solution array
    t_SolnSystem.m_Utemp.makeGhostCopy();// we always use u_temp as u-current in FormBulkFE !!!
    t_SolnSystem.m_Uold.makeGhostCopy();
    t_SolnSystem.m_Uolder.makeGhostCopy();

    // for the velocity and acceleration
    t_SolnSystem.m_V.makeGhostCopy();
    t_SolnSystem.m_A.makeGhostCopy();


    double xi,eta,zeta,w,J,JxW;

    int SubElmtBlockID;
    int GlobalDofID,GlobalNodeID;
    int GlobalI,GlobalJ;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    vector<SingleMeshCell> MyLocalCellVec;
    MyLocalCellVec=t_FECell.getLocalBulkFECellVecCopy();
    m_LocalElmtInfo.m_ElmtsNum=t_FECell.getLocalFECellBulkElmtsNum();
    for(int e=1;e<=m_LocalElmtInfo.m_ElmtsNum;e++){
        m_LocalElmtInfo.m_Dim=MyLocalCellVec[e-1].Dim;
        m_LocalElmtInfo.m_DofsNum=static_cast<int>(MyLocalCellVec[e-1].ElmtDofIDs.size());
        m_LocalElmtInfo.m_Dt=Dt;
        m_LocalElmtInfo.m_T=T;
        m_LocalElmtInfo.m_ElmtID=e;
        m_LocalElmtInfo.m_NodesNum=MyLocalCellVec[e-1].NodesNumPerElmt;
        m_BulkElmtNodesNum=MyLocalCellVec[e-1].NodesNumPerElmt;

        m_Nodes=MyLocalCellVec[e-1].ElmtNodeCoords;
        m_ElmtConn=MyLocalCellVec[e-1].ElmtConn;
        t_DofHandler.getIthLocalBulkElmtDofIDs(e,m_ElmtDofIDs);
        for(int i=0;i<t_DofHandler.getMaxDofsPerElmt();i++){
            m_ElmtUolder[i]=t_SolnSystem.m_Uolder.getIthValueFromGhost(m_ElmtDofIDs[i]);
            m_ElmtUold[i]  =t_SolnSystem.m_Uold.getIthValueFromGhost(m_ElmtDofIDs[i]);
            m_ElmtU[i]     =t_SolnSystem.m_Utemp.getIthValueFromGhost(m_ElmtDofIDs[i]);
            
            m_ElmtV[i]=t_SolnSystem.m_V.getIthValueFromGhost(m_ElmtDofIDs[i]);
            m_ElmtA[i]=t_SolnSystem.m_A.getIthValueFromGhost(m_ElmtDofIDs[i]);
        }

        /**
         * Now, we execute the qpoint integration
         */
        m_LocalElmtInfo.m_QpointsNum=t_FE.m_BulkQpoints.getQPointsNum();
        MyLocalCellVec[e-1].Volume=0.0;
        for(int qp=1;qp<=m_LocalElmtInfo.m_QpointsNum;qp++){
            m_LocalElmtInfo.m_QpointID=qp;
            w=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,0);
            xi=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,1);
            if (m_LocalElmtInfo.m_Dim==1){
                eta=0.0;
                zeta=0.0;
            }
            else if(m_LocalElmtInfo.m_Dim==2) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=0.0;
            }
            else if(m_LocalElmtInfo.m_Dim==3) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,3);
            }

            t_FE.m_BulkShp.calc(xi,eta,zeta,m_Nodes,true);
            J=t_FE.m_BulkShp.getJacDet();
            JxW=J*w;
            MyLocalCellVec[e-1].Volume+=1.0*JxW;/**< for the volume of the current fe cell */

            m_LocalElmtInfo.m_QpCoords=0.0;
            for (int i=1;i<=m_LocalElmtInfo.m_NodesNum;i++) {
                m_LocalElmtInfo.m_QpCoords(1)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes(i,1);
                m_LocalElmtInfo.m_QpCoords(2)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes(i,2);
                m_LocalElmtInfo.m_QpCoords(3)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes(i,3);
            }

            //*********************************************************************
            //*** loop over all the sub element/modulus of current bulk element
            //*********************************************************************
            if (t_CalcType==FECalcType::COMPUTERESIDUAL) {
                m_LocalR.setToZero();
            }
            else if (t_CalcType==FECalcType::COMPUTEJACOBIAN) {
                m_LocalK.setToZero();
            }
            else if (t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                m_LocalR.setToZero();
                m_LocalK.setToZero();
            }

            for (int SubElmt=1;SubElmt<=t_ElmtSystem.getLocalIthBulkElmtSubElmtsNum(e);SubElmt++) {
                if (t_CalcType==FECalcType::COMPUTERESIDUAL) {
                    m_SubR.setToZero();
                }
                else if (t_CalcType==FECalcType::COMPUTEJACOBIAN) {
                    m_SubK.setToZero();
                }
                else if (t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                    m_SubR.setToZero();
                    m_SubK.setToZero();
                }

                SubElmtBlockID=t_ElmtSystem.getLocalIthBulkElmtJthSubElmtID(e,SubElmt);
                m_LocalElmtInfo.m_DofsNum=static_cast<int>(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_DofIDs.size());
                m_SubElmtDofs=m_LocalElmtInfo.m_DofsNum;

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
                        m_LocalElmtSoln.m_QpU[i]     +=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Utemp.getIthValueFromGhost(GlobalDofID);

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
                        m_LocalElmtSoln.m_QpGradU[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Utemp.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Utemp.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Utemp.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpGradV[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                    }
                }// end-of-sub-element-dofs-loop

                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                if (t_CalcType==FECalcType::INITMATERIAL) {
                    t_MateSystem.initBulkMateLibs(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_MateType,
                                                  t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_JsonParams,
                                                  m_LocalElmtInfo,
                                                  m_LocalElmtSoln);


                    t_SolnSystem.m_QpointsScalarMaterials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getScalarMaterialsCopy();
                    t_SolnSystem.m_QpointsVectorMaterials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getVectorMaterialsCopy();
                    t_SolnSystem.m_QpointsRank2Materials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank2MaterialsCopy();
                    t_SolnSystem.m_QpointsRank4Materials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank4MaterialsCopy();

                    t_SolnSystem.m_QpointsScalarMaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getScalarMaterialsCopy();
                    t_SolnSystem.m_QpointsVectorMaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getVectorMaterialsCopy();
                    t_SolnSystem.m_QpointsRank2MaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank2MaterialsCopy();
                    t_SolnSystem.m_QpointsRank4MaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank4MaterialsCopy();

                }
                else {
                    t_MateSystem.m_MaterialContainerOld.getScalarMaterialsRef()=t_SolnSystem.m_QpointsScalarMaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1];
                    t_MateSystem.m_MaterialContainerOld.getVectorMaterialsRef()=t_SolnSystem.m_QpointsVectorMaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1];
                    t_MateSystem.m_MaterialContainerOld.getRank2MaterialsRef()=t_SolnSystem.m_QpointsRank2MaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1];
                    t_MateSystem.m_MaterialContainerOld.getRank4MaterialsRef()=t_SolnSystem.m_QpointsRank4MaterialsOld_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1];

                    t_MateSystem.runBulkMateLibs(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_MateType,
                                                 t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_JsonParams,
                                                 m_LocalElmtInfo,
                                                 m_LocalElmtSoln);

                    t_SolnSystem.m_QpointsScalarMaterials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getScalarMaterialsCopy();
                    t_SolnSystem.m_QpointsVectorMaterials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getVectorMaterialsCopy();
                    t_SolnSystem.m_QpointsRank2Materials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank2MaterialsCopy();
                    t_SolnSystem.m_QpointsRank4Materials_Local[(e-1)*m_LocalElmtInfo.m_QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank4MaterialsCopy();
                }

                if (t_CalcType==FECalcType::COMPUTERESIDUAL||
                    t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN){
                    for (int i=1;i<=m_BulkElmtNodesNum;i++) {
                        m_LocalShp.m_Test=t_FE.m_BulkShp.shape_value(i);
                        m_LocalShp.m_GradTest=t_FE.m_BulkShp.shape_grad(i);
                        GlobalI=MyLocalCellVec[e-1].ElmtConn[i-1];

                        t_ElmtSystem.runBulkElmtLibs(t_CalcType,
                                                     Ctan,
                                                     SubElmtBlockID,
                                                     t_MateSystem.m_MaterialContainerOld,
                                                     t_MateSystem.m_MaterialContainer,
                                                     m_LocalElmtInfo,
                                                     m_LocalElmtSoln,
                                                     m_LocalShp,
                                                     m_SubK,
                                                     m_SubR);
                        assembleLocalResidual2GlobalR(m_SubElmtDofs,
                                                      m_SubElmtDofIDs,
                                                      GlobalI,
                                                      t_DofHandler,
                                                      JxW,m_SubR,
                                                      RHS);
                        if (t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                            for (int j=1;j<=m_BulkElmtNodesNum;j++) {
                                m_LocalShp.m_Trial=t_FE.m_BulkShp.shape_value(j);
                                m_LocalShp.m_GradTrial=t_FE.m_BulkShp.shape_grad(j);
                                GlobalJ=MyLocalCellVec[e-1].ElmtConn[j-1];
                                t_ElmtSystem.runBulkElmtLibs(t_CalcType,
                                                             Ctan,
                                                             SubElmtBlockID,
                                                             t_MateSystem.m_MaterialContainerOld,
                                                             t_MateSystem.m_MaterialContainer,
                                                             m_LocalElmtInfo,
                                                             m_LocalElmtSoln,
                                                             m_LocalShp,
                                                             m_SubK,
                                                             m_SubR);
                                assembleLocalJacobian2GlobalK(m_SubElmtDofs,m_SubElmtDofIDs,GlobalI,GlobalJ,JxW,t_DofHandler,m_SubK,AMATRIX);
                            } // end-of-J-loop
                        }
                    } // end-of-I-loop
                }// end-of-residual-calculation
                else if (t_CalcType==FECalcType::COMPUTEJACOBIAN) {
                    for (int i=1;i<=m_BulkElmtNodesNum;i++) {
                        m_LocalShp.m_Test=t_FE.m_BulkShp.shape_value(i);
                        m_LocalShp.m_GradTest=t_FE.m_BulkShp.shape_grad(i);
                        GlobalI=MyLocalCellVec[e-1].ElmtConn[i-1];
                        for (int j=1;j<=m_BulkElmtNodesNum;j++) {
                            m_LocalShp.m_Trial=t_FE.m_BulkShp.shape_value(j);
                            m_LocalShp.m_GradTrial=t_FE.m_BulkShp.shape_grad(j);
                            GlobalJ=MyLocalCellVec[e-1].ElmtConn[j-1];
                            t_ElmtSystem.runBulkElmtLibs(t_CalcType,
                                                         Ctan,
                                                         SubElmtBlockID,
                                                         t_MateSystem.m_MaterialContainerOld,
                                                         t_MateSystem.m_MaterialContainer,
                                                         m_LocalElmtInfo,
                                                         m_LocalElmtSoln,
                                                         m_LocalShp,
                                                         m_SubK,
                                                         m_SubR);
                            assembleLocalJacobian2GlobalK(m_SubElmtDofs,m_SubElmtDofIDs,GlobalI,GlobalJ,JxW,t_DofHandler,m_SubK,AMATRIX);
                        } // end-of-J-loop
                    } // end-of-I-loop
                }// end-of-jacobian-calculation
            } // end-of-sub-element-loop
        }// end-of-qpoints-loop
    }// end-of-local-element-loop


    // finish the final assemble
    if(t_CalcType==FECalcType::COMPUTERESIDUAL) RHS.assemble();
    if(t_CalcType==FECalcType::COMPUTEJACOBIAN) AMATRIX.assemble();
    if (t_CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        RHS.assemble();
        AMATRIX.assemble();
    }


    // release the ghost copy
    t_SolnSystem.m_Utemp.destroyGhostCopy();
    t_SolnSystem.m_Uold.destroyGhostCopy();
    t_SolnSystem.m_Uolder.destroyGhostCopy();
    t_SolnSystem.m_V.destroyGhostCopy();
    t_SolnSystem.m_A.destroyGhostCopy();

}
