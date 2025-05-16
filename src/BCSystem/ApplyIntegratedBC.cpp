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
//+++ Purpose: apply the integrated boundary conditions with 'surface'
//+++          integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyIntegratedBC(const FECalcType &CalcType,
                                 const double &BCValue,
                                 const BCType &t_BCType,
                                 const double (&Ctan)[3],
                                 const nlohmann::json &Parameters,
                                 const vector<int> &DofIDs,
                                 const vector<string> &BCNameList,
                                 const FECell &t_FECell,const DofHandler &t_DofHandler,
                                 FE &t_FE,
                                 Vector &U,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS){
    // for other types of boundary conditions
    int nNodesPerBCElmt;
    int GlobalDofID,i,j,k,iInd,jInd,gpInd;
    double xi,eta,w,JxW,dist;

    if(t_FECell.getFECellBulkElmtsNum()){}

    m_LocalElmtInfo.m_DofsNum=static_cast<int>(DofIDs.size());

    U.makeGhostCopy();
    Uold.makeGhostCopy();
    Uolder.makeGhostCopy();
    V.makeGhostCopy();    

    for(const auto &name:BCNameList){
        for (const auto &cell:t_FECell.getLocalMeshCellVectorCopyViaPhyName(name)) {
            nNodesPerBCElmt=cell.NodesNumPerElmt;
            m_LocalElmtInfo.m_Dim=cell.Dim;
            if(m_LocalElmtInfo.m_Dim>2 || m_LocalElmtInfo.m_Dim<0){
                MessagePrinter::printErrorTxt("Invalid dim (="+to_string(m_LocalElmtInfo.m_Dim)+") for integrated boundary condition");
                MessagePrinter::exitAsFem();
            }
            if(m_LocalElmtInfo.m_Dim==0){
                 // for 'point' case
                for(i=1;i<=nNodesPerBCElmt;i++){
                    iInd=cell.ElmtConn[i-1];
                    m_LocalElmtInfo.m_QpCoords0(1)=cell.ElmtNodeCoords0(i,1);
                    m_LocalElmtInfo.m_QpCoords0(2)=cell.ElmtNodeCoords0(i,2);
                    m_LocalElmtInfo.m_QpCoords0(3)=cell.ElmtNodeCoords0(i,3);
                    for(k=0;k<m_LocalElmtInfo.m_DofsNum;k++){
                        GlobalDofID=t_DofHandler.getIthNodeJthDofID(iInd,k+1);
                        // get the local solutions
                        m_LocalElmtSoln.m_QpU[k+1]=U.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpUold[k+1]=Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpUolder[k+1]=Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpV[k+1]=V.getIthValueFromGhost(GlobalDofID);
                    }
                }

                JxW=1.0;

                // do the residual and jacobian calculation
                if(CalcType==FECalcType::COMPUTERESIDUAL||
                   CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN){
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_LocalShp.m_Test=1.0;
                        m_LocalShp.m_GradTest=0.0;
                        m_LocalShp.m_Trial=0.0;
                        m_LocalShp.m_GradTrial=0.0;
                        iInd=cell.ElmtConn[i-1];
                        runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                  m_Normal,
                                  m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                  m_LocalK,m_LocalR);
                        // for assemble
                        assembleLocalResidual2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,JxW,t_DofHandler,m_LocalR,RHS);
                        if (CalcType == FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                            for(j=1;j<=nNodesPerBCElmt;j++){
                                m_LocalShp.m_Trial=0.0;
                                m_LocalShp.m_GradTrial=0.0;
                                jInd=cell.ElmtConn[j-1];
                                runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                          m_Normal,
                                          m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                          m_LocalK,m_LocalR);
                                // for jacobian assemble
                                assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                            }
                        }
                    }
                }
                else if(CalcType==FECalcType::COMPUTEJACOBIAN){
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_LocalShp.m_Test=1.0;
                        m_LocalShp.m_GradTest=0.0;
                        iInd=cell.ElmtConn[i-1];
                        for(j=1;j<=nNodesPerBCElmt;j++){
                            m_LocalShp.m_Trial=0.0;
                            m_LocalShp.m_GradTrial=0.0;
                            jInd=cell.ElmtConn[j-1];
                            runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                      m_Normal,
                                      m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                      m_LocalK,m_LocalR);
                            // for jacobian assemble
                            assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                        }
                    }
                }// end-of-K-and-R-calculation

            }// end-of-dim=0-case
            else if(m_LocalElmtInfo.m_Dim==1){
                // for 'line' element case
                // do the gauss point integration loop
                for(gpInd=1;gpInd<=t_FE.m_LineQpoints.getQPointsNum();gpInd++){
                    xi =t_FE.m_LineQpoints.getIthPointJthCoord(gpInd,1);
                    eta=0.0;
                    w  =t_FE.m_LineQpoints.getIthPointJthCoord(gpInd,0);
                    t_FE.m_LineShp.calc(xi,eta,0.0,cell.ElmtNodeCoords0,false);

                    // for the normal vector of current qpoint
                    m_XS.setToZeros();m_Normal=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_XS(1,1)+=t_FE.m_LineShp.shape_grad(i)(1)*cell.ElmtNodeCoords0(i,1);
                        m_XS(2,1)+=t_FE.m_LineShp.shape_grad(i)(1)*cell.ElmtNodeCoords0(i,2);
                    }
                    dist=sqrt(m_XS(1,1)*m_XS(1,1)+m_XS(2,1)*m_XS(2,1));
                    m_Normal(1)= m_XS(2,1)/dist;// dy/dxi
                    m_Normal(2)=-m_XS(1,1)/dist;// dx/dxi
                    m_Normal(3)= 0.0;

                    //********************************************************
                    //*** for the physical quantities on current qpoint
                    //********************************************************
                    t_FE.m_LineShp.calc(xi,eta,0.0,cell.ElmtNodeCoords0,true);
                    JxW=t_FE.m_LineShp.getJacDet()*w;

                    for(k=1;k<=m_LocalElmtInfo.m_DofsNum;k++){
                        m_LocalElmtSoln.m_QpU[k]=0.0;
                        m_LocalElmtSoln.m_QpUold[k]=0.0;
                        m_LocalElmtSoln.m_QpUolder[k]=0.0;
                        m_LocalElmtSoln.m_QpV[k]=0.0;

                        m_LocalElmtSoln.m_QpGradU[k]=0.0;
                        m_LocalElmtSoln.m_QpGradUold[k]=0.0;
                        m_LocalElmtSoln.m_QpGradUolder[k]=0.0;
                    }
                    m_LocalElmtInfo.m_QpCoords0=0.0;
                            
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        iInd=cell.ElmtConn[i-1];
                        for(k=1;k<=m_LocalElmtInfo.m_DofsNum;k++){
                            GlobalDofID=t_DofHandler.getIthNodeJthDofID(iInd,DofIDs[k-1]);

                            m_LocalElmtSoln.m_QpU[k]     +=t_FE.m_LineShp.shape_value(i)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpUold[k]  +=t_FE.m_LineShp.shape_value(i)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpUolder[k]+=t_FE.m_LineShp.shape_value(i)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpV[k]     +=t_FE.m_LineShp.shape_value(i)*V.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradU[k](1)+=t_FE.m_LineShp.shape_grad(i)(1)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradU[k](2)+=t_FE.m_LineShp.shape_grad(i)(2)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradU[k](3)+=t_FE.m_LineShp.shape_grad(i)(3)*U.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradUold[k](1)+=t_FE.m_LineShp.shape_grad(i)(1)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUold[k](2)+=t_FE.m_LineShp.shape_grad(i)(2)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUold[k](3)+=t_FE.m_LineShp.shape_grad(i)(3)*Uold.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradUolder[k](1)+=t_FE.m_LineShp.shape_grad(i)(1)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUolder[k](2)+=t_FE.m_LineShp.shape_grad(i)(2)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUolder[k](3)+=t_FE.m_LineShp.shape_grad(i)(3)*Uolder.getIthValueFromGhost(GlobalDofID);
                        }
                                
                        m_LocalElmtInfo.m_QpCoords0(1)+=t_FE.m_LineShp.shape_value(i)*cell.ElmtNodeCoords0(i,1);
                        m_LocalElmtInfo.m_QpCoords0(2)+=t_FE.m_LineShp.shape_value(i)*cell.ElmtNodeCoords0(i,2);
                        m_LocalElmtInfo.m_QpCoords0(3)+=t_FE.m_LineShp.shape_value(i)*cell.ElmtNodeCoords0(i,3);
                            
                    }// end-of-node-loop-for-phy-quantities

                    if(CalcType==FECalcType::COMPUTERESIDUAL||
                       CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_LocalShp.m_Test=t_FE.m_LineShp.shape_value(i);
                            m_LocalShp.m_GradTest=t_FE.m_LineShp.shape_grad(i);
                            m_LocalShp.m_Trial=0.0;
                            m_LocalShp.m_GradTrial=0.0;
                            iInd=cell.ElmtConn[i-1];
                            runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                      m_Normal,
                                      m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                      m_LocalK,m_LocalR);
                            // assemble local residual to global
                            assembleLocalResidual2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,JxW,t_DofHandler,m_LocalR,RHS);
                            if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                                for(j=1;j<=nNodesPerBCElmt;j++){
                                    m_LocalShp.m_Trial    =t_FE.m_LineShp.shape_value(j);
                                    m_LocalShp.m_GradTrial=t_FE.m_LineShp.shape_grad(j);
                                    jInd=cell.ElmtConn[j-1];

                                    runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                              m_Normal,
                                              m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                              m_LocalK,m_LocalR);

                                    // assemble local jacobian to global
                                    assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                                }
                            }
                        }
                    }
                    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_LocalShp.m_Test    =t_FE.m_LineShp.shape_value(i);
                            m_LocalShp.m_GradTest=t_FE.m_LineShp.shape_grad(i);
                            iInd=cell.ElmtConn[i-1];
                            for(j=1;j<=nNodesPerBCElmt;j++){
                                m_LocalShp.m_Trial    =t_FE.m_LineShp.shape_value(j);
                                m_LocalShp.m_GradTrial=t_FE.m_LineShp.shape_grad(j);
                                jInd=cell.ElmtConn[j-1];

                                runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                          m_Normal,
                                          m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                          m_LocalK,m_LocalR);

                                // assemble local jacobian to global
                                assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                            }
                        }
                    }// end-of-K-and-R-calculation
                }// end-of-line-qpoints-loop (1D case)

            } // end-of-dim=1-case
            else if(m_LocalElmtInfo.m_Dim==2){
                // for 'surf' element case
                for(gpInd=1;gpInd<=t_FE.m_SurfaceQpoints.getQPointsNum();gpInd++){
                    xi =t_FE.m_SurfaceQpoints.getIthPointJthCoord(gpInd,1);
                    eta=t_FE.m_SurfaceQpoints.getIthPointJthCoord(gpInd,2);
                    w  =t_FE.m_SurfaceQpoints.getIthPointJthCoord(gpInd,0);

                    t_FE.m_SurfaceShp.calc(xi,eta,0.0,cell.ElmtNodeCoords0,false);

                    // for normal vector calculation
                    m_XS.setToZeros();m_Normal=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_XS(1,1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*cell.ElmtNodeCoords0(i,1);
                        m_XS(1,2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*cell.ElmtNodeCoords0(i,1);

                        m_XS(2,1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*cell.ElmtNodeCoords0(i,2);
                        m_XS(2,2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*cell.ElmtNodeCoords0(i,2);

                        m_XS(3,1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*cell.ElmtNodeCoords0(i,3);
                        m_XS(3,2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*cell.ElmtNodeCoords0(i,3);
                    }
                    m_Normal(1)=m_XS(2,1)*m_XS(3,2)-m_XS(3,1)*m_XS(2,2);
                    m_Normal(2)=m_XS(3,1)*m_XS(1,2)-m_XS(1,1)*m_XS(3,2);
                    m_Normal(3)=m_XS(1,1)*m_XS(2,2)-m_XS(2,1)*m_XS(1,2);

                    dist=sqrt(m_Normal(1)*m_Normal(1)+m_Normal(2)*m_Normal(2)+m_Normal(3)*m_Normal(3));
                    m_Normal(1)=m_Normal(1)/dist;
                    m_Normal(2)=m_Normal(2)/dist;
                    m_Normal(3)=m_Normal(3)/dist;

                    t_FE.m_SurfaceShp.calc(xi,eta,0.0,cell.ElmtNodeCoords0,true);
                    JxW=t_FE.m_SurfaceShp.getJacDet()*w;

                    for(k=1;k<=m_LocalElmtInfo.m_DofsNum;k++){
                        m_LocalElmtSoln.m_QpU[k]=0.0;
                        m_LocalElmtSoln.m_QpUold[k]=0.0;
                        m_LocalElmtSoln.m_QpUolder[k]=0.0;
                        m_LocalElmtSoln.m_QpV[k]=0.0;

                        m_LocalElmtSoln.m_QpGradU[k]=0.0;
                        m_LocalElmtSoln.m_QpGradUold[k]=0.0;
                        m_LocalElmtSoln.m_QpGradUolder[k]=0.0;
                    }
                    m_LocalElmtInfo.m_QpCoords0=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        iInd=cell.ElmtConn[i-1];
                        for(k=1;k<=m_LocalElmtInfo.m_DofsNum;k++){
                            GlobalDofID=t_DofHandler.getIthNodeJthDofID(iInd,DofIDs[k-1]);

                            m_LocalElmtSoln.m_QpU[k]     +=t_FE.m_SurfaceShp.shape_value(i)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpUold[k]  +=t_FE.m_SurfaceShp.shape_value(i)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpUolder[k]+=t_FE.m_SurfaceShp.shape_value(i)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpV[k]     +=t_FE.m_SurfaceShp.shape_value(i)*V.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradU[k](1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradU[k](2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*U.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradU[k](3)+=t_FE.m_SurfaceShp.shape_grad(i)(3)*U.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradUold[k](1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUold[k](2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*Uold.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUold[k](3)+=t_FE.m_SurfaceShp.shape_grad(i)(3)*Uold.getIthValueFromGhost(GlobalDofID);

                            m_LocalElmtSoln.m_QpGradUolder[k](1)+=t_FE.m_SurfaceShp.shape_grad(i)(1)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUolder[k](2)+=t_FE.m_SurfaceShp.shape_grad(i)(2)*Uolder.getIthValueFromGhost(GlobalDofID);
                            m_LocalElmtSoln.m_QpGradUolder[k](3)+=t_FE.m_SurfaceShp.shape_grad(i)(3)*Uolder.getIthValueFromGhost(GlobalDofID);
                        }
                                
                        m_LocalElmtInfo.m_QpCoords0(1)+=t_FE.m_SurfaceShp.shape_value(i)*cell.ElmtNodeCoords0(i,1);
                        m_LocalElmtInfo.m_QpCoords0(2)+=t_FE.m_SurfaceShp.shape_value(i)*cell.ElmtNodeCoords0(i,2);
                        m_LocalElmtInfo.m_QpCoords0(3)+=t_FE.m_SurfaceShp.shape_value(i)*cell.ElmtNodeCoords0(i,3);
                    
                    }// end-of-node-loop-for-phy-quantities

                    if(CalcType==FECalcType::COMPUTERESIDUAL||
                       CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_LocalShp.m_Test    =t_FE.m_SurfaceShp.shape_value(i);
                            m_LocalShp.m_GradTest=t_FE.m_SurfaceShp.shape_grad(i);
                            m_LocalShp.m_Trial=0.0;
                            m_LocalShp.m_GradTrial=0.0;
                            iInd=cell.ElmtConn[i-1];
                            runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                      m_Normal,
                                      m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                      m_LocalK,m_LocalR);
                                    
                            // assemble local residual to global
                            assembleLocalResidual2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,JxW,t_DofHandler,m_LocalR,RHS);
                            if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
                                for(j=1;j<=nNodesPerBCElmt;j++){
                                    m_LocalShp.m_Trial    =t_FE.m_SurfaceShp.shape_value(j);
                                    m_LocalShp.m_GradTrial=t_FE.m_SurfaceShp.shape_grad(j);
                                    jInd=cell.ElmtConn[j-1];

                                    runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                              m_Normal,
                                              m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                              m_LocalK,m_LocalR);

                                    // assemble local jacobian to global
                                    assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                                }
                            }
                        }
                    }
                    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_LocalShp.m_Test    =t_FE.m_SurfaceShp.shape_value(i);
                            m_LocalShp.m_GradTest=t_FE.m_SurfaceShp.shape_grad(i);
                            iInd=cell.ElmtConn[i-1];
                            for(j=1;j<=nNodesPerBCElmt;j++){
                                m_LocalShp.m_Trial    =t_FE.m_SurfaceShp.shape_value(j);
                                m_LocalShp.m_GradTrial=t_FE.m_SurfaceShp.shape_grad(j);
                                jInd=cell.ElmtConn[j-1];

                                runBCLibs(CalcType,t_BCType,BCValue,Ctan,Parameters,
                                          m_Normal,
                                          m_LocalElmtInfo,m_LocalElmtSoln,m_LocalShp,
                                          m_LocalK,m_LocalR);

                                // assemble local jacobian to global
                                assembleLocalJacobian2Global(m_LocalElmtInfo.m_DofsNum,DofIDs,iInd,jInd,JxW,t_DofHandler,m_LocalK,AMATRIX);
                            }
                        }
                                
                    }// end-of-K-and-R-calculation

                }// end-of-qpoints-loop

            }// end-of-dim=2-case

        }// end-of-element-loop
    }// end-of-boundary-name-list-loop

    U.destroyGhostCopy();
    Uold.destroyGhostCopy();
    Uolder.destroyGhostCopy();
    V.destroyGhostCopy();


    if(CalcType==FECalcType::COMPUTERESIDUAL) RHS.assemble();
    if(CalcType==FECalcType::COMPUTEJACOBIAN) AMATRIX.assemble();
    if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        RHS.assemble();
        AMATRIX.assemble();
    }

}