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
//+++ Date   : 2020.12.26
//+++ Purpose: apply the integrated boundary conditions with 'surface'
//+++          integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyIntegratedBC(const FECalcType &calctype,
                                 const double &bcvalue,
                                 const BCType &bctype,
                                 const double (&ctan)[3],
                                 const nlohmann::json &parameters,
                                 const vector<int> &dofids,
                                 const vector<string> &bcnamelist,
                                 const Mesh &mesh,const DofHandler &dofhandler,
                                 FE &fe,
                                 Vector &U,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS){
    // for other types of boundary conditions
    // for other type boundary conditions
    int rankne,eStart,eEnd,nNodesPerBCElmt;
    int e,i,j,k,iInd,jInd,gpInd;
    double xi,eta,w,JxW,dist;

    m_local_elmtinfo.m_dofsnum=static_cast<int>(dofids.size());

    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    U.makeGhostCopy();
    Uold.makeGhostCopy();
    Uolder.makeGhostCopy();
    V.makeGhostCopy();    

    for(const auto &name:bcnamelist){
        rankne=mesh.getBulkMeshElmtsNumViaPhyName(name)/m_size;
        eStart=m_rank*rankne;
        eEnd=(m_rank+1)*rankne;
        if(m_rank==m_size-1) eEnd=mesh.getBulkMeshElmtsNumViaPhyName(name);
        m_local_elmtinfo.m_dim=mesh.getBulkMeshElmtDimViaPhyName(name);
        nNodesPerBCElmt=0;

        if(m_local_elmtinfo.m_dim>2 || m_local_elmtinfo.m_dim<0){
            MessagePrinter::printErrorTxt("Invalid dim (="+to_string(m_local_elmtinfo.m_dim)+") for integrated boundary condition");
            MessagePrinter::exitAsFem();
        }

        for(e=eStart;e<eEnd;e++){
            nNodesPerBCElmt=mesh.getBulkMeshIthElmtNodesNumViaPhyName(name,e+1);      
            if(m_local_elmtinfo.m_dim==0){
                 // for 'point' case
                for(i=1;i<=nNodesPerBCElmt;i++){
                    j=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);// global id
                    m_local_elmtinfo.m_gpCoords0(1)=mesh.getBulkMeshIthNodeJthCoord0(j,1);
                    m_local_elmtinfo.m_gpCoords0(2)=mesh.getBulkMeshIthNodeJthCoord0(j,2);
                    m_local_elmtinfo.m_gpCoords0(3)=mesh.getBulkMeshIthNodeJthCoord0(j,3);
                    for(k=0;k<m_local_elmtinfo.m_dofsnum;k++){
                        iInd=dofhandler.getIthNodeJthDofID(j,k+1);
                        // get the local solutions
                        m_local_elmtsoln.m_gpU[k+1]=U.getIthValueFromGhost(iInd);
                        m_local_elmtsoln.m_gpUold[k+1]=Uold.getIthValueFromGhost(iInd);
                        m_local_elmtsoln.m_gpUolder[k+1]=Uolder.getIthValueFromGhost(iInd);
                        m_local_elmtsoln.m_gpV[k+1]=V.getIthValueFromGhost(iInd);
                    }
                }

                JxW=1.0;

                // do the residual and jacobian calculation
                if(calctype==FECalcType::COMPUTERESIDUAL){
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_local_shp.m_test=1.0;
                        m_local_shp.m_grad_test=0.0;
                        m_local_shp.m_trial=0.0;
                        m_local_shp.m_grad_trial=0.0;
                        iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                        runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                  m_normal,
                                  m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                  m_localK,m_localR);
                        // for assemble
                        assembleLocalResidual2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,JxW,dofhandler,m_localR,RHS);
                    }
                }
                else if(calctype==FECalcType::COMPUTEJACOBIAN){
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_local_shp.m_test=1.0;
                        m_local_shp.m_grad_test=0.0;
                        iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                        for(j=1;j<=nNodesPerBCElmt;j++){
                            m_local_shp.m_trial=0.0;
                            m_local_shp.m_grad_trial=0.0;
                            jInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,j);
                            runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                      m_normal,
                                      m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                      m_localK,m_localR);
                            // for jacobian assemble
                            assembleLocalJacobian2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,jInd,JxW,dofhandler,m_localK,AMATRIX);
                        }
                    }
                }// end-of-K-and-R-calculation

            }// end-of-dim=0-case
            else if(m_local_elmtinfo.m_dim==1){
                // for 'line' element case
                mesh.getBulkMeshIthElmtNodeCoords0ViaPhyName(name,e+1,m_nodes0);
                
                // do the gauss point integration loop
                for(gpInd=1;gpInd<=fe.m_line_qpoints.getQPointsNum();gpInd++){
                    xi =fe.m_line_qpoints.getIthPointJthCoord(gpInd,1);
                    eta=0.0;
                    w  =fe.m_line_qpoints.getIthPointJthCoord(gpInd,0);
                    fe.m_line_shp.calc(xi,eta,0.0,m_nodes0,false);

                    // for the normal vector of current qpoint
                    m_xs.setToZeros();m_normal=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_xs(1,1)+=fe.m_line_shp.shape_grad(i)(1)*m_nodes0(i,1);
                        m_xs(2,1)+=fe.m_line_shp.shape_grad(i)(1)*m_nodes0(i,2);
                    }
                    dist=sqrt(m_xs(1,1)*m_xs(1,1)+m_xs(2,1)*m_xs(2,1));
                    m_normal(1)= m_xs(2,1)/dist;// dy/dxi
                    m_normal(2)=-m_xs(1,1)/dist;// dx/dxi
                    m_normal(3)= 0.0;

                    //********************************************************
                    //*** for the physical quantities on current qpoint
                    //********************************************************
                    fe.m_line_shp.calc(xi,eta,0.0,m_nodes0,true);
                    JxW=fe.m_line_shp.getJacDet()*w;

                    for(k=1;k<=m_local_elmtinfo.m_dofsnum;k++){
                        m_local_elmtsoln.m_gpU[k]=0.0;
                        m_local_elmtsoln.m_gpUold[k]=0.0;
                        m_local_elmtsoln.m_gpUolder[k]=0.0;
                        m_local_elmtsoln.m_gpV[k]=0.0;

                        m_local_elmtsoln.m_gpGradU[k]=0.0;
                        m_local_elmtsoln.m_gpGradUold[k]=0.0;
                        m_local_elmtsoln.m_gpGradUolder[k]=0.0;
                    }
                    m_local_elmtinfo.m_gpCoords0=0.0;
                            
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        j=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);//global id
                        for(k=1;k<=m_local_elmtinfo.m_dofsnum;k++){
                            iInd=dofhandler.getIthNodeJthDofID(j,dofids[k-1]);

                            m_local_elmtsoln.m_gpU[k]+=fe.m_line_shp.shape_value(i)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpUold[k]+=fe.m_line_shp.shape_value(i)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpUolder[k]+=fe.m_line_shp.shape_value(i)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpV[k]+=fe.m_line_shp.shape_value(i)*V.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradU[k](1)+=fe.m_line_shp.shape_grad(i)(1)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradU[k](2)+=fe.m_line_shp.shape_grad(i)(2)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradU[k](3)+=fe.m_line_shp.shape_grad(i)(3)*U.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradUold[k](1)+=fe.m_line_shp.shape_grad(i)(1)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUold[k](2)+=fe.m_line_shp.shape_grad(i)(2)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUold[k](3)+=fe.m_line_shp.shape_grad(i)(3)*Uold.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradUolder[k](1)+=fe.m_line_shp.shape_grad(i)(1)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUolder[k](2)+=fe.m_line_shp.shape_grad(i)(2)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUolder[k](3)+=fe.m_line_shp.shape_grad(i)(3)*Uolder.getIthValueFromGhost(iInd);
                        }
                                
                        m_local_elmtinfo.m_gpCoords0(1)+=fe.m_line_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,1);
                        m_local_elmtinfo.m_gpCoords0(2)+=fe.m_line_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,2);
                        m_local_elmtinfo.m_gpCoords0(3)+=fe.m_line_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,3);
                            
                    }// end-of-node-loop-for-phy-quantities

                    if(calctype==FECalcType::COMPUTERESIDUAL){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_local_shp.m_test=fe.m_line_shp.shape_value(i);
                            m_local_shp.m_grad_test=fe.m_line_shp.shape_grad(i);
                            m_local_shp.m_trial=0.0;
                            m_local_shp.m_grad_trial=0.0;
                            iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                            runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                      m_normal,
                                      m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                      m_localK,m_localR);
                            // assemble local residual to global
                            assembleLocalResidual2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,JxW,dofhandler,m_localR,RHS);
                        }
                    }
                    else if(calctype==FECalcType::COMPUTEJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_local_shp.m_test=fe.m_line_shp.shape_value(i);
                            m_local_shp.m_grad_test=fe.m_line_shp.shape_grad(i);
                            iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                            for(j=1;j<=nNodesPerBCElmt;j++){
                                m_local_shp.m_trial=fe.m_line_shp.shape_value(j);
                                m_local_shp.m_grad_trial=fe.m_line_shp.shape_grad(j);
                                jInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,j);

                                runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                          m_normal,
                                          m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                          m_localK,m_localR);

                                // assemble local jacobian to global
                                assembleLocalJacobian2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,jInd,JxW,dofhandler,m_localK,AMATRIX);
                            }
                        }
                    }// end-of-K-and-R-calculation
                }// end-of-line-qpoints-loop (1D case)

            } // end-of-dim=1-case
            else if(m_local_elmtinfo.m_dim==2){
                // for 'surf' element case
                mesh.getBulkMeshIthElmtNodeCoords0ViaPhyName(name,e+1,m_nodes0);

                for(gpInd=1;gpInd<=fe.m_surface_qpoints.getQPointsNum();gpInd++){
                    xi =fe.m_surface_qpoints.getIthPointJthCoord(gpInd,1);
                    eta=fe.m_surface_qpoints.getIthPointJthCoord(gpInd,2);
                    w  =fe.m_surface_qpoints.getIthPointJthCoord(gpInd,0);

                    fe.m_surface_shp.calc(xi,eta,0.0,m_nodes0,false);

                    // for normal vector calculation
                    m_xs.setToZeros();m_normal=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        m_xs(1,1)+=fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,1);
                        m_xs(1,2)+=fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,1);

                        m_xs(2,1)+=fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,2);
                        m_xs(2,2)+=fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,2);

                        m_xs(3,1)+=fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,3);
                        m_xs(3,2)+=fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,3);
                    }
                    m_normal(1)=m_xs(2,1)*m_xs(3,2)-m_xs(3,1)*m_xs(2,2);
                    m_normal(2)=m_xs(3,1)*m_xs(1,2)-m_xs(1,1)*m_xs(3,2);
                    m_normal(3)=m_xs(1,1)*m_xs(2,2)-m_xs(2,1)*m_xs(1,2);

                    dist=sqrt(m_normal(1)*m_normal(1)+m_normal(2)*m_normal(2)+m_normal(3)*m_normal(3));
                    m_normal(1)=m_normal(1)/dist;
                    m_normal(2)=m_normal(2)/dist;
                    m_normal(3)=m_normal(3)/dist;

                    fe.m_surface_shp.calc(xi,eta,0.0,m_nodes0,true);
                    JxW=fe.m_surface_shp.getJacDet()*w;

                    for(k=1;k<=m_local_elmtinfo.m_dofsnum;k++){
                        m_local_elmtsoln.m_gpU[k]=0.0;
                        m_local_elmtsoln.m_gpUold[k]=0.0;
                        m_local_elmtsoln.m_gpUolder[k]=0.0;
                        m_local_elmtsoln.m_gpV[k]=0.0;

                        m_local_elmtsoln.m_gpGradU[k]=0.0;
                        m_local_elmtsoln.m_gpGradUold[k]=0.0;
                        m_local_elmtsoln.m_gpGradUolder[k]=0.0;
                    }
                    m_local_elmtinfo.m_gpCoords0=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        j=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);//global id
                        for(k=1;k<=m_local_elmtinfo.m_dofsnum;k++){
                            iInd=dofhandler.getIthNodeJthDofID(j,dofids[k-1]);

                            m_local_elmtsoln.m_gpU[k]+=fe.m_surface_shp.shape_value(i)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpUold[k]+=fe.m_surface_shp.shape_value(i)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpUolder[k]+=fe.m_surface_shp.shape_value(i)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpV[k]+=fe.m_surface_shp.shape_value(i)*V.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradU[k](1)+=fe.m_surface_shp.shape_grad(i)(1)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradU[k](2)+=fe.m_surface_shp.shape_grad(i)(2)*U.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradU[k](3)+=fe.m_surface_shp.shape_grad(i)(3)*U.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradUold[k](1)+=fe.m_surface_shp.shape_grad(i)(1)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUold[k](2)+=fe.m_surface_shp.shape_grad(i)(2)*Uold.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUold[k](3)+=fe.m_surface_shp.shape_grad(i)(3)*Uold.getIthValueFromGhost(iInd);

                            m_local_elmtsoln.m_gpGradUolder[k](1)+=fe.m_surface_shp.shape_grad(i)(1)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUolder[k](2)+=fe.m_surface_shp.shape_grad(i)(2)*Uolder.getIthValueFromGhost(iInd);
                            m_local_elmtsoln.m_gpGradUolder[k](3)+=fe.m_surface_shp.shape_grad(i)(3)*Uolder.getIthValueFromGhost(iInd);
                        }
                                
                        m_local_elmtinfo.m_gpCoords0(1)+=fe.m_surface_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,1);
                        m_local_elmtinfo.m_gpCoords0(2)+=fe.m_surface_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,2);
                        m_local_elmtinfo.m_gpCoords0(3)+=fe.m_surface_shp.shape_value(i)*mesh.getBulkMeshIthNodeJthCoord(j,3);
                    
                    }// end-of-node-loop-for-phy-quantities

                    if(calctype==FECalcType::COMPUTERESIDUAL){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_local_shp.m_test=fe.m_surface_shp.shape_value(i);
                            m_local_shp.m_grad_test=fe.m_surface_shp.shape_grad(i);
                            m_local_shp.m_trial=0.0;
                            m_local_shp.m_grad_trial=0.0;
                            iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                            runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                      m_normal,
                                      m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                      m_localK,m_localR);
                                    
                            // assemble local residual to global
                            assembleLocalResidual2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,JxW,dofhandler,m_localR,RHS);
                        }
                    }
                    else if(calctype==FECalcType::COMPUTEJACOBIAN){
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_local_shp.m_test=fe.m_surface_shp.shape_value(i);
                            m_local_shp.m_grad_test=fe.m_surface_shp.shape_grad(i);
                            iInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                            for(j=1;j<=nNodesPerBCElmt;j++){
                                m_local_shp.m_trial=fe.m_surface_shp.shape_value(j);
                                m_local_shp.m_grad_trial=fe.m_surface_shp.shape_grad(j);
                                jInd=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,j);

                                runBCLibs(calctype,bctype,bcvalue,ctan,parameters,
                                          m_normal,
                                          m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                          m_localK,m_localR);

                                // assemble local jacobian to global
                                assembleLocalJacobian2Global(m_local_elmtinfo.m_dofsnum,dofids,iInd,jInd,JxW,dofhandler,m_localK,AMATRIX);
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


    if(calctype==FECalcType::COMPUTERESIDUAL) RHS.assemble();
    if(calctype==FECalcType::COMPUTEJACOBIAN) AMATRIX.assemble();

}