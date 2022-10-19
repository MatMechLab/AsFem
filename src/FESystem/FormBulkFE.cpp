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
//+++ Date   : 2020.11.30
//+++ Purpose: Loop over all the bulk elements, calculate the system
//+++          residual, and K matrix
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/BulkFESystem.h"

void BulkFESystem::formBulkFE(const FECalcType &t_calctype,const double &t,const double &dt,const double (&ctan)[3],
                              Mesh &t_mesh,const DofHandler &t_dofhandler,FE &t_fe,
                              ElmtSystem &t_elmtsystem,MateSystem &t_matesystem,
                              SolutionSystem &t_solutionsystem,
                              SparseMatrix &AMATRIX,Vector &RHS){
    if(t_calctype==FECalcType::COMPUTERESIDUAL){
        RHS.setToZero();
    }
    else if(t_calctype==FECalcType::COMPUTEJACOBIAN){
        AMATRIX.setToZero();
    }

    // for the current and the previous steps' solution array
    t_solutionsystem.m_u_temp.makeGhostCopy();// we always use u_temp as u-current in FormBulkFE !!!
    t_solutionsystem.m_u_old.makeGhostCopy();
    t_solutionsystem.m_u_older.makeGhostCopy();

    // for the velocity and acceleration
    t_solutionsystem.m_v.makeGhostCopy();
    t_solutionsystem.m_a.makeGhostCopy();

    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);

    int rankne=t_mesh.getBulkMeshBulkElmtsNum()/m_size;
    int eStart=m_rank*rankne;
    int eEnd=(m_rank+1)*rankne;
    if(m_rank==m_size-1) eEnd=t_mesh.getBulkMeshBulkElmtsNum();

    int nDim,e,qpoints_num;
    int ndofs_per_elmt;
    double xi,eta,zeta,w,J,JxW,elVolume;
    nDim=t_mesh.getBulkMeshMaxDim();

    int subelmtid;
    int globaldofid,globalnodeid;
    int globalI,globalJ;

    m_local_elmtinfo.m_dim=nDim;
    m_local_elmtinfo.m_nodesnum=t_mesh.getBulkMeshNodesNumPerBulkElmt();
    m_local_elmtinfo.m_t=t;
    m_local_elmtinfo.m_dt=dt;
    

    for(int ee=eStart;ee<eEnd;ee++){
        e=ee+1;

        t_mesh.getBulkMeshIthBulkElmtNodeCoords0(e,m_nodes0);// for nodal coordinates in reference configuration
        t_mesh.getBulkMeshIthBulkElmtNodeCoords(e,m_nodes);// for nodal coordinates in current configuration

        t_mesh.getBulkMeshIthBulkElmtConnectivity(e,m_elmtconn);// for current element's connectivity

        t_dofhandler.getIthBulkElmtDofIDs(e,m_elmtdofsid);// the global dofs id, start from 1
        ndofs_per_elmt=t_dofhandler.getIthBulkElmtDofsNum(e);

        // for the local solution
        for(int i=0;i<ndofs_per_elmt;i++){
            m_elmtUolder[i]=t_solutionsystem.m_u_older.getIthValueFromGhost(m_elmtdofsid[i]);
            m_elmtUold[i]=t_solutionsystem.m_u_old.getIthValueFromGhost(m_elmtdofsid[i]);
            m_elmtU[i]=t_solutionsystem.m_u_temp.getIthValueFromGhost(m_elmtdofsid[i]);

            m_elmtV[i]=t_solutionsystem.m_v.getIthValueFromGhost(m_elmtdofsid[i]);
            m_elmtA[i]=t_solutionsystem.m_a.getIthValueFromGhost(m_elmtdofsid[i]);
        }

        //***********************************************************
        //*** now we do the gauss point integration(qpoints loop)
        //***********************************************************
        elVolume=0.0;
        qpoints_num=t_fe.m_bulk_qpoints.getQPointsNum();
        for(int qpInd=1;qpInd<=qpoints_num;qpInd++){
            w =t_fe.m_bulk_qpoints.getIthPointJthCoord(qpInd,0);
            xi=t_fe.m_bulk_qpoints.getIthPointJthCoord(qpInd,1);
            if(nDim==1){
                eta=0.0;zeta=0.0;
            }
            else if(nDim==2){
                eta=t_fe.m_bulk_qpoints.getIthPointJthCoord(qpInd,2);
                zeta=0.0;
            }
            else if(nDim==3){
                eta =t_fe.m_bulk_qpoints.getIthPointJthCoord(qpInd,2);
                zeta=t_fe.m_bulk_qpoints.getIthPointJthCoord(qpInd,3);
            }
            t_fe.m_bulk_shp.calc(xi,eta,zeta,m_nodes0,true);
            J=t_fe.m_bulk_shp.getJacDet();
            JxW=J*w;
            elVolume+=1.0*JxW;// for the 'volume' of current element

            m_local_elmtinfo.m_gpCoords0=0.0;
            for(int i=1;i<=m_bulkelmt_nodesnum;i++){
                m_local_elmtinfo.m_gpCoords0(1)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,1);
                m_local_elmtinfo.m_gpCoords0(2)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,2);
                m_local_elmtinfo.m_gpCoords0(3)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,3);
            }


            //*********************************************************************
            //*** loop over all the sub element/modulus of current bulk element
            //*********************************************************************
            if(t_calctype==FECalcType::COMPUTERESIDUAL){
                m_localR.setToZero();// its size is the maximum dofs per node
            }
            else if(t_calctype==FECalcType::COMPUTEJACOBIAN){
                m_localK.setToZero();// its size is the maximum dofs per node
            }

            if(t_calctype!=FECalcType::INITMATERIAL){
                // get the old material properties on each qpoint
                t_matesystem.m_materialcontainer_old.getScalarMaterialsRef()=t_solutionsystem.m_qpoints_scalarmaterials[(e-1)*qpoints_num+qpInd-1];
                t_matesystem.m_materialcontainer_old.getVectorMaterialsRef()=t_solutionsystem.m_qpoints_vectormaterials[(e-1)*qpoints_num+qpInd-1];
                t_matesystem.m_materialcontainer_old.getRank2MaterialsRef()=t_solutionsystem.m_qpoints_rank2materials[(e-1)*qpoints_num+qpInd-1];
                t_matesystem.m_materialcontainer_old.getRank4MaterialsRef()=t_solutionsystem.m_qpoints_rank4materials[(e-1)*qpoints_num+qpInd-1];
            }

            for(int subelmt=1;subelmt<=t_elmtsystem.getIthBulkElmtSubElmtsNum(e);subelmt++){
                if(t_calctype==FECalcType::COMPUTERESIDUAL){
                    m_subR.setToZero();
                }
                else if(t_calctype==FECalcType::COMPUTEJACOBIAN){
                    m_subK.setToZero();
                }

                subelmtid=t_elmtsystem.getIthBulkElmtJthSubElmtID(e,subelmt);
                m_subelmt_dofs=static_cast<int>(t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_dof_ids.size());
                m_local_elmtinfo.m_dofsnum=m_subelmt_dofs;
                
                for(int i=0;i<m_subelmt_dofs;i++){
                    m_subelmtdofsid[i]=t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_dof_ids[i];// start from 1

                    //*********************************************************************
                    //*** prepare physical quantities on each integration point
                    //*********************************************************************
                    // here, the U/V/A vector's index should start from 1, to make it consistent with
                    // 1st dof, 2nd dof, 3rd one , and so on. So, please let the index starts from 1 !!!
                    m_local_elmtsoln.m_gpU[i+1]=0.0;
                    m_local_elmtsoln.m_gpUold[i+1]=0.0;
                    m_local_elmtsoln.m_gpUolder[i+1]=0.0;

                    m_local_elmtsoln.m_gpV[i+1]=0.0;
                    m_local_elmtsoln.m_gpA[i+1]=0.0;

                    m_local_elmtsoln.m_gpGradU[i+1]=0.0;
                    m_local_elmtsoln.m_gpGradUold[i+1]=0.0;
                    m_local_elmtsoln.m_gpGradUolder[i+1]=0.0;

                    m_local_elmtsoln.m_gpGradV[i+1]=0.0;
                    for(int j=1;j<=m_bulkelmt_nodesnum;j++){
                        globalnodeid=t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,j);
                        globaldofid=t_dofhandler.getIthNodeJthDofID(globalnodeid,m_subelmtdofsid[i]);
                        m_local_elmtsoln.m_gpUolder[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solutionsystem.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpUold[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solutionsystem.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpU[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solutionsystem.m_u_temp.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpV[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solutionsystem.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpA[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solutionsystem.m_a.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpGradUolder[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solutionsystem.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUolder[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solutionsystem.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUolder[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solutionsystem.m_u_older.getIthValueFromGhost(globaldofid);
                        
                        m_local_elmtsoln.m_gpGradUold[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solutionsystem.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUold[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solutionsystem.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUold[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solutionsystem.m_u_old.getIthValueFromGhost(globaldofid);
                        
                        m_local_elmtsoln.m_gpGradU[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solutionsystem.m_u_temp.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradU[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solutionsystem.m_u_temp.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradU[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solutionsystem.m_u_temp.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpGradV[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solutionsystem.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradV[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solutionsystem.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradV[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solutionsystem.m_v.getIthValueFromGhost(globaldofid);
                    
                    }

                }// end-of-sub-element-dofs-loop

                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                if(t_calctype==FECalcType::INITMATERIAL){
                    t_matesystem.initBulkMateLibs(t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_matetype,
                                                  t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_json_params,
                                                  m_local_elmtinfo,m_local_elmtsoln);
                }
                else{
                    t_matesystem.runBulkMateLibs(t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_matetype,
                                                 t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_json_params,
                                                 m_local_elmtinfo,
                                                 m_local_elmtsoln);
                }
                if(t_calctype==FECalcType::COMPUTERESIDUAL){
                    for(int i=1;i<=m_bulkelmt_nodesnum;i++){
                        m_local_shp.m_test=t_fe.m_bulk_shp.shape_value(i);
                        m_local_shp.m_grad_test=t_fe.m_bulk_shp.shape_grad(i);
                        m_local_shp.m_trial=0.0;
                        m_local_shp.m_grad_trial=0.0;

                        t_elmtsystem.runBulkElmtLibs(t_calctype,ctan,subelmtid,
                                         t_matesystem.m_materialcontainer_old,
                                         t_matesystem.m_materialcontainer,
                                         m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                         m_subK,m_subR);

                        globalI=t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,i);
                        
                        assembleLocalResidual2GlobalR(m_subelmt_dofs,m_subelmtdofsid,globalI,t_dofhandler,JxW,m_subR,RHS);
                    
                    }
                }
                else if(t_calctype==FECalcType::COMPUTEJACOBIAN){
                    for(int i=1;i<=m_bulkelmt_nodesnum;i++){
                        m_local_shp.m_test=t_fe.m_bulk_shp.shape_value(i);
                        m_local_shp.m_grad_test=t_fe.m_bulk_shp.shape_grad(i);
                        globalI=t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,i);
                        for(int j=1;j<=m_bulkelmt_nodesnum;j++){
                            m_local_shp.m_trial=t_fe.m_bulk_shp.shape_value(j);
                            m_local_shp.m_grad_trial=t_fe.m_bulk_shp.shape_grad(j);
                            globalJ=t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,j);
                            
                            t_elmtsystem.runBulkElmtLibs(t_calctype,ctan,subelmtid,
                                         t_matesystem.m_materialcontainer_old,
                                         t_matesystem.m_materialcontainer,
                                         m_local_elmtinfo,m_local_elmtsoln,m_local_shp,
                                         m_subK,m_subR);

                            assembleLocalJacobian2GlobalK(m_subelmt_dofs,m_subelmtdofsid,globalI,globalJ,JxW,t_dofhandler,m_subK,AMATRIX);
                            
                        }
                    }
                }// end-of-residual-jacobian-calc-in-subElement                          

            }// end-of-sub-element-loop

            if(t_calctype==FECalcType::INITMATERIAL){
                // store the material properties on each qpoint
                t_solutionsystem.m_qpoints_scalarmaterials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getScalarMaterialsRef();
                t_solutionsystem.m_qpoints_vectormaterials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getVectorMaterialsRef();
                t_solutionsystem.m_qpoints_rank2materials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getRank2MaterialsRef();
                t_solutionsystem.m_qpoints_rank4materials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getRank4MaterialsRef();
            }
            else if(t_calctype==FECalcType::UPDATEMATERIAL){
                // store the material properties on each qpoint
                t_solutionsystem.m_qpoints_scalarmaterials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getScalarMaterialsRef();
                t_solutionsystem.m_qpoints_vectormaterials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getVectorMaterialsRef();
                t_solutionsystem.m_qpoints_rank2materials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getRank2MaterialsRef();
                t_solutionsystem.m_qpoints_rank4materials[(e-1)*qpoints_num+qpInd-1]=t_matesystem.m_materialcontainer.getRank4MaterialsRef();

                // update the old material
                t_solutionsystem.m_qpoints_scalarmaterials_old[(e-1)*qpoints_num+qpInd-1]=t_solutionsystem.m_qpoints_scalarmaterials[(e-1)*qpoints_num+qpInd-1];
                t_solutionsystem.m_qpoints_vectormaterials_old[(e-1)*qpoints_num+qpInd-1]=t_solutionsystem.m_qpoints_vectormaterials[(e-1)*qpoints_num+qpInd-1];
                t_solutionsystem.m_qpoints_rank2materials_old[(e-1)*qpoints_num+qpInd-1]=t_solutionsystem.m_qpoints_rank2materials[(e-1)*qpoints_num+qpInd-1];
                t_solutionsystem.m_qpoints_rank4materials_old[(e-1)*qpoints_num+qpInd-1]=t_solutionsystem.m_qpoints_rank4materials[(e-1)*qpoints_num+qpInd-1];
            }

        }// end-of-qpoints-loop

    }// end-of-element-loop

    // finish the final assemble
    if(t_calctype==FECalcType::COMPUTERESIDUAL){
        RHS.assemble();
    }
    else if(t_calctype==FECalcType::COMPUTEJACOBIAN){
        AMATRIX.assemble();
    }

    // release the ghost copy
    t_solutionsystem.m_u_temp.destroyGhostCopy();
    t_solutionsystem.m_u_old.destroyGhostCopy();
    t_solutionsystem.m_u_older.destroyGhostCopy();
    t_solutionsystem.m_v.destroyGhostCopy();
    t_solutionsystem.m_a.destroyGhostCopy();

}