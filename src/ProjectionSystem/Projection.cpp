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
//+++ Date   : 2022.08.22
//+++ Purpose: execute the projection from integration points to
//+++          nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/ProjectionSystem.h"

void ProjectionSystem::executeProjection(const Mesh &t_mesh,const DofHandler &t_dofhandler,
                           const ElmtSystem &t_elmtsystem,
                           MateSystem &t_matesystem,
                           FE &t_fe,
                           SolutionSystem &t_solution,
                           const FEControlInfo &t_fectrlinfo){
    // if no projection is required, then return back
    if(getScalarMaterialNum()<1 &&
       getVectorMaterialNum()<1 &&
       getRank2MaterialNum()<1 &&
       getRank4MaterialNum()<1){
        return;
    }
    
    // before we do the projection, set all the vector to zeros
    m_data.m_proj_scalarmate_vec.setToZero();
    m_data.m_proj_vectormate_vec.setToZero();
    m_data.m_proj_rank2mate_vec.setToZero();
    m_data.m_proj_rank4mate_vec.setToZero();

    // for the current and the previous steps' solution array
    t_solution.m_u_current.makeGhostCopy();//
    t_solution.m_u_old.makeGhostCopy();
    t_solution.m_u_older.makeGhostCopy();

    // for the velocity and acceleration
    t_solution.m_v.makeGhostCopy();
    t_solution.m_a.makeGhostCopy();

    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);

    int rankne=t_mesh.getBulkMeshBulkElmtsNum()/m_size;
    int eStart=m_rank*rankne;
    int eEnd=(m_rank+1)*rankne;
    if(m_rank==m_size-1) eEnd=t_mesh.getBulkMeshBulkElmtsNum();

    int nDim,e,qpoints_num;
    double xi,eta,zeta,w,J,JxW;
    nDim=t_mesh.getBulkMeshMaxDim();

    int subelmtid;
    int globaldofid,globalnodeid;

    m_local_elmtinfo.m_dim=nDim;
    m_local_elmtinfo.m_nodesnum=t_mesh.getBulkMeshNodesNumPerBulkElmt();
    m_local_elmtinfo.m_t=t_fectrlinfo.t;
    m_local_elmtinfo.m_dt=t_fectrlinfo.dt;

    m_bulkelmt_nodesnum=t_mesh.getBulkMeshNodesNumPerBulkElmt();

    for(int ee=eStart;ee<eEnd;ee++){
        e=ee+1;

        t_mesh.getBulkMeshIthBulkElmtNodeCoords0(e,m_nodes0);// for nodal coordinates in reference configuration
        t_mesh.getBulkMeshIthBulkElmtNodeCoords(e,m_nodes);// for nodal coordinates in current configuration

        t_mesh.getBulkMeshIthBulkElmtConnectivity(e,m_elmtconn);// for current element's connectivity

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

            m_local_elmtinfo.m_gpCoords0=0.0;
            for(int i=1;i<=m_bulkelmt_nodesnum;i++){
                m_local_elmtinfo.m_gpCoords0(1)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,1);
                m_local_elmtinfo.m_gpCoords0(2)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,2);
                m_local_elmtinfo.m_gpCoords0(3)+=t_fe.m_bulk_shp.shape_value(i)*m_nodes0(i,3);
            }

            t_matesystem.m_materialcontainer_old.getScalarMaterialsRef()=t_solution.m_qpoints_scalarmaterials[(e-1)*qpoints_num+qpInd-1];
            t_matesystem.m_materialcontainer_old.getVectorMaterialsRef()=t_solution.m_qpoints_vectormaterials[(e-1)*qpoints_num+qpInd-1];
            t_matesystem.m_materialcontainer_old.getRank2MaterialsRef()=t_solution.m_qpoints_rank2materials[(e-1)*qpoints_num+qpInd-1];
            t_matesystem.m_materialcontainer_old.getRank4MaterialsRef()=t_solution.m_qpoints_rank4materials[(e-1)*qpoints_num+qpInd-1];

            for(int subelmt=1;subelmt<=t_elmtsystem.getIthBulkElmtSubElmtsNum(e);subelmt++){

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
                        m_local_elmtsoln.m_gpUolder[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solution.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpUold[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solution.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpU[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solution.m_u_current.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpV[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solution.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpA[i+1]+=t_fe.m_bulk_shp.shape_value(j)*t_solution.m_a.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpGradUolder[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solution.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUolder[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solution.m_u_older.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUolder[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solution.m_u_older.getIthValueFromGhost(globaldofid);
                        
                        m_local_elmtsoln.m_gpGradUold[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solution.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUold[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solution.m_u_old.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradUold[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solution.m_u_old.getIthValueFromGhost(globaldofid);
                        
                        m_local_elmtsoln.m_gpGradU[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solution.m_u_current.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradU[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solution.m_u_current.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradU[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solution.m_u_current.getIthValueFromGhost(globaldofid);

                        m_local_elmtsoln.m_gpGradV[i+1](1)+=t_fe.m_bulk_shp.shape_grad(j)(1)*t_solution.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradV[i+1](2)+=t_fe.m_bulk_shp.shape_grad(j)(2)*t_solution.m_v.getIthValueFromGhost(globaldofid);
                        m_local_elmtsoln.m_gpGradV[i+1](3)+=t_fe.m_bulk_shp.shape_grad(j)(3)*t_solution.m_v.getIthValueFromGhost(globaldofid);
                    
                    }

                }// end-of-sub-element-dofs-loop

                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                t_matesystem.runBulkMateLibs(t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_matetype,
                                             t_elmtsystem.getIthBulkElmtBlock(subelmtid).m_json_params,
                                             m_local_elmtinfo,
                                             m_local_elmtsoln);                        

            }// end-of-sub-element-loop
            // true--> for local projection
            runProjectionLibs(true,t_mesh,m_bulkelmt_nodesnum,m_elmtconn,JxW,t_fe.m_bulk_shp,t_matesystem.m_materialcontainer,m_data);

        }// end-of-qpoints-loop

    }// end-of-element libs

    t_solution.m_u_current.destroyGhostCopy();//
    t_solution.m_u_old.destroyGhostCopy();
    t_solution.m_u_older.destroyGhostCopy();

    // for the velocity and acceleration
    t_solution.m_v.destroyGhostCopy();
    t_solution.m_a.destroyGhostCopy();


    //******************************************************
    //*** after the local projection(in the element loop)
    //*** one needs to execute the global projection behavior
    //******************************************************
    // false--> for global projection action
    runProjectionLibs(false,t_mesh,m_bulkelmt_nodesnum,m_elmtconn,JxW,t_fe.m_bulk_shp,t_matesystem.m_materialcontainer,m_data);

}