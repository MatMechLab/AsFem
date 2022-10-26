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
//+++ Date   : 2022.09.28
//+++ Purpose: execute the side integral postprocess and return a scalar value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

double Postprocessor::executeSideIntegralPostprocess(const PostprocessorType &pps_type,
                                          const int &dofid,
                                          const vector<string> &sidenames,
                                          const nlohmann::json &t_parameters,
                                          const Mesh &t_mesh,
                                          const DofHandler &t_dofhandler,
                                          FE &t_fe,
                                          SolutionSystem &t_soln,
                                          ProjectionSystem &t_projsystem){
    double pps_value,side_area;
    double pps_value_global,side_area_global;
    string sidename;
    int i,j,iInd,e,rankne,eStart,eEnd,nElmts,nNodesPerBCElmt;
    int nqpoints;
    double dist;
    double xi,eta,w,JxW;

    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    t_soln.m_u_current.makeGhostCopy();
    t_soln.m_u_old.makeGhostCopy();
    t_soln.m_u_older.makeGhostCopy();
    t_soln.m_v.makeGhostCopy();

    pps_value=0.0;
    side_area=0.0;
    for(i=0;i<static_cast<int>(sidenames.size());i++){
        sidename=sidenames[i];
        nElmts=t_mesh.getBulkMeshElmtsNumViaPhyName(sidename);

        rankne=nElmts/m_size;
        eStart=m_rank*rankne;
        eEnd=(m_rank+1)*rankne;
        if(m_rank==m_size-1) eEnd=nElmts;

        m_local_elmtinfo.m_dim=t_mesh.getBulkMeshElmtDimViaPhyName(sidename);
        nNodesPerBCElmt=0;

        for(e=eStart;e<eEnd;e++){
            nNodesPerBCElmt=t_mesh.getBulkMeshIthElmtNodesNumViaPhyName(sidename,e+1);
            m_local_elmtinfo.m_nodesnum=nNodesPerBCElmt;
            JxW=0.0;
            if(m_local_elmtinfo.m_dim==0){
                // for 'point' case
                JxW=1.0;
                side_area=1.0;
                for(i=1;i<=nNodesPerBCElmt;i++){
                    j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(sidename,e+1,i);// global node id
                    iInd=t_dofhandler.getIthNodeJthDofID(j,dofid);
                    m_local_elmtinfo.m_gpCoords0(1)=t_mesh.getBulkMeshIthNodeJthCoord0(j,1);
                    m_local_elmtinfo.m_gpCoords0(2)=t_mesh.getBulkMeshIthNodeJthCoord0(j,2);
                    m_local_elmtinfo.m_gpCoords0(3)=t_mesh.getBulkMeshIthNodeJthCoord0(j,3);
                    m_local_shp.m_test=1.0;
                    m_local_shp.m_grad_test=0.0;
                    m_local_shp.m_trial=0.0;
                    m_local_shp.m_grad_trial=0.0;
                    pps_value+=JxW*runSideIntegralPostprocessLibs(pps_type,iInd,j,t_parameters,m_local_elmtinfo,m_local_shp,t_soln,t_projsystem);
                }
            }// end-of-dim=0-case
            else{
                t_mesh.getBulkMeshIthElmtNodeCoords0ViaPhyName(sidename,e+1,m_nodes0);

                nqpoints=0;
                if(m_local_elmtinfo.m_dim==1){
                    // for 'line' element case
                    nqpoints=t_fe.m_line_qpoints.getQPointsNum();
                }
                else if(m_local_elmtinfo.m_dim==2){
                    // for 'surface' element case
                    nqpoints=t_fe.m_surface_qpoints.getQPointsNum();
                }
                else{
                    MessagePrinter::printErrorTxt("dim>=3 is invalid for a side integral postprocess, please check your code or your input file");
                    MessagePrinter::exitAsFem();
                }

                // do the gauss point integration loop
                for(int gpInd=1;gpInd<=nqpoints;gpInd++){
                    if(m_local_elmtinfo.m_dim==1){
                        xi =t_fe.m_line_qpoints.getIthPointJthCoord(gpInd,1);
                        eta=0.0;
                        w  =t_fe.m_line_qpoints.getIthPointJthCoord(gpInd,0);
                        t_fe.m_line_shp.calc(xi,eta,0.0,m_nodes0,false);
                        // for the normal vector of current qpoint
                        m_xs.setToZeros();m_normal=0.0;
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_xs(1,1)+=t_fe.m_line_shp.shape_grad(i)(1)*m_nodes0(i,1);
                            m_xs(2,1)+=t_fe.m_line_shp.shape_grad(i)(1)*m_nodes0(i,2);
                        }
                        dist=sqrt(m_xs(1,1)*m_xs(1,1)+m_xs(2,1)*m_xs(2,1));
                        m_normal(1)= m_xs(2,1)/dist;// dy/dxi
                        m_normal(2)=-m_xs(1,1)/dist;// dx/dxi
                        m_normal(3)= 0.0;

                        t_fe.m_line_shp.calc(xi,eta,0.0,m_nodes0,true);
                        JxW=t_fe.m_line_shp.getJacDet()*w;
                    }
                    else if(m_local_elmtinfo.m_dim==2){
                        xi =t_fe.m_surface_qpoints.getIthPointJthCoord(gpInd,1);
                        eta=t_fe.m_surface_qpoints.getIthPointJthCoord(gpInd,2);
                        w  =t_fe.m_surface_qpoints.getIthPointJthCoord(gpInd,0);

                        t_fe.m_surface_shp.calc(xi,eta,0.0,m_nodes0,false);
                        // for normal vector calculation
                        m_xs.setToZeros();m_normal=0.0;
                        for(i=1;i<=nNodesPerBCElmt;i++){
                            m_xs(1,1)+=t_fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,1);
                            m_xs(1,2)+=t_fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,1);

                            m_xs(2,1)+=t_fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,2);
                            m_xs(2,2)+=t_fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,2);

                            m_xs(3,1)+=t_fe.m_surface_shp.shape_grad(i)(1)*m_nodes0(i,3);
                            m_xs(3,2)+=t_fe.m_surface_shp.shape_grad(i)(2)*m_nodes0(i,3);
                        }
                        m_normal(1)=m_xs(2,1)*m_xs(3,2)-m_xs(3,1)*m_xs(2,2);
                        m_normal(2)=m_xs(3,1)*m_xs(1,2)-m_xs(1,1)*m_xs(3,2);
                        m_normal(3)=m_xs(1,1)*m_xs(2,2)-m_xs(2,1)*m_xs(1,2);

                        dist=sqrt(m_normal(1)*m_normal(1)+m_normal(2)*m_normal(2)+m_normal(3)*m_normal(3));
                        m_normal(1)=m_normal(1)/dist;
                        m_normal(2)=m_normal(2)/dist;
                        m_normal(3)=m_normal(3)/dist;

                        t_fe.m_surface_shp.calc(xi,eta,0.0,m_nodes0,true);
                        JxW=t_fe.m_surface_shp.getJacDet()*w;
                    }

                    side_area+=1.0*JxW;

                    //********************************************************
                    //*** for the physical quantities on current qpoint
                    //********************************************************
                    
                    m_local_elmtinfo.m_gpCoords0=0.0;
                    for(i=1;i<=nNodesPerBCElmt;i++){
                        j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(sidename,e+1,i);//global id

                        if(m_local_elmtinfo.m_dim==1){
                            m_local_elmtinfo.m_gpCoords0(1)+=t_fe.m_line_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,1);
                            m_local_elmtinfo.m_gpCoords0(2)+=t_fe.m_line_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,2);
                            m_local_elmtinfo.m_gpCoords0(3)+=t_fe.m_line_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,3);
                        }
                        else if(m_local_elmtinfo.m_dim==2){   
                            m_local_elmtinfo.m_gpCoords0(1)+=t_fe.m_surface_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,1);
                            m_local_elmtinfo.m_gpCoords0(2)+=t_fe.m_surface_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,2);
                            m_local_elmtinfo.m_gpCoords0(3)+=t_fe.m_surface_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,3);
                        }   
                    }// end-of-node-loop-for-phy-quantities

                    for(i=1;i<=nNodesPerBCElmt;i++){
                        if(m_local_elmtinfo.m_dim==1){
                            m_local_shp.m_test=t_fe.m_line_shp.shape_value(i);
                            m_local_shp.m_grad_test=t_fe.m_line_shp.shape_grad(i);
                            m_local_shp.m_trial=0.0;
                            m_local_shp.m_grad_trial=0.0;
                        }
                        else{
                            m_local_shp.m_test=t_fe.m_surface_shp.shape_value(i);
                            m_local_shp.m_grad_test=t_fe.m_surface_shp.shape_grad(i);
                            m_local_shp.m_trial=0.0;
                            m_local_shp.m_grad_trial=0.0;
                        }
                        j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(sidename,e+1,i);//global id
                        if(dofid<1){
                            // this means no dof is given in the postprocess block, it is processing the material properties
                            iInd=t_dofhandler.getIthNodeJthDofID(j,1);
                        }
                        else{
                            iInd=t_dofhandler.getIthNodeJthDofID(j,dofid);
                        }

                        pps_value+=JxW*runSideIntegralPostprocessLibs(pps_type,iInd,j,t_parameters,m_local_elmtinfo,m_local_shp,t_soln,t_projsystem);
                    
                    }// end-of-node-loop-for-qpoint-quantities-accumulation

                }// end-of-qpoints-loop (1D+2D case)
            } // end-of-dim-selection
        }// end-of-element-loop

    }// end-of-side-name-loop

    t_soln.m_u_current.destroyGhostCopy();
    t_soln.m_u_old.destroyGhostCopy();
    t_soln.m_u_older.destroyGhostCopy();
    t_soln.m_v.destroyGhostCopy();

    // collect all the results from different cpus
    MPI_Allreduce(&pps_value,&pps_value_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&side_area,&side_area_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(pps_type==PostprocessorType::SIDEAVERAGEVALUE||
       pps_type==PostprocessorType::SIDEAVERAGESCALARMATERIALVALUE||
       pps_type==PostprocessorType::SIDEAVERAGEVECTORMATERIALVALUE||
       pps_type==PostprocessorType::SIDEAVERAGERANK2MATERIALVALUE||
       pps_type==PostprocessorType::SIDEAVERAGERANK4MATERIALVALUE){
        pps_value_global/=side_area_global;
    }
    
    return pps_value_global;
}