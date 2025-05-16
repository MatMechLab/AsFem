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
//+++ Date   : 2022.09.28
//+++ Purpose: execute the side integral postprocess and return a scalar value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

double Postprocessor::executeVolumeIntegralPostprocess(const PostprocessorType &pps_type,
                                            const int &dofid,
                                            const vector<string> &domainnames,
                                            const nlohmann::json &t_parameters,
                                            const FECell &t_fecell,
                                            const DofHandler &t_dofhandler,
                                            FE &t_fe,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem){
    double pps_value,domain_volume;
    double pps_value_global,domain_volume_global;
    string domainname;
    int i,j,iInd,e,rankne,eStart,eEnd,nElmts,nNodesPerElmt;
    int nqpoints;
    double xi,eta,zeta,w,JxW;

    if(t_fecell.getFECellBulkElmtsNum()){}

    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    t_soln.m_Ucurrent.makeGhostCopy();
    t_soln.m_Uold.makeGhostCopy();
    t_soln.m_Uolder.makeGhostCopy();
    t_soln.m_V.makeGhostCopy();

    pps_value=0.0;
    domain_volume=0.0;
    for(i=0;i<static_cast<int>(domainnames.size());i++){
        domainname=domainnames[i];
        nElmts=2.0;//t_mesh.getBulkMeshElmtsNumViaPhyName(domainname);
        
        rankne=nElmts/m_size;
        eStart=m_rank*rankne;
        eEnd=(m_rank+1)*rankne;
        if(m_rank==m_size-1) eEnd=nElmts;

        // m_local_elmtinfo.m_dim=t_mesh.getBulkMeshElmtDimViaPhyName(domainname);
        nNodesPerElmt=0;

        for(e=eStart;e<eEnd;e++){
            // nNodesPerElmt=t_mesh.getBulkMeshIthElmtNodesNumViaPhyName(domainname,e+1);
            m_local_elmtinfo.m_NodesNum=nNodesPerElmt;
            JxW=0.0;
            if(m_local_elmtinfo.m_Dim<=0){
                MessagePrinter::printErrorTxt("Invalid dim(<=0) for volume integral postprocess, please check your input file or your code");
                MessagePrinter::exitAsFem();
            }// end-of-dim=0-case
            else{
                // t_mesh.getBulkMeshIthElmtNodeCoords0ViaPhyName(domainname,e+1,m_nodes0);

                nqpoints=t_fe.m_BulkQpoints.getQPointsNum();

                // do the gauss point integration loop
                for(int gpInd=1;gpInd<=nqpoints;gpInd++){
                    xi=t_fe.m_BulkQpoints.getIthPointJthCoord(gpInd,1);
                    if(m_local_elmtinfo.m_Dim==2){
                        eta=t_fe.m_BulkQpoints.getIthPointJthCoord(gpInd,2);
                    }
                    else if(m_local_elmtinfo.m_Dim==3){
                        eta =t_fe.m_BulkQpoints.getIthPointJthCoord(gpInd,2);
                        zeta=t_fe.m_BulkQpoints.getIthPointJthCoord(gpInd,3);
                    }
                    w  =t_fe.m_BulkQpoints.getIthPointJthCoord(gpInd,0);
                    t_fe.m_BulkShp.calc(xi,eta,zeta,m_nodes0,true);
                    JxW=w*t_fe.m_BulkShp.getJacDet();
                    domain_volume+=1.0*JxW;

                    //********************************************************
                    //*** for the physical quantities on current qpoint
                    //********************************************************
                    
                    m_local_elmtinfo.m_QpCoords0=0.0;
                    for(i=1;i<=nNodesPerElmt;i++){
                        // j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(domainname,e+1,i);//global id
                
                        // m_local_elmtinfo.m_gpCoords0(1)+=t_fe.m_bulk_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,1);
                        // m_local_elmtinfo.m_gpCoords0(2)+=t_fe.m_bulk_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,2);
                        // m_local_elmtinfo.m_gpCoords0(3)+=t_fe.m_bulk_shp.shape_value(i)*t_mesh.getBulkMeshIthNodeJthCoord(j,3);
                       
                    }// end-of-node-loop-for-phy-quantities

                    for(i=1;i<=nNodesPerElmt;i++){
                        m_local_shp.m_Test    =t_fe.m_BulkShp.shape_value(i);
                        m_local_shp.m_GradTest=t_fe.m_BulkShp.shape_grad(i);
                        m_local_shp.m_Trial    =0.0;
                        m_local_shp.m_GradTrial=0.0;
                        // j=t_mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(domainname,e+1,i);//global id
                        if(dofid<1){
                            // if no dofid is given, then we use the first one
                            iInd=t_dofhandler.getIthNodeJthDofID(j,1);
                        }
                        else{
                            iInd=t_dofhandler.getIthNodeJthDofID(j,dofid);
                        }

                        pps_value+=JxW*runVolumeIntegralPostprocessLibs(pps_type,iInd,j,t_parameters,m_local_elmtinfo,m_local_shp,t_soln,t_projsystem);
                    
                    }// end-of-node-loop-for-qpoint-quantities-accumulation

                }// end-of-qpoints-loop (1D+2D case)
            } // end-of-dim-selection
        }// end-of-element-loop
    }// end-of-side-name-loop


    t_soln.m_Ucurrent.destroyGhostCopy();
    t_soln.m_Uold.destroyGhostCopy();
    t_soln.m_Uolder.destroyGhostCopy();
    t_soln.m_V.destroyGhostCopy();

    // collect all the values from all the cpus
    MPI_Allreduce(&pps_value,&pps_value_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&domain_volume,&domain_volume_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    if(pps_type==PostprocessorType::VOLUMEAVERAGEVALUE||
       pps_type==PostprocessorType::VOLUMEAVERAGESCALARMATERIALVALUE||
       pps_type==PostprocessorType::VOLUMEAVERAGEVECTORMATERIALVALUE||
       pps_type==PostprocessorType::VOLUMEAVERAGERANK2MATERIALVALUE||
       pps_type==PostprocessorType::VOLUMEAVERAGERANK4MATERIALVALUE){
        pps_value_global/=domain_volume_global;
    }
    return pps_value_global;
}