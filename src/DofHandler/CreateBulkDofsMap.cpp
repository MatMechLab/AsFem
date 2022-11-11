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
//+++ Date   : 2022.07.23
//+++ Purpose: generate the dofs map for bulk elements
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

void BulkDofHandler::createBulkDofsMap(const Mesh &t_mesh,const ElmtSystem &t_elmtSystem){
    m_bulkelmts=t_mesh.getBulkMeshBulkElmtsNum();
    m_nodes=t_mesh.getBulkMeshNodesNum();
    m_total_dofs=m_nodes*m_maxdofs_pernode;// maxdofs_pernode dependes on the dofs name 
                                           // one defined in the input file
    m_maxdofs_perelmt=m_maxdofs_pernode*t_mesh.getBulkMeshNodesNumPerBulkElmt();

    // allocate memory for nodal and elemental dofs map
    m_nodal_dofids.resize(m_nodes,vector<int>(m_maxdofs_pernode,0));
    m_elmt_dofids.resize(m_bulkelmts,vector<int>(m_maxdofs_perelmt,0));


    m_active_dofs=0;
    int dofid,elmtid,nodeid;
    for(const auto &block:t_elmtSystem.getBulkElmtBlockList()){
        for(const auto &name:block.m_domain_namelist){
            for(int e=1;e<=t_mesh.getBulkMeshBulkElmtsNumViaPhyName(name);e++){
                elmtid=t_mesh.getBulkMeshIthBulkElmtIDViaPhyName(name,e);//global element id
                for(int i=1;i<=t_mesh.getBulkMeshIthBulkElmtNodesNum(elmtid);i++){
                    nodeid=t_mesh.getBulkMeshIthBulkElmtJthNodeID(elmtid,i);
                    for(int k=0;k<static_cast<int>(block.m_dof_ids.size());k++){
                        dofid=block.m_dof_ids[k];
                        if(m_nodal_dofids[nodeid-1][dofid-1]==0){
                            // if current nodal dof flag is zero, then we assign a non-zero value to it
                            // otherwise, it is duplicated, then skip it.
                            // Be careful about the ordering, the node id is not 1 to 1 mapping to the dof id !!!
                            m_active_dofs+=1;
                            m_nodal_dofids[nodeid-1][dofid-1]=m_active_dofs;
                        }
                    }
                }
            }
        }
    }

    // check the dofs status for each node
    bool HasDofID=false;
    for(int i=0;i<m_nodes;i++){
        HasDofID=false;
        for(int j=0;j<m_maxdofs_pernode;j++){
            if(m_nodal_dofids[i][j]){
                HasDofID=true;break;
            }
        }
        if(!HasDofID){
            MessagePrinter::printErrorTxt("Node-"+to_string(i+1)+" hasen\'t been assigned by the dof, please check your code");
            MessagePrinter::exitAsFem();
        }
    }

    // now we can create the elemental dofs map
    vector<int> row_maxnnz;/*< for the maximum nonzeros of each row */
    m_maxnnz=0;
    row_maxnnz.resize(m_active_dofs,0);
    for(int e=1;e<=m_bulkelmts;e++){
        for(int j=1;j<=t_mesh.getBulkMeshIthBulkElmtNodesNum(e);j++){
            nodeid=t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,j);
            for(int k=1;k<=m_maxdofs_pernode;k++){
                dofid=m_nodal_dofids[nodeid-1][k-1];
                m_elmt_dofids[e-1][(j-1)*m_maxdofs_pernode+k-1]=dofid;
                if(dofid!=0){
                    row_maxnnz[dofid-1]+=m_maxdofs_perelmt;
                    if(row_maxnnz[dofid-1]>m_maxnnz) m_maxnnz=row_maxnnz[dofid-1];
                }
            }
        }
        // remove the zero elements, only keep the nonzero part
        m_elmt_dofids[e-1].erase(
            remove(m_elmt_dofids[e-1].begin(),m_elmt_dofids[e-1].end(),0),
            m_elmt_dofids[e-1].end()
        );
        m_elmt_dofids[e-1].shrink_to_fit();
    }
    row_maxnnz.clear();

}