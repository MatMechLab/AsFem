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
//+++ Date   : 2022.04.20
//+++ Purpose: the lagrange 1d mesh generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Lagrange1DMeshGenerator.h"

#include "Utils/MessagePrinter.h"

Lagrange1DMeshGenerator::Lagrange1DMeshGenerator(){
    m_mesh_generated=false;
}
Lagrange1DMeshGenerator::~Lagrange1DMeshGenerator(){
    leftconn.clear();
    rightconn.clear();
    leftnodes.clear();
    rightnodes.clear();
}

bool Lagrange1DMeshGenerator::generateMesh(const MeshType &t_meshtype,MeshData &t_meshdata){
    t_meshdata.m_maxdim=1;
    t_meshdata.m_mindim=0;
    t_meshdata.m_bulkelmt_connectivity.clear();
    t_meshdata.m_bulkelmt_volume.clear();

    t_meshdata.m_bulkelmt_type=t_meshtype;

    t_meshdata.m_surfaceelmt_connectivity.clear();
    t_meshdata.m_surfaceelmt_type=MeshType::NULLTYPE;
    t_meshdata.m_surfaceelmt_volume.clear();
    t_meshdata.m_surfaceelmt_vtktype=0;
    t_meshdata.m_surfaceelmts=0;

    t_meshdata.m_lineelmt_connectivity.clear();
    t_meshdata.m_lineelmt_type=MeshType::NULLTYPE;
    t_meshdata.m_lineelmt_volume.clear();
    t_meshdata.m_lineelmt_vtktype=0;
    t_meshdata.m_lineelmts=0;

    t_meshdata.m_nodesperlineelmt=0;
    t_meshdata.m_nodespersurfaceelmt=0;

    m_mesh_generated=false;

    if(t_meshtype==MeshType::EDGE2){
        t_meshdata.m_bulkelmt_vtktype=3;
        t_meshdata.m_nodesperbulkelmt=2;
        t_meshdata.m_order=1;
        t_meshdata.m_bulkelmt_typename="edge2";
    }
    else if(t_meshtype==MeshType::EDGE3){
        t_meshdata.m_bulkelmt_vtktype=4;
        t_meshdata.m_nodesperbulkelmt=3;
        t_meshdata.m_order=2;
        t_meshdata.m_bulkelmt_typename="edge3";
    }
    else if(t_meshtype==MeshType::EDGE4){
        t_meshdata.m_bulkelmt_vtktype=4;
        t_meshdata.m_nodesperbulkelmt=4;
        t_meshdata.m_order=3;
        t_meshdata.m_bulkelmt_typename="edge4";
    }
    else if(t_meshtype==MeshType::EDGE5){
        t_meshdata.m_bulkelmt_vtktype=4;
        t_meshdata.m_nodesperbulkelmt=5;
        t_meshdata.m_order=4;
        t_meshdata.m_bulkelmt_typename="edge5";
    }
    else{
        MessagePrinter::printErrorTxt("unsupported mesh type for 1D case. Currently, only edge2,edge3,edge4 and edge5 are supported.");
        m_mesh_generated=false;
        return m_mesh_generated;
    }

    // for bulk elements generation
    t_meshdata.m_elements=t_meshdata.m_nx+2;
    t_meshdata.m_bulkelmts=t_meshdata.m_nx;
    t_meshdata.m_nodes=t_meshdata.m_bulkelmts*t_meshdata.m_order+1;

    double dx,dy,dz;
    dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/(t_meshdata.m_nodes-1);
    dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/(t_meshdata.m_nodes-1);
    dz=(t_meshdata.m_zmax-t_meshdata.m_zmin)/(t_meshdata.m_nodes-1);

    t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_nx,vector<int>(t_meshdata.m_nodesperbulkelmt,0));
    t_meshdata.m_nodecoords.resize( t_meshdata.m_nodes*3,0.0);
    t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);

    for(int i=0;i<t_meshdata.m_nodes;i++){
        t_meshdata.m_nodecoords0[i*3+1-1]=t_meshdata.m_xmin+i*dx;
        t_meshdata.m_nodecoords0[i*3+2-1]=t_meshdata.m_ymin+i*dy;
        t_meshdata.m_nodecoords0[i*3+3-1]=t_meshdata.m_zmin+i*dz;

        t_meshdata.m_nodecoords[i*3+1-1]=t_meshdata.m_xmin+i*dx;
        t_meshdata.m_nodecoords[i*3+2-1]=t_meshdata.m_ymin+i*dy;
        t_meshdata.m_nodecoords[i*3+3-1]=t_meshdata.m_zmin+i*dz;
    }

    // generate the element connectivity information
    t_meshdata.m_pointelmts=2;
    t_meshdata.m_pointelmt_connectivity.resize(2);
    // for volume
    t_meshdata.m_pointelmt_volume.resize(2);
    t_meshdata.m_pointelmt_volume[0]=0.0;
    t_meshdata.m_pointelmt_volume[1]=0.0;
    // for left point
    leftconn.resize(1);
    leftconn[0].clear();
    leftconn[0].push_back(1);
    leftnodes.clear();
    leftnodes.push_back(1);
    t_meshdata.m_pointelmt_connectivity[0].clear();
    t_meshdata.m_pointelmt_connectivity[0].push_back(1);
    // for right point
    rightconn.resize(1);
    rightconn[0].clear();
    rightconn[0].push_back(t_meshdata.m_nodes);
    rightnodes.clear();
    rightnodes.push_back(t_meshdata.m_nodes);
    t_meshdata.m_pointelmt_connectivity[1].clear();
    t_meshdata.m_pointelmt_connectivity[1].push_back(t_meshdata.m_nodes);

    vector<int> tempconn;
    tempconn.clear();
    for(int e=0;e<t_meshdata.m_bulkelmts;e++){
        t_meshdata.m_bulkelmt_connectivity[e].clear();
        tempconn.push_back(e+1);
        for(int j=1;j<=t_meshdata.m_nodesperbulkelmt;j++){
            t_meshdata.m_bulkelmt_connectivity[e].push_back(e*t_meshdata.m_order+j);
        }
    }
    //***********************************************
    // set up the physical group information
    //***********************************************
    t_meshdata.m_phygroups=1+2;// two bc point+bulk elements
    t_meshdata.m_phygroup_dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyidvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phynamevec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_elmtnumvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_nodesnumperelmtvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2elmtconnvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyidvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phynamevec.resize(t_meshdata.m_phygroups);
    // for phy id and phy name vector
    t_meshdata.m_phygroup_phyidvec[0]=1;
    t_meshdata.m_phygroup_phyidvec[1]=2;
    t_meshdata.m_phygroup_phyidvec[2]=3;
    t_meshdata.m_phygroup_phynamevec[0]="left";
    t_meshdata.m_phygroup_phynamevec[1]="right";
    t_meshdata.m_phygroup_phynamevec[2]="alldomain";
    // now we can add physical group information
    t_meshdata.m_phygroup_dimvec[0]=0;// for left  point
    t_meshdata.m_phygroup_dimvec[1]=0;// for right point
    t_meshdata.m_phygroup_dimvec[2]=1;// for bulk element (1d-line)
    //*** element number for each phy group
    t_meshdata.m_phygroup_elmtnumvec[0]=1;
    t_meshdata.m_phygroup_elmtnumvec[1]=1;
    t_meshdata.m_phygroup_elmtnumvec[2]=t_meshdata.m_bulkelmts;
    //*** nodes number per elmt for each phy group
    t_meshdata.m_phygroup_nodesnumperelmtvec[0]=1;
    t_meshdata.m_phygroup_nodesnumperelmtvec[1]=1;
    t_meshdata.m_phygroup_nodesnumperelmtvec[2]=t_meshdata.m_nodesperbulkelmt;
    //*** for phy name to dim vec map
    t_meshdata.m_phygroup_name2dimvec[0]=make_pair("left",     0);
    t_meshdata.m_phygroup_name2dimvec[1]=make_pair("right",    0);
    t_meshdata.m_phygroup_name2dimvec[2]=make_pair("alldomain",1);
    //*** for phy name to element connectivity map
    t_meshdata.m_phygroup_name2elmtconnvec[0]=make_pair("left",                               leftconn);
    t_meshdata.m_phygroup_name2elmtconnvec[1]=make_pair("right",                             rightconn);
    t_meshdata.m_phygroup_name2elmtconnvec[2]=make_pair("alldomain",t_meshdata.m_bulkelmt_connectivity);
    //*** for phy name to bulk element id map
    t_meshdata.m_phygroup_name2bulkelmtidvec.clear();
    t_meshdata.m_phygroup_name2bulkelmtidvec.push_back(make_pair("alldomain",tempconn));
    tempconn.clear();
    //*** for phy name to id, and phyid to name map
    t_meshdata.m_phygroup_name2phyidvec.resize(1+2);
    t_meshdata.m_phygroup_name2phyidvec[0]=make_pair("left",     1);
    t_meshdata.m_phygroup_name2phyidvec[1]=make_pair("right",    2);
    t_meshdata.m_phygroup_name2phyidvec[2]=make_pair("alldomain",3);
    //***
    t_meshdata.m_phygroup_phyid2namevec.resize(1+2);
    t_meshdata.m_phygroup_phyid2namevec[0]=make_pair(1     ,"left");
    t_meshdata.m_phygroup_phyid2namevec[1]=make_pair(2    ,"right");
    t_meshdata.m_phygroup_phyid2namevec[2]=make_pair(3,"alldomain");

    //***********************************************
    // set up the node set physical group information
    //***********************************************
    t_meshdata.m_nodal_phygroups=2;
    t_meshdata.m_nodephygroup_name2nodeidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2nodeidvec[0]=make_pair("left",  leftnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[1]=make_pair("right",rightnodes);
    t_meshdata.m_nodephygroup_phynamevec.resize(2);
    t_meshdata.m_nodephygroup_phyidvec.resize(2);
    //*** for phy name to id map
    t_meshdata.m_nodephygroup_name2phyidvec.resize(2);
    t_meshdata.m_nodephygroup_name2phyidvec[0]=make_pair("left", 100001);
    t_meshdata.m_nodephygroup_name2phyidvec[1]=make_pair("right",100002);
    //*** for phyid to name map
    t_meshdata.m_nodephygroup_phyid2namevec.resize(2);
    t_meshdata.m_nodephygroup_phyid2namevec[0]=make_pair(100001 ,"left");
    t_meshdata.m_nodephygroup_phyid2namevec[1]=make_pair(100002,"right");

    t_meshdata.m_nodephygroup_phynamevec[0]="left";
    t_meshdata.m_nodephygroup_phynamevec[1]="right";
    t_meshdata.m_nodephygroup_phyidvec[0]=100001;
    t_meshdata.m_nodephygroup_phyidvec[1]=100002;

    m_mesh_generated=true;

    return m_mesh_generated;

}