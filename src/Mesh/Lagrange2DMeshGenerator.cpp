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
//+++ Date   : 2022.04.21
//+++ Purpose: the lagrange 2d mesh generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Lagrange2DMeshGenerator.h"

#include "Utils/MessagePrinter.h"

Lagrange2DMeshGenerator::Lagrange2DMeshGenerator(){
    m_mesh_generated=false;
    leftconn.clear();
    rightconn.clear();
    bottomconn.clear();
    topconn.clear();
    leftnodes.clear();
    rightnodes.clear();
    bottomnodes.clear();
    topnodes.clear();
}
Lagrange2DMeshGenerator::~Lagrange2DMeshGenerator(){
    leftconn.clear();
    rightconn.clear();
    bottomconn.clear();
    topconn.clear();
    leftnodes.clear();
    rightnodes.clear();
    bottomnodes.clear();
    topnodes.clear();
}

bool Lagrange2DMeshGenerator::generateMesh(const MeshType &t_meshtype,MeshData &t_meshdata){
    t_meshdata.m_maxdim=2;
    t_meshdata.m_mindim=1;

    t_meshdata.m_bulkelmt_connectivity.clear();
    t_meshdata.m_bulkelmt_volume.clear();

    t_meshdata.m_bulkelmt_type=t_meshtype;

    // init surface elements
    t_meshdata.m_surfaceelmt_connectivity.clear();
    t_meshdata.m_surfaceelmt_type=MeshType::NULLTYPE;
    t_meshdata.m_surfaceelmt_volume.clear();
    t_meshdata.m_surfaceelmt_vtktype=0;
    t_meshdata.m_surfaceelmts=0;
    // init line elements
    t_meshdata.m_lineelmt_connectivity.clear();
    t_meshdata.m_lineelmt_type=MeshType::NULLTYPE;
    t_meshdata.m_lineelmt_volume.clear();
    t_meshdata.m_lineelmt_vtktype=0;
    t_meshdata.m_lineelmts=0;
    // init point elements
    t_meshdata.m_pointelmts=0;
    t_meshdata.m_pointelmt_connectivity.clear();
    t_meshdata.m_pointelmt_volume.clear();
    // init nodes
    t_meshdata.m_nodesperlineelmt=0;
    t_meshdata.m_nodespersurfaceelmt=0;
    t_meshdata.m_nodesperbulkelmt=0;

    m_mesh_generated=false;

    double dx,dy;
    int i,j,k,e;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;

    vector<int> tempconn;

    if(t_meshtype==MeshType::QUAD4){
        t_meshdata.m_order=1;
        t_meshdata.m_bulkelmt_typename="quad4";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/t_meshdata.m_nx;
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/t_meshdata.m_ny;

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny;
        t_meshdata.m_lineelmts=2*(t_meshdata.m_nx+t_meshdata.m_ny);
        t_meshdata.m_surfaceelmts=0;
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_lineelmts;

        t_meshdata.m_nodes=(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
        t_meshdata.m_nodesperbulkelmt=4;
        t_meshdata.m_nodesperlineelmt=2;

        t_meshdata.m_bulkelmt_vtktype=9;
        t_meshdata.m_lineelmt_vtktype=3;
        t_meshdata.m_lineelmt_type=MeshType::EDGE2;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize(t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        for(j=1;j<=t_meshdata.m_ny+1;j++){
            for(i=1;i<=t_meshdata.m_nx+1;i++){
                k=(j-1)*(t_meshdata.m_nx+1)+i;
                t_meshdata.m_nodecoords[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                t_meshdata.m_nodecoords[(k-1)*3+3-1]=0.0;

                t_meshdata.m_nodecoords0[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords0[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                t_meshdata.m_nodecoords0[(k-1)*3+3-1]=0.0;

                if(i==1){
                    // for left side nodes
                    leftnodes.push_back(k);// global id, start from 1
                }
                if(i==t_meshdata.m_nx+1){
                    // for right side nodes
                    rightnodes.push_back(k);
                }
                if(j==1){
                    // for bottom side nodes
                    bottomnodes.push_back(k);
                }
                if(j==t_meshdata.m_ny+1){
                    // for top side nodes
                    topnodes.push_back(k);
                }
            }
        }
        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny);rightconn.resize(t_meshdata.m_ny);
        bottomconn.resize(t_meshdata.m_nx);topconn.resize(t_meshdata.m_nx);
        tempconn.clear();
        for(j=1;j<=t_meshdata.m_ny;j++){
            for(i=1;i<=t_meshdata.m_nx;i++){
                e=(j-1)*t_meshdata.m_nx+i;
                i1=(j-1)*(t_meshdata.m_nx+1)+i;
                i2=i1+1;
                i3=i2+t_meshdata.m_nx+1;
                i4=i3-1;

                tempconn.push_back(e);

                t_meshdata.m_bulkelmt_connectivity[e-1].clear();
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i1);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i2);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i3);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i4);

                // for the boundary element
                // the layout of your quad4 should be:
                // 4-----3
                // |     |
                // |     |
                // 1-----2
                if(j==1){
                    // for bottom bc elements
                    bottomconn[i-1].clear();
                    bottomconn[i-1].push_back(i1);
                    bottomconn[i-1].push_back(i2);
                }
                if(j==t_meshdata.m_ny){
                    // for top bc elements
                    topconn[i-1].clear();
                    topconn[i-1].push_back(i3);
                    topconn[i-1].push_back(i4);
                }
                if(i==1){
                    // for left bc elements
                    leftconn[j-1].clear();
                    leftconn[j-1].push_back(i4);
                    leftconn[j-1].push_back(i1);
                }
                if(i==t_meshdata.m_nx){
                    // for right bc elements
                    rightconn[j-1].clear();
                    rightconn[j-1].push_back(i2);
                    rightconn[j-1].push_back(i3);
                }
            }
        }
    }
    else if(t_meshtype==MeshType::QUAD8){
        t_meshdata.m_order=2;
        t_meshdata.m_bulkelmt_typename="quad8";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/(2.0*t_meshdata.m_nx);
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/(2.0*t_meshdata.m_ny);

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny;
        t_meshdata.m_lineelmts=2*(t_meshdata.m_nx+t_meshdata.m_ny);
        t_meshdata.m_surfaceelmts=0;
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_lineelmts;

        t_meshdata.m_nodes=(2*t_meshdata.m_nx+1)*(2*t_meshdata.m_ny+1)-t_meshdata.m_bulkelmts;
        t_meshdata.m_nodesperbulkelmt=8;
        t_meshdata.m_nodesperlineelmt=3;

        t_meshdata.m_bulkelmt_vtktype=23;
        t_meshdata.m_lineelmt_vtktype=4;
        t_meshdata.m_lineelmt_type=MeshType::EDGE3;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize( t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        for(j=1;j<=t_meshdata.m_ny;j++){
            // for bottom line of each element
            for(i=1;i<=2*t_meshdata.m_nx+1;i++){
                k=(j-1)*(2*t_meshdata.m_nx+1+t_meshdata.m_nx+1)+i;
                t_meshdata.m_nodecoords[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
                t_meshdata.m_nodecoords[(k-1)*3+3-1]=0.0;
            }
            // for middle line of each element
            for(i=1;i<=t_meshdata.m_nx+1;i++){
                k=(j-1)*(2*t_meshdata.m_nx+1+t_meshdata.m_nx+1)+2*t_meshdata.m_nx+1+i;
                t_meshdata.m_nodecoords[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*2*dx;
                t_meshdata.m_nodecoords[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy+dy;
                t_meshdata.m_nodecoords[(k-1)*3+3-1]=0.0;
            }
        }
        // for the last top line
        j=t_meshdata.m_ny+1;
        for(i=1;i<=2*t_meshdata.m_nx+1;i++){
            k=(j-1)*(2*t_meshdata.m_nx+1+t_meshdata.m_nx+1)+i;
            t_meshdata.m_nodecoords[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
            t_meshdata.m_nodecoords[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
            t_meshdata.m_nodecoords[(k-1)*3+3-1]=0.0;
        }
        // make a copy for another nodecoords
        t_meshdata.m_nodecoords0=t_meshdata.m_nodecoords;

        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny);rightconn.resize(t_meshdata.m_ny);
        bottomconn.resize(t_meshdata.m_nx);topconn.resize(t_meshdata.m_nx);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();topnodes.clear();
        tempconn.clear();
        for(j=1;j<=t_meshdata.m_ny;j++){
            for(i=1;i<=t_meshdata.m_nx;i++){
                e=(j-1)*t_meshdata.m_nx+i;
                i1=(j-1)*(2*t_meshdata.m_nx+1+t_meshdata.m_nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+(2*t_meshdata.m_nx+1+t_meshdata.m_nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*t_meshdata.m_nx+1)-i;
                i7=i3-1;
                i8=i1+(2*t_meshdata.m_nx+1)-(i-1);

                tempconn.push_back(e);

                t_meshdata.m_bulkelmt_connectivity[e-1].clear();
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i1);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i2);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i3);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i4);

                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i5);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i6);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i7);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i8);

                // for the boundary element and nodes
                // the layout of your quad8 should be:
                // 4--7--3
                // |     |
                // 8     6
                // |     |
                // 1--5--2
                if(j==1){
                    // for bottom bc elements
                    bottomconn[i-1].clear();
                    bottomconn[i-1].push_back(i1);
                    bottomconn[i-1].push_back(i5);
                    bottomconn[i-1].push_back(i2);

                    bottomnodes.push_back(i1);
                    bottomnodes.push_back(i5);
                    bottomnodes.push_back(i2);
                }
                if(j==t_meshdata.m_ny){
                    // for top bc elements
                    topconn[i-1].clear();
                    topconn[i-1].push_back(i3);
                    topconn[i-1].push_back(i7);
                    topconn[i-1].push_back(i4);

                    topnodes.push_back(i3);
                    topnodes.push_back(i7);
                    topnodes.push_back(i4);
                }
                if(i==1){
                    // for left bc elements
                    leftconn[j-1].clear();
                    leftconn[j-1].push_back(i4);
                    leftconn[j-1].push_back(i8);
                    leftconn[j-1].push_back(i1);

                    leftnodes.push_back(i4);
                    leftnodes.push_back(i8);
                    leftnodes.push_back(i1);
                }
                if(i==t_meshdata.m_nx){
                    // for right bc elements
                    rightconn[j-1].clear();
                    rightconn[j-1].push_back(i2);
                    rightconn[j-1].push_back(i6);
                    rightconn[j-1].push_back(i3);

                    rightnodes.push_back(i2);
                    rightnodes.push_back(i6);
                    rightnodes.push_back(i3);
                }
            }
        }
    }
    else if(t_meshtype==MeshType::QUAD9){
        t_meshdata.m_order=2;
        t_meshdata.m_bulkelmt_typename="quad9";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/(2.0*t_meshdata.m_nx);
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/(2.0*t_meshdata.m_ny);

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny;
        t_meshdata.m_lineelmts=2*(t_meshdata.m_nx+t_meshdata.m_ny);
        t_meshdata.m_surfaceelmts=0;
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_lineelmts;

        t_meshdata.m_nodes=(2*t_meshdata.m_nx+1)*(2*t_meshdata.m_ny+1);
        t_meshdata.m_nodesperbulkelmt=9;
        t_meshdata.m_nodesperlineelmt=3;

        t_meshdata.m_bulkelmt_vtktype=28;
        t_meshdata.m_lineelmt_vtktype=4;
        t_meshdata.m_lineelmt_type=MeshType::EDGE3;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize( t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        for(j=1;j<=2*t_meshdata.m_ny+1;j++){
            for(i=1;i<=2*t_meshdata.m_nx+1;i++){
                k=(j-1)*(2*t_meshdata.m_nx+1)+i;
                t_meshdata.m_nodecoords[(k-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords[(k-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                t_meshdata.m_nodecoords[(k-1)*3+3-1]=0.0;
                if(i==1){
                    // for left nodes
                    leftnodes.push_back(k);
                }
                if(i==2*t_meshdata.m_nx+1){
                    // for right nodes
                    rightnodes.push_back(k);
                }
                if(j==1){
                    // for bottom nodes
                    bottomnodes.push_back(k);
                }
                if(j==2*t_meshdata.m_ny+1){
                    // for top nodes
                    topnodes.push_back(k);
                }
            }
        }
        // make a copy for coordinate0
        t_meshdata.m_nodecoords0=t_meshdata.m_nodecoords;

        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny);rightconn.resize(t_meshdata.m_ny);
        bottomconn.resize(t_meshdata.m_nx);topconn.resize(t_meshdata.m_nx);
        tempconn.clear();
        for(j=1;j<=t_meshdata.m_ny;j++){
            for(i=1;i<=t_meshdata.m_nx;i++){
                e=(j-1)*t_meshdata.m_nx+i;
                i1=(j-1)*2*(2*t_meshdata.m_nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+2*(2*t_meshdata.m_nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*t_meshdata.m_nx+1);
                i7=i3-1;
                i8=i1+(2*t_meshdata.m_nx+1);
                i9=i8+1;

                tempconn.push_back(e);

                t_meshdata.m_bulkelmt_connectivity[e-1].clear();
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i1);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i2);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i3);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i4);

                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i5);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i6);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i7);
                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i8);

                t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i9);

                // for the boundary element
                // the layout of your quad8 should be:
                // 4--7--3
                // |     |
                // 8  9  6
                // |     |
                // 1--5--2
                if(j==1){
                    // for bottom bc elements
                    bottomconn[i-1].clear();
                    bottomconn[i-1].push_back(i1);
                    bottomconn[i-1].push_back(i5);
                    bottomconn[i-1].push_back(i2);
                }
                if(j==t_meshdata.m_ny){
                    // for top bc elements
                    topconn[i-1].clear();
                    topconn[i-1].push_back(i3);
                    topconn[i-1].push_back(i7);
                    topconn[i-1].push_back(i4);
                }
                if(i==1){
                    // for left bc elements
                    leftconn[j-1].clear();
                    leftconn[j-1].push_back(i4);
                    leftconn[j-1].push_back(i8);
                    leftconn[j-1].push_back(i1);
                }
                if(i==t_meshdata.m_nx){
                    // for right bc elements
                    rightconn[j-1].clear();
                    rightconn[j-1].push_back(i2);
                    rightconn[j-1].push_back(i6);
                    rightconn[j-1].push_back(i3);
                }
            }
        }
    }
    else{
        MessagePrinter::printErrorTxt("unsupported mesh generation in 2d case. Currently, only quad4,quad8, and quad9 are supported.");
        m_mesh_generated=false;
        return m_mesh_generated;
    }

    // remove the duplicate node id, for nodal type set, you don't need these duplicate node id!
    // for leftnodes
    sort(leftnodes.begin(),leftnodes.end());
    leftnodes.erase(unique(leftnodes.begin(),leftnodes.end()),leftnodes.end());
    // for rightnodes
    sort(rightnodes.begin(),rightnodes.end());
    rightnodes.erase(unique(rightnodes.begin(),rightnodes.end()),rightnodes.end());
    // for bottomnodes
    sort(bottomnodes.begin(),bottomnodes.end());
    bottomnodes.erase(unique(bottomnodes.begin(),bottomnodes.end()),bottomnodes.end());
    // for topnodes
    sort(topnodes.begin(),topnodes.end());
    topnodes.erase(unique(topnodes.begin(),topnodes.end()),topnodes.end());

    
    //***********************************************
    // set up the physical group information
    //***********************************************
    t_meshdata.m_phygroups=1+4;// two bc point+bulk elements
    t_meshdata.m_phygroup_dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyidvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phynamevec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_elmtnumvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_nodesnumperelmtvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2elmtconnvec.resize(t_meshdata.m_phygroups);
    // now we can add physical group information
    t_meshdata.m_phygroup_dimvec[0]=1;// for left   edge
    t_meshdata.m_phygroup_dimvec[1]=1;// for right  edge
    t_meshdata.m_phygroup_dimvec[2]=1;// for bottom edge
    t_meshdata.m_phygroup_dimvec[3]=1;// for top    edge
    t_meshdata.m_phygroup_dimvec[4]=2;// for bulk element (2d-surface)
    //*** for phy name vector
    t_meshdata.m_phygroup_phynamevec[0]="left";
    t_meshdata.m_phygroup_phynamevec[1]="right";
    t_meshdata.m_phygroup_phynamevec[2]="bottom";
    t_meshdata.m_phygroup_phynamevec[3]="top";
    t_meshdata.m_phygroup_phynamevec[4]="alldomain";
    //*** for phy id vector
    t_meshdata.m_phygroup_phyidvec[0]=1;
    t_meshdata.m_phygroup_phyidvec[1]=2;
    t_meshdata.m_phygroup_phyidvec[2]=3;
    t_meshdata.m_phygroup_phyidvec[3]=4;
    t_meshdata.m_phygroup_phyidvec[4]=5;
    //*** element number for each phy group
    t_meshdata.m_phygroup_elmtnumvec[0]=t_meshdata.m_ny;// for left   edge
    t_meshdata.m_phygroup_elmtnumvec[1]=t_meshdata.m_ny;// for right  edge
    t_meshdata.m_phygroup_elmtnumvec[2]=t_meshdata.m_nx;// for bottom edge
    t_meshdata.m_phygroup_elmtnumvec[3]=t_meshdata.m_nx;// for top    edge
    t_meshdata.m_phygroup_elmtnumvec[4]=t_meshdata.m_bulkelmts;
    //*** nodes number per elmt for each phy group
    t_meshdata.m_phygroup_nodesnumperelmtvec[0]=t_meshdata.m_nodesperlineelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[1]=t_meshdata.m_nodesperlineelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[2]=t_meshdata.m_nodesperlineelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[3]=t_meshdata.m_nodesperlineelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[4]=t_meshdata.m_nodesperbulkelmt;
    //*** for phy name to dim vec map
    t_meshdata.m_phygroup_name2dimvec[0]=make_pair("left",     1);
    t_meshdata.m_phygroup_name2dimvec[1]=make_pair("right",    1);
    t_meshdata.m_phygroup_name2dimvec[2]=make_pair("bottom",   1);
    t_meshdata.m_phygroup_name2dimvec[3]=make_pair("top",      1);
    t_meshdata.m_phygroup_name2dimvec[4]=make_pair("alldomain",2);
    //*** for phy name to element connectivity map
    t_meshdata.m_phygroup_name2elmtconnvec[0]=make_pair("left",                               leftconn);
    t_meshdata.m_phygroup_name2elmtconnvec[1]=make_pair("right",                             rightconn);
    t_meshdata.m_phygroup_name2elmtconnvec[2]=make_pair("bottom",                           bottomconn);
    t_meshdata.m_phygroup_name2elmtconnvec[3]=make_pair("top",                                 topconn);
    t_meshdata.m_phygroup_name2elmtconnvec[4]=make_pair("alldomain",t_meshdata.m_bulkelmt_connectivity);
    //*** for phy name to bulk element id map
    t_meshdata.m_phygroup_name2bulkelmtidvec.clear();
    t_meshdata.m_phygroup_name2bulkelmtidvec.push_back(make_pair("alldomain",tempconn));
    tempconn.clear();
    //*** for phy name to id, and phyid to name map
    t_meshdata.m_phygroup_name2phyidvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2phyidvec[0]=make_pair("left",     1);
    t_meshdata.m_phygroup_name2phyidvec[1]=make_pair("right",    2);
    t_meshdata.m_phygroup_name2phyidvec[2]=make_pair("bottom",   3);
    t_meshdata.m_phygroup_name2phyidvec[3]=make_pair("top",      4);
    t_meshdata.m_phygroup_name2phyidvec[4]=make_pair("alldomain",5);
    //***
    t_meshdata.m_phygroup_phyid2namevec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyid2namevec[0]=make_pair(1     ,"left");
    t_meshdata.m_phygroup_phyid2namevec[1]=make_pair(2    ,"right");
    t_meshdata.m_phygroup_phyid2namevec[2]=make_pair(3   ,"bottom");
    t_meshdata.m_phygroup_phyid2namevec[3]=make_pair(4      ,"top");
    t_meshdata.m_phygroup_phyid2namevec[4]=make_pair(5,"alldomain");

    //***********************************************
    // set up the node set physical group information
    //***********************************************
    t_meshdata.m_nodal_phygroups=4;
    t_meshdata.m_nodephygroup_phynamevec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_phyidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2nodeidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2nodeidvec[0]=make_pair("left",    leftnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[1]=make_pair("right",  rightnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[2]=make_pair("bottom",bottomnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[3]=make_pair("top",      topnodes);
    
    //*** phy id vector
    t_meshdata.m_nodephygroup_phyidvec[0]=10001;
    t_meshdata.m_nodephygroup_phyidvec[1]=20001;
    t_meshdata.m_nodephygroup_phyidvec[2]=30001;
    t_meshdata.m_nodephygroup_phyidvec[3]=40001;
    //*** phy name vector
    t_meshdata.m_nodephygroup_phynamevec[0]="left";
    t_meshdata.m_nodephygroup_phynamevec[1]="right";
    t_meshdata.m_nodephygroup_phynamevec[2]="bottom";
    t_meshdata.m_nodephygroup_phynamevec[3]="top";
    //*** for phy name to id map
    t_meshdata.m_nodephygroup_name2phyidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2phyidvec[0]=make_pair("left",  10001);
    t_meshdata.m_nodephygroup_name2phyidvec[1]=make_pair("right", 20001);
    t_meshdata.m_nodephygroup_name2phyidvec[2]=make_pair("bottom",30001);
    t_meshdata.m_nodephygroup_name2phyidvec[3]=make_pair("top",   40001);
    //*** for phyid to name map
    t_meshdata.m_nodephygroup_phyid2namevec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_phyid2namevec[0]=make_pair(10001  ,"left");
    t_meshdata.m_nodephygroup_phyid2namevec[1]=make_pair(20001 ,"right");
    t_meshdata.m_nodephygroup_phyid2namevec[2]=make_pair(30001,"bottom");
    t_meshdata.m_nodephygroup_phyid2namevec[3]=make_pair(40001   ,"top");

    m_mesh_generated=true;

    return m_mesh_generated;

}