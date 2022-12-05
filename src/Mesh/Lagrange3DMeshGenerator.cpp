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
//+++ Date   : 2022.05.04
//+++ Purpose: the lagrange 3d mesh generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Lagrange3DMeshGenerator.h"

#include "Utils/MessagePrinter.h"

Lagrange3DMeshGenerator::Lagrange3DMeshGenerator(){
    m_mesh_generated=false;
    leftconn.clear();
    rightconn.clear();
    bottomconn.clear();
    topconn.clear();
    backconn.clear();
    frontconn.clear();
    //
    leftnodes.clear();
    rightnodes.clear();
    bottomnodes.clear();
    topnodes.clear();
    backnodes.clear();
    frontnodes.clear();
}
Lagrange3DMeshGenerator::~Lagrange3DMeshGenerator(){
    leftconn.clear();
    rightconn.clear();
    bottomconn.clear();
    topconn.clear();
    backconn.clear();
    frontconn.clear();
    //
    leftnodes.clear();
    rightnodes.clear();
    bottomnodes.clear();
    topnodes.clear();
    backnodes.clear();
    frontnodes.clear();
}

bool Lagrange3DMeshGenerator::generateMesh(const MeshType &t_meshtype,MeshData &t_meshdata){
    t_meshdata.m_maxdim=3;
    t_meshdata.m_mindim=2;

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

    double dx,dy,dz;
    int i,j,k,kk,e;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    int i10,i11,i12,i13,i14,i15,i16,i17,i18,i19;
    int i20,i21,i22,i23,i24,i25,i26,i27;

    vector<int> tempconn;

    if(t_meshtype==MeshType::HEX8){
        t_meshdata.m_order=1;
        t_meshdata.m_bulkelmt_typename="hex8";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/t_meshdata.m_nx;
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/t_meshdata.m_ny;
        dz=(t_meshdata.m_zmax-t_meshdata.m_zmin)/t_meshdata.m_nz;

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny*t_meshdata.m_nz;
        t_meshdata.m_lineelmts=0;
        t_meshdata.m_surfaceelmts=2*(t_meshdata.m_nx*t_meshdata.m_ny
                                    +t_meshdata.m_nx*t_meshdata.m_nz
                                    +t_meshdata.m_ny*t_meshdata.m_nz);
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_surfaceelmts;

        t_meshdata.m_nodes=(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1)*(t_meshdata.m_nz+1);
        t_meshdata.m_nodesperbulkelmt=8;
        t_meshdata.m_nodespersurfaceelmt=4;

        t_meshdata.m_bulkelmt_vtktype=12;
        
        t_meshdata.m_lineelmt_vtktype=3;
        t_meshdata.m_lineelmt_type=MeshType::EDGE2;
        t_meshdata.m_surfaceelmt_vtktype=9;
        t_meshdata.m_surfaceelmt_type=MeshType::QUAD4;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize(t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
        for(k=1;k<=t_meshdata.m_nz+1;k++){
            for(j=1;j<=t_meshdata.m_ny+1;j++){
                for(i=1;i<=t_meshdata.m_nx+1;i++){
                    kk=(j-1)*(t_meshdata.m_nx+1)+i+(k-1)*(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*dz;

                    if(i==1){
                        // for left side nodes
                        leftnodes.push_back(kk);// global id, start from 1
                    }
                    if(i==t_meshdata.m_nx+1){
                        // for right side nodes
                        rightnodes.push_back(kk);
                    }
                    if(j==1){
                        // for bottom side nodes
                        bottomnodes.push_back(kk);
                    }
                    if(j==t_meshdata.m_ny+1){
                        // for top side nodes
                        topnodes.push_back(kk);
                    }
                    if(k==1){
                        // for back side nodes
                        backnodes.push_back(kk);
                    }
                    if(k==t_meshdata.m_nz+1){
                        // for front side nodes
                        frontnodes.push_back(kk);
                    }
                }
            }
        }
        // make a copy for nodcal coordinates
        t_meshdata.m_nodecoords=t_meshdata.m_nodecoords0;
        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        rightconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        //
        bottomconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        topconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        //
        backconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);
        frontconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);

        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();

        tempconn.clear();

        for(k=1;k<=t_meshdata.m_nz;k++){
            for(j=1;j<=t_meshdata.m_ny;j++){
                for(i=1;i<=t_meshdata.m_nx;i++){
                    e=(j-1)*t_meshdata.m_nx+i+(k-1)*t_meshdata.m_nx*t_meshdata.m_ny;
                    i1=(j-1)*(t_meshdata.m_nx+1)+i+(k-1)*(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
                    i2=i1+1;
                    i3=i2+t_meshdata.m_nx+1;
                    i4=i3-1;
                    i5=i1+(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
                    i6=i2+(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
                    i7=i3+(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);
                    i8=i4+(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);

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

                    t_meshdata.m_bulkelmt_volume[e-1]=dx*dy*dz;// for the volume of e-th bulk element

                    if(i==1){
                        // for left bc elements
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i1);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i5);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i8);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i4);

                        leftnodes.push_back(i1);
                        leftnodes.push_back(i5);
                        leftnodes.push_back(i8);
                        leftnodes.push_back(i4);
                    }
                    if(i==t_meshdata.m_nx){
                        // for right bc elements
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i2);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i3);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i7);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i6);

                        rightnodes.push_back(i2);
                        rightnodes.push_back(i3);
                        rightnodes.push_back(i7);
                        rightnodes.push_back(i6);
                    }
                    if(j==1){
                        // for bottom bc elements
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i2);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i5);

                        bottomnodes.push_back(i1);
                        bottomnodes.push_back(i2);
                        bottomnodes.push_back(i6);
                        bottomnodes.push_back(i5);
                    }
                    if(j==t_meshdata.m_ny){
                        // for bottom bc elements
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i8);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i3);

                        topnodes.push_back(i4);
                        topnodes.push_back(i8);
                        topnodes.push_back(i7);
                        topnodes.push_back(i3);
                    }
                    if(k==1){
                        // for back bc elements
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i3);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i2);

                        backnodes.push_back(i1);
                        backnodes.push_back(i4);
                        backnodes.push_back(i3);
                        backnodes.push_back(i2);
                    }
                    if(k==t_meshdata.m_nz){
                        // for front bc elements
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i5);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i8);

                        frontnodes.push_back(i5);
                        frontnodes.push_back(i6);
                        frontnodes.push_back(i7);
                        frontnodes.push_back(i8);
                    }
                }
            }
        }
    }
    else if(t_meshtype==MeshType::HEX20){
        t_meshdata.m_order=2;
        t_meshdata.m_bulkelmt_typename="hex20";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/(2.0*t_meshdata.m_nx);
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/(2.0*t_meshdata.m_ny);
        dz=(t_meshdata.m_zmax-t_meshdata.m_zmin)/(2.0*t_meshdata.m_nz);

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny*t_meshdata.m_nz;
        t_meshdata.m_lineelmts=0;
        t_meshdata.m_surfaceelmts=2*(t_meshdata.m_nx*t_meshdata.m_ny
                                    +t_meshdata.m_nx*t_meshdata.m_nz
                                    +t_meshdata.m_ny*t_meshdata.m_nz);
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_surfaceelmts;

        int nLayer1Nodes=(2*t_meshdata.m_nx+1)*(2*t_meshdata.m_ny+1)-t_meshdata.m_nx*t_meshdata.m_ny;// for norm layer
        int nLayer2Nodes=(t_meshdata.m_nx+1)*(t_meshdata.m_ny+1);          // for middle layer

        t_meshdata.m_nodes=nLayer1Nodes*(t_meshdata.m_nz+1)+nLayer2Nodes*t_meshdata.m_nz;
        t_meshdata.m_nodesperbulkelmt=20;
        t_meshdata.m_nodespersurfaceelmt=8;

        t_meshdata.m_bulkelmt_vtktype=25;

        t_meshdata.m_lineelmt_vtktype=4;
        t_meshdata.m_lineelmt_type=MeshType::EDGE3;
        t_meshdata.m_surfaceelmt_vtktype=23;
        t_meshdata.m_surfaceelmt_type=MeshType::QUAD9;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize(t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
        int _Nz,_Ny,_Nx;
        _Nz=t_meshdata.m_nz;_Ny=t_meshdata.m_ny;_Nx=t_meshdata.m_nx;
        for(k=1;k<=_Nz;++k){
            // First for normal layer
            for(j=1;j<=_Ny;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;i++){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
                }
                // for middle line of each element
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*2*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy+dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
                }
            }
            // for top line
            j=_Ny+1;
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
                t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
            }
            // Then for middle type layer
            for(j=1;j<=_Ny+1;++j){
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*2*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=_Nz+1;
        for(j=1;j<=_Ny;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
                t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
            }
            // for middle line of each element
            for(i=1;i<=_Nx+1;++i){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*2*dx;
                t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy+dy;
                t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
            }
        }
        // for top line
        j=_Ny+1;
        for(i=1;i<=2*_Nx+1;i++){
            kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
            t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
            t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*2*dy;
            t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
        }
        // make a copy for coordniate
        t_meshdata.m_nodecoords=t_meshdata.m_nodecoords0;
        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        rightconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        //
        bottomconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        topconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        //
        backconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);
        frontconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);

        tempconn.clear();

        for(k=1;k<=_Nz;++k){
            for(j=1;j<=_Ny;++j){
                for(i=1;i<=_Nx;++i){
                    e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
                    i1=(j-1)*(2*_Nx+1+_Nx+1)+2*i-1+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    i2=i1+2;
                    i3=i2+(2*_Nx+1+_Nx+1);
                    i4=i3-2;

                    i5=i1+nLayer1Nodes+nLayer2Nodes;
                    i6=i2+nLayer1Nodes+nLayer2Nodes;
                    i7=i3+nLayer1Nodes+nLayer2Nodes;
                    i8=i4+nLayer1Nodes+nLayer2Nodes;

                    i9 =i1+1;
                    i10=i2+(2*_Nx+1-i);
                    i11=i3-1;
                    i12=i10-1;

                    i13= i9+nLayer1Nodes+nLayer2Nodes;
                    i14=i10+nLayer1Nodes+nLayer2Nodes;
                    i15=i11+nLayer1Nodes+nLayer2Nodes;
                    i16=i12+nLayer1Nodes+nLayer2Nodes;

                    i17=i1+nLayer1Nodes-(i-1+(j-1)*(_Nx+_Nx+1));
                    i18=i17+1;
                    i19=i18+_Nx+1;
                    i20=i19-1;

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
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i10);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i11);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i12);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i13);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i14);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i15);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i16);

                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i17);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i18);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i19);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i20);
                   
                    t_meshdata.m_bulkelmt_volume[e-1]=2.0*dx*2.0*dy*2.0*dz;

                    if(i==1){
                        // for left bc elements
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i1);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i5);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i8);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i4);

                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i17);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i16);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i20);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i12);

                        leftnodes.push_back(i1);
                        leftnodes.push_back(i5);
                        leftnodes.push_back(i8);
                        leftnodes.push_back(i4);
                        leftnodes.push_back(i17);
                        leftnodes.push_back(i16);
                        leftnodes.push_back(i20);
                        leftnodes.push_back(i12);
                    }
                    if(i==t_meshdata.m_nx){
                        // for right bc elements
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i2);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i3);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i7);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i6);

                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i10);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i19);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i14);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i18);

                        rightnodes.push_back(i2);
                        rightnodes.push_back(i3);
                        rightnodes.push_back(i7);
                        rightnodes.push_back(i6);
                        rightnodes.push_back(i10);
                        rightnodes.push_back(i19);
                        rightnodes.push_back(i14);
                        rightnodes.push_back(i18);
                    }
                    if(j==1){
                        // for bottom bc elements
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i2);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i5);

                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i9 );
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i18);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i13);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i17);

                        bottomnodes.push_back(i1);
                        bottomnodes.push_back(i2);
                        bottomnodes.push_back(i6);
                        bottomnodes.push_back(i5);
                        bottomnodes.push_back(i9 );
                        bottomnodes.push_back(i18);
                        bottomnodes.push_back(i13);
                        bottomnodes.push_back(i17);
                    }
                    if(j==t_meshdata.m_ny){
                        // for top bc elements
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i8);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i3);

                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i20);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i15);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i19);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i11);

                        topnodes.push_back(i4);
                        topnodes.push_back(i8);
                        topnodes.push_back(i7);
                        topnodes.push_back(i3);
                        topnodes.push_back(i20);
                        topnodes.push_back(i15);
                        topnodes.push_back(i19);
                        topnodes.push_back(i11);
                    }
                    if(k==1){
                        // for back bc elements
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i3);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i2);

                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i12);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i11);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i10);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i9 );

                        backnodes.push_back(i1);
                        backnodes.push_back(i4);
                        backnodes.push_back(i3);
                        backnodes.push_back(i2);
                        backnodes.push_back(i12);
                        backnodes.push_back(i11);
                        backnodes.push_back(i10);
                        backnodes.push_back(i9 );
                    }
                    if(k==t_meshdata.m_nz){
                        // for front bc elements
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i5);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i8);

                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i13);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i14);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i15);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i16);

                        frontnodes.push_back(i5);
                        frontnodes.push_back(i6);
                        frontnodes.push_back(i7);
                        frontnodes.push_back(i8);
                        frontnodes.push_back(i13);
                        frontnodes.push_back(i14);
                        frontnodes.push_back(i15);
                        frontnodes.push_back(i16);
                    }

                }
            }
        }
    }
    else if(t_meshtype==MeshType::HEX27){
        t_meshdata.m_order=2;
        t_meshdata.m_bulkelmt_typename="hex27";

        dx=(t_meshdata.m_xmax-t_meshdata.m_xmin)/(2.0*t_meshdata.m_nx);
        dy=(t_meshdata.m_ymax-t_meshdata.m_ymin)/(2.0*t_meshdata.m_ny);
        dz=(t_meshdata.m_zmax-t_meshdata.m_zmin)/(2.0*t_meshdata.m_nz);

        t_meshdata.m_bulkelmts=t_meshdata.m_nx*t_meshdata.m_ny*t_meshdata.m_nz;
        t_meshdata.m_lineelmts=0;
        t_meshdata.m_surfaceelmts=2*(t_meshdata.m_nx*t_meshdata.m_ny
                                    +t_meshdata.m_nx*t_meshdata.m_nz
                                    +t_meshdata.m_ny*t_meshdata.m_nz);
        t_meshdata.m_elements=t_meshdata.m_bulkelmts+t_meshdata.m_surfaceelmts;

        int nLayerNodes=(2*t_meshdata.m_nx+1)*(2*t_meshdata.m_ny+1);// for norm layer

        t_meshdata.m_nodes=(2*t_meshdata.m_nz+1)*nLayerNodes;
        t_meshdata.m_nodesperbulkelmt=27;
        t_meshdata.m_nodespersurfaceelmt=9;

        t_meshdata.m_bulkelmt_vtktype=29;

        t_meshdata.m_lineelmt_vtktype=4;
        t_meshdata.m_lineelmt_type=MeshType::EDGE3;
        t_meshdata.m_surfaceelmt_vtktype=23;
        t_meshdata.m_surfaceelmt_type=MeshType::QUAD9;

        // for the coordinates of each node
        t_meshdata.m_nodecoords.resize(t_meshdata.m_nodes*3,0.0);
        t_meshdata.m_nodecoords0.resize(t_meshdata.m_nodes*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
        int _Nz,_Ny,_Nx;
        _Nz=t_meshdata.m_nz;_Ny=t_meshdata.m_ny;_Nx=t_meshdata.m_nx;
        for(k=1;k<=_Nz;++k){
            // For first layer
            for(j=1;j<=2*_Ny+1;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
                }
            }
            // Then for second layer
            for(j=1;j<=2*_Ny+1;++j){
                for(i=1;i<=2*_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes+nLayerNodes;
                    t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                    t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                    t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=_Nz+1;
        for(j=1;j<=2*_Ny+1;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;++i){
                kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                t_meshdata.m_nodecoords0[(kk-1)*3+1-1]=t_meshdata.m_xmin+(i-1)*dx;
                t_meshdata.m_nodecoords0[(kk-1)*3+2-1]=t_meshdata.m_ymin+(j-1)*dy;
                t_meshdata.m_nodecoords0[(kk-1)*3+3-1]=t_meshdata.m_zmin+(k-1)*2*dz;
            }
        }
        // make a copy for coordniate
        t_meshdata.m_nodecoords=t_meshdata.m_nodecoords0;
        // for the connectivity information of bulk elements
        t_meshdata.m_bulkelmt_connectivity.resize(t_meshdata.m_bulkelmts,vector<int>(0));
        t_meshdata.m_bulkelmt_volume.resize(t_meshdata.m_bulkelmts,0.0);
        leftconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        rightconn.resize(t_meshdata.m_ny*t_meshdata.m_nz);
        //
        bottomconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        topconn.resize(t_meshdata.m_nx*t_meshdata.m_nz);
        //
        backconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);
        frontconn.resize(t_meshdata.m_nx*t_meshdata.m_ny);

        tempconn.clear();

        for(k=1;k<=_Nz;++k){
            for(j=1;j<=_Ny;++j){
                for(i=1;i<=_Nx;++i){
                    e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
                    i1=(j-1)*2*(2*_Nx+1)+2*i-1+(k-1)*2*nLayerNodes;
                    i2=i1+2;
                    i3=i2+(2*_Nx+1)*2;
                    i4=i3-2;

                    i5=i1+2*nLayerNodes;
                    i6=i2+2*nLayerNodes;
                    i7=i3+2*nLayerNodes;
                    i8=i4+2*nLayerNodes;

                    i9 =i1+1;
                    i10=i2+(2*_Nx+1);
                    i11=i3-1;
                    i12=i1+(2*_Nx+1);

                    i13=i5+1;
                    i14=i6+(2*_Nx+1);
                    i15=i7-1;
                    i16=i5+(2*_Nx+1);

                    i17=i1+nLayerNodes;
                    i18=i2+nLayerNodes;
                    i19=i3+nLayerNodes;
                    i20=i4+nLayerNodes;

                    i21=i17+(2*_Nx+1);
                    i22=i21+2;

                    //i23=i20+1;
                    //i24=i17+1;

                    i23=i17+1;
                    i24=i20+1;

                    i25=i12+1;
                    i26=i16+1;

                    i27=i21+1;

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
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i10);

                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i11);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i12);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i13);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i14);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i15);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i16);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i17);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i18);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i19);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i20);

                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i21);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i22);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i23);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i24);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i25);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i26);
                    t_meshdata.m_bulkelmt_connectivity[e-1].push_back(i27);
                   
                    t_meshdata.m_bulkelmt_volume[e-1]=2.0*dx*2.0*dy*2.0*dz;

                    if(i==1){
                        // for left bc elements
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i1);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i5);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i8);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i4);

                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i17);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i16);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i20);
                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i12);

                        leftconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i21);

                        leftnodes.push_back(i1);
                        leftnodes.push_back(i5);
                        leftnodes.push_back(i8);
                        leftnodes.push_back(i4);
                        leftnodes.push_back(i17);
                        leftnodes.push_back(i16);
                        leftnodes.push_back(i20);
                        leftnodes.push_back(i12);
                        leftnodes.push_back(i21);
                    }
                    if(i==t_meshdata.m_nx){
                        // for right bc elements
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i2);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i3);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i7);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i6);

                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i10);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i19);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i14);
                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i18);

                        rightconn[(k-1)*t_meshdata.m_ny+j-1].push_back(i22);

                        rightnodes.push_back(i2);
                        rightnodes.push_back(i3);
                        rightnodes.push_back(i7);
                        rightnodes.push_back(i6);
                        rightnodes.push_back(i10);
                        rightnodes.push_back(i19);
                        rightnodes.push_back(i14);
                        rightnodes.push_back(i18);
                        rightnodes.push_back(i22);
                    }
                    if(j==1){
                        // for bottom bc elements
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i2);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i5);

                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i9 );
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i18);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i13);
                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i17);

                        bottomconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i23);

                        bottomnodes.push_back(i1);
                        bottomnodes.push_back(i2);
                        bottomnodes.push_back(i6);
                        bottomnodes.push_back(i5);
                        bottomnodes.push_back(i9 );
                        bottomnodes.push_back(i18);
                        bottomnodes.push_back(i13);
                        bottomnodes.push_back(i17);
                        bottomnodes.push_back(i23);
                    }
                    if(j==t_meshdata.m_ny){
                        // for top bc elements
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i8);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i3);

                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i20);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i15);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i19);
                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i11);

                        topconn[(k-1)*t_meshdata.m_nx+i-1].push_back(i24);

                        topnodes.push_back(i4);
                        topnodes.push_back(i8);
                        topnodes.push_back(i7);
                        topnodes.push_back(i3);
                        topnodes.push_back(i20);
                        topnodes.push_back(i15);
                        topnodes.push_back(i19);
                        topnodes.push_back(i11);
                        topnodes.push_back(i24);
                    }
                    if(k==1){
                        // for back bc elements
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i1);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i4);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i3);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i2);

                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i12);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i11);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i10);
                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i9 );

                        backconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i25);

                        backnodes.push_back(i1);
                        backnodes.push_back(i4);
                        backnodes.push_back(i3);
                        backnodes.push_back(i2);
                        backnodes.push_back(i12);
                        backnodes.push_back(i11);
                        backnodes.push_back(i10);
                        backnodes.push_back(i9 );
                        backnodes.push_back(i25);
                    }
                    if(k==t_meshdata.m_nz){
                        // for front bc elements
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i5);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i6);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i7);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i8);

                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i13);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i14);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i15);
                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i16);

                        frontconn[(j-1)*t_meshdata.m_nx+i-1].push_back(i26);

                        frontnodes.push_back(i5);
                        frontnodes.push_back(i6);
                        frontnodes.push_back(i7);
                        frontnodes.push_back(i8);
                        frontnodes.push_back(i13);
                        frontnodes.push_back(i14);
                        frontnodes.push_back(i15);
                        frontnodes.push_back(i16);
                        frontnodes.push_back(i26);
                    }
                }
            }
        }
    }//end-of-hex27-generation
    else{
        MessagePrinter::printErrorTxt("unsupported mesh generation in 2d case. Currently, only quad4,quad8, and quad9 are supported.");
        m_mesh_generated=false;
        return m_mesh_generated;
    }

    // remove the duplicate node id, for nodal type set, you don't need these duplicate node ids!
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
    // for backnodes
    sort(backnodes.begin(),backnodes.end());
    backnodes.erase(unique(backnodes.begin(),backnodes.end()),backnodes.end());
    // for frontnodes
    sort(frontnodes.begin(),frontnodes.end());
    frontnodes.erase(unique(frontnodes.begin(),frontnodes.end()),frontnodes.end());

    //***********************************************
    // set up the physical group information
    //***********************************************
    t_meshdata.m_phygroups=1+6;// two bc point+bulk elements
    t_meshdata.m_phygroup_dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phynamevec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyidvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_elmtnumvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_nodesnumperelmtvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2dimvec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_name2elmtconnvec.resize(t_meshdata.m_phygroups);
    // now we can add physical group information
    t_meshdata.m_phygroup_dimvec[0]=2;// for left    surface
    t_meshdata.m_phygroup_dimvec[1]=2;// for right   surface
    t_meshdata.m_phygroup_dimvec[2]=2;// for bottom  surface
    t_meshdata.m_phygroup_dimvec[3]=2;// for top     surface
    t_meshdata.m_phygroup_dimvec[4]=2;// for back    surface
    t_meshdata.m_phygroup_dimvec[5]=2;// for front   surface
    t_meshdata.m_phygroup_dimvec[6]=3;// for bulk element (3d-volume)
    //*** element number for each phy group
    t_meshdata.m_phygroup_elmtnumvec[0]=t_meshdata.m_ny*t_meshdata.m_nz;// for left   surface
    t_meshdata.m_phygroup_elmtnumvec[1]=t_meshdata.m_ny*t_meshdata.m_nz;// for right  surface
    t_meshdata.m_phygroup_elmtnumvec[2]=t_meshdata.m_nx*t_meshdata.m_nz;// for bottom surface
    t_meshdata.m_phygroup_elmtnumvec[3]=t_meshdata.m_nx*t_meshdata.m_nz;// for top    surface
    t_meshdata.m_phygroup_elmtnumvec[4]=t_meshdata.m_nx*t_meshdata.m_ny;// for back   surface
    t_meshdata.m_phygroup_elmtnumvec[5]=t_meshdata.m_nx*t_meshdata.m_ny;// for front  surface
    t_meshdata.m_phygroup_elmtnumvec[6]=t_meshdata.m_bulkelmts;
    //*** nodes number per elmt for each phy group
    t_meshdata.m_phygroup_nodesnumperelmtvec[0]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[1]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[2]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[3]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[4]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[5]=t_meshdata.m_nodespersurfaceelmt;
    t_meshdata.m_phygroup_nodesnumperelmtvec[6]=t_meshdata.m_nodesperbulkelmt;
    //*** for phy name vector
    t_meshdata.m_phygroup_phynamevec[0]="left";
    t_meshdata.m_phygroup_phynamevec[1]="right";
    t_meshdata.m_phygroup_phynamevec[2]="bottom";
    t_meshdata.m_phygroup_phynamevec[3]="top";
    t_meshdata.m_phygroup_phynamevec[4]="back";
    t_meshdata.m_phygroup_phynamevec[5]="front";
    t_meshdata.m_phygroup_phynamevec[6]="alldomain";
    //*** for phy id vector
    t_meshdata.m_phygroup_phyidvec[0]=1;
    t_meshdata.m_phygroup_phyidvec[1]=2;
    t_meshdata.m_phygroup_phyidvec[2]=3;
    t_meshdata.m_phygroup_phyidvec[3]=4;
    t_meshdata.m_phygroup_phyidvec[4]=5;
    t_meshdata.m_phygroup_phyidvec[5]=6;
    t_meshdata.m_phygroup_phyidvec[6]=7;
    //*** for phy name 2 dim map
    t_meshdata.m_phygroup_name2dimvec[0]=make_pair("left",     2);
    t_meshdata.m_phygroup_name2dimvec[1]=make_pair("right",    2);
    t_meshdata.m_phygroup_name2dimvec[2]=make_pair("bottom",   2);
    t_meshdata.m_phygroup_name2dimvec[3]=make_pair("top",      2);
    t_meshdata.m_phygroup_name2dimvec[4]=make_pair("back",     2);
    t_meshdata.m_phygroup_name2dimvec[5]=make_pair("front",    2);
    t_meshdata.m_phygroup_name2dimvec[6]=make_pair("alldomain",3);
    //*** for phy name to element connectivity map
    t_meshdata.m_phygroup_name2elmtconnvec[0]=make_pair("left",                               leftconn);
    t_meshdata.m_phygroup_name2elmtconnvec[1]=make_pair("right",                             rightconn);
    t_meshdata.m_phygroup_name2elmtconnvec[2]=make_pair("bottom",                           bottomconn);
    t_meshdata.m_phygroup_name2elmtconnvec[3]=make_pair("top",                                 topconn);
    t_meshdata.m_phygroup_name2elmtconnvec[4]=make_pair("back",                               backconn);
    t_meshdata.m_phygroup_name2elmtconnvec[5]=make_pair("front",                             frontconn);
    t_meshdata.m_phygroup_name2elmtconnvec[6]=make_pair("alldomain",t_meshdata.m_bulkelmt_connectivity);
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
    t_meshdata.m_phygroup_name2phyidvec[4]=make_pair("back",     5);
    t_meshdata.m_phygroup_name2phyidvec[5]=make_pair("front",    6);
    t_meshdata.m_phygroup_name2phyidvec[6]=make_pair("alldomain",7);
    //***
    t_meshdata.m_phygroup_phyid2namevec.resize(t_meshdata.m_phygroups);
    t_meshdata.m_phygroup_phyid2namevec[0]=make_pair(1     ,"left");
    t_meshdata.m_phygroup_phyid2namevec[1]=make_pair(2    ,"right");
    t_meshdata.m_phygroup_phyid2namevec[2]=make_pair(3   ,"bottom");
    t_meshdata.m_phygroup_phyid2namevec[3]=make_pair(4      ,"top");
    t_meshdata.m_phygroup_phyid2namevec[4]=make_pair(5     ,"back");
    t_meshdata.m_phygroup_phyid2namevec[5]=make_pair(6    ,"front");
    t_meshdata.m_phygroup_phyid2namevec[6]=make_pair(7,"alldomain");

    //***********************************************
    // set up the node set physical group information
    //***********************************************
    t_meshdata.m_nodal_phygroups=6;
    t_meshdata.m_nodephygroup_phynamevec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_phyidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2nodeidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2nodeidvec[0]=make_pair("left",    leftnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[1]=make_pair("right",  rightnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[2]=make_pair("bottom",bottomnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[3]=make_pair("top",      topnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[4]=make_pair("back",    backnodes);
    t_meshdata.m_nodephygroup_name2nodeidvec[5]=make_pair("front",  frontnodes);
    //*** phy id vector
    t_meshdata.m_nodephygroup_phyidvec[0]=1000;
    t_meshdata.m_nodephygroup_phyidvec[1]=2000;
    t_meshdata.m_nodephygroup_phyidvec[2]=3000;
    t_meshdata.m_nodephygroup_phyidvec[3]=4000;
    t_meshdata.m_nodephygroup_phyidvec[4]=5000;
    t_meshdata.m_nodephygroup_phyidvec[5]=6000;
    //*** phy name vector
    t_meshdata.m_nodephygroup_phynamevec[0]="left";
    t_meshdata.m_nodephygroup_phynamevec[1]="right";
    t_meshdata.m_nodephygroup_phynamevec[2]="bottom";
    t_meshdata.m_nodephygroup_phynamevec[3]="top";
    t_meshdata.m_nodephygroup_phynamevec[4]="back";
    t_meshdata.m_nodephygroup_phynamevec[5]="front";
    //*** for phy name to id map
    t_meshdata.m_nodephygroup_name2phyidvec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_name2phyidvec[0]=make_pair("left",  1000);
    t_meshdata.m_nodephygroup_name2phyidvec[1]=make_pair("right", 2000);
    t_meshdata.m_nodephygroup_name2phyidvec[2]=make_pair("bottom",3000);
    t_meshdata.m_nodephygroup_name2phyidvec[3]=make_pair("top",   4000);
    t_meshdata.m_nodephygroup_name2phyidvec[4]=make_pair("back",  5000);
    t_meshdata.m_nodephygroup_name2phyidvec[5]=make_pair("front", 6000);
    //*** for phyid to name map
    t_meshdata.m_nodephygroup_phyid2namevec.resize(t_meshdata.m_nodal_phygroups);
    t_meshdata.m_nodephygroup_phyid2namevec[0]=make_pair(1000  ,"left");
    t_meshdata.m_nodephygroup_phyid2namevec[1]=make_pair(2000 ,"right");
    t_meshdata.m_nodephygroup_phyid2namevec[2]=make_pair(3000,"bottom");
    t_meshdata.m_nodephygroup_phyid2namevec[3]=make_pair(4000   ,"top");
    t_meshdata.m_nodephygroup_phyid2namevec[4]=make_pair(5000,  "back");
    t_meshdata.m_nodephygroup_phyid2namevec[5]=make_pair(6000, "front");

    m_mesh_generated=true;

    return m_mesh_generated;

}