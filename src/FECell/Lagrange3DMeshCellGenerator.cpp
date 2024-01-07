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
//+++ Date   : 2024.01.06
//+++ Purpose: the 3D lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Lagrange3DMeshCellGenerator.h"

Lagrange3DMeshCellGenerator::Lagrange3DMeshCellGenerator(){
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
Lagrange3DMeshCellGenerator::~Lagrange3DMeshCellGenerator(){
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
//*********************************************************
bool Lagrange3DMeshCellGenerator::generateFECell(const MeshType &t_meshtype,FECellData &t_celldata){
    int rank,size;
    t_celldata.MaxDim=3;
    t_celldata.MinDim=2;

    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank==0){
        // only create mesh on the master rank, then distributed them into different ranks !
        double dx,dy,dz;
        int i,j,k,kk,e;
        int i1,i2,i3,i4,i5,i6,i7,i8,i9;
        int i10,i11,i12,i13,i14,i15,i16,i17,i18,i19;
        int i20,i21,i22,i23,i24,i25,i26,i27;
        vector<int> tempconn;
        vector<SingleMeshCell> TotalMeshCell;
        if(t_meshtype==MeshType::HEX8){
            // generate the mesh cell for hex8 mesh
            t_celldata.MeshOrder=1;
            t_celldata.BulkMeshTypeName="hex8";
            dx=(t_celldata.Xmax-t_celldata.Xmin)/t_celldata.Nx;
            dy=(t_celldata.Ymax-t_celldata.Ymin)/t_celldata.Ny;
            dz=(t_celldata.Zmax-t_celldata.Zmin)/t_celldata.Nz;
            
            t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
            t_celldata.LineElmtsNum=0;
            t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                      +t_celldata.Nx*t_celldata.Nz
                                      +t_celldata.Ny*t_celldata.Nz);
            t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                               +t_celldata.SurfElmtsNum
                               +t_celldata.LineElmtsNum;
                               
            t_celldata.NodesNum=(t_celldata.Nx+1)*(t_celldata.Ny+1)*(t_celldata.Nz+1);
            t_celldata.NodesNumPerBulkElmt=8;
            t_celldata.NodesNumPerSurfElmt=4;
            t_celldata.NodesNumPerLineElmt=0;
            
            t_celldata.BulkElmtVTKCellType=12;
            t_celldata.BulkElmtMeshType=MeshType::HEX8;

            t_celldata.LineElmtVTKCellType=3;
            t_celldata.LineElmtMeshType=MeshType::EDGE2;
            
            t_celldata.SurfElmtVTKCellType=9;
            t_celldata.SurfElmtMeshType=MeshType::QUAD4;

            vector<double> nodecoords;
            nodecoords.resize(t_celldata.NodesNum*3,0.0);
            leftnodes.clear();rightnodes.clear();
            bottomnodes.clear();rightnodes.clear();
            backnodes.clear();frontnodes.clear();
            for(k=1;k<=t_celldata.Nz+1;k++){
                for(j=1;j<=t_celldata.Ny+1;j++){
                    for(i=1;i<=t_celldata.Nx+1;i++){
                        kk=(j-1)*(t_celldata.Nx+1)+i+(k-1)*(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*dz;
                        
                        if(i==1){
                            // for left side nodes
                            leftnodes.push_back(kk);// global id, start from 1
                        }
                        if(i==t_celldata.Nx+1){
                            // for right side nodes
                            rightnodes.push_back(kk);
                        }
                        if(j==1){
                            // for bottom side nodes
                            bottomnodes.push_back(kk);
                        }
                        if(j==t_celldata.Ny+1){
                            // for top side nodes
                            topnodes.push_back(kk);
                        }
                        if(k==1){
                            // for back side nodes
                            backnodes.push_back(kk);
                        }
                        if(k==t_celldata.Nz+1){
                            // for front side nodes
                            frontnodes.push_back(kk);
                        }
                    }
                }
            }// end-of-node-generation

            // for the connectivity information of bulk elements
            t_celldata.MeshCell_Global.resize(t_celldata.BulkElmtsNum);
            leftconn.resize(t_celldata.Ny*t_celldata.Nz);
            rightconn.resize(t_celldata.Ny*t_celldata.Nz);
            //
            bottomconn.resize(t_celldata.Nx*t_celldata.Nz);
            topconn.resize(t_celldata.Nx*t_celldata.Nz);
            //
            backconn.resize(t_celldata.Nx*t_celldata.Ny);
            frontconn.resize(t_celldata.Nx*t_celldata.Ny);
            
            leftnodes.clear();rightnodes.clear();
            bottomnodes.clear();rightnodes.clear();
            backnodes.clear();frontnodes.clear();
            
            tempconn.clear();
            
            for(k=1;k<=t_celldata.Nz;k++){
                for(j=1;j<=t_celldata.Ny;j++){
                    for(i=1;i<=t_celldata.Nx;i++){
                        e=(j-1)*t_celldata.Nx+i+(k-1)*t_celldata.Nx*t_celldata.Ny;
                        i1=(j-1)*(t_celldata.Nx+1)+i+(k-1)*(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        i2=i1+1;
                        i3=i2+t_celldata.Nx+1;
                        i4=i3-1;
                        i5=i1+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        i6=i2+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        i7=i3+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        i8=i4+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        
                        tempconn.push_back(e);

                        t_celldata.MeshCell_Total[e-1].Dim=3;
                        t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=8;
                        t_celldata.MeshCell_Total[e-1].VTKCellType=t_celldata.BulkElmtVTKCellType;
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(8);
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0.resize(8);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);
                        // assign node coords
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,1)=nodecoords[(i2-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,2)=nodecoords[(i2-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,3)=nodecoords[(i2-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,1)=nodecoords[(i3-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,2)=nodecoords[(i3-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,1)=nodecoords[(i4-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,2)=nodecoords[(i4-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,3)=nodecoords[(i4-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,1)=nodecoords[(i5-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,2)=nodecoords[(i5-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,3)=nodecoords[(i5-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,1)=nodecoords[(i6-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,2)=nodecoords[(i6-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,3)=nodecoords[(i6-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,1)=nodecoords[(i7-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,2)=nodecoords[(i7-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,1)=nodecoords[(i8-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,2)=nodecoords[(i8-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,3)=nodecoords[(i8-1)*3+3-1];
                        
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;
                        
                        if(i==1){
                            // for left bc elements
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i1);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i5);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i8);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i4);

                            leftnodes.push_back(i1);
                            leftnodes.push_back(i5);
                            leftnodes.push_back(i8);
                            leftnodes.push_back(i4);
                        }
                        if(i==t_celldata.Nx){
                            // for right bc elements
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i2);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i3);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i7);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i6);

                            rightnodes.push_back(i2);
                            rightnodes.push_back(i3);
                            rightnodes.push_back(i7);
                            rightnodes.push_back(i6);
                        }
                        if(j==1){
                            // for bottom bc elements
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i1);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i2);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i6);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i5);

                            bottomnodes.push_back(i1);
                            bottomnodes.push_back(i2);
                            bottomnodes.push_back(i6);
                            bottomnodes.push_back(i5);
                        }
                        if(j==t_celldata.Ny){
                            // for bottom bc elements
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i4);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i8);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i7);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i3);

                            topnodes.push_back(i4);
                            topnodes.push_back(i8);
                            topnodes.push_back(i7);
                            topnodes.push_back(i3);
                        }
                        if(k==1){
                            // for back bc elements
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i1);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i4);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i3);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i2);

                            backnodes.push_back(i1);
                            backnodes.push_back(i4);
                            backnodes.push_back(i3);
                            backnodes.push_back(i2);
                        }
                        if(k==t_celldata.Nz){
                            // for front bc elements
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i5);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i6);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i7);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i8);

                            frontnodes.push_back(i5);
                            frontnodes.push_back(i6);
                            frontnodes.push_back(i7);
                            frontnodes.push_back(i8);
                        }
                    }
                }
            }// end-of-element-generation-loop
        } // end-of-hex8-generation
        else if(t_meshtype==MeshType::HEX20){
            t_celldata.MeshOrder=2;
            t_celldata.BulkMeshTypeName="hex20";

            dx=(t_celldata.Xmax-t_celldata.Xmin)/(2.0*t_celldata.Nx);
            dy=(t_celldata.Ymax-t_celldata.Ymin)/(2.0*t_celldata.Ny);
            dz=(t_celldata.Zmax-t_celldata.Zmin)/(2.0*t_celldata.Nz);
            
            t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
            t_celldata.LineElmtsNum=0;
            t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                      +t_celldata.Nx*t_celldata.Nz
                                      +t_celldata.Ny*t_celldata.Nz);
            t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                               +t_celldata.SurfElmtsNum
                               +t_celldata.LineElmtsNum;
            
            int nLayer1Nodes=(2*t_celldata.Nx+1)*(2*t_celldata.Ny+1)-t_celldata.Nx*t_celldata.Ny;// for norm layer
            int nLayer2Nodes=(t_celldata.Nx+1)*(t_celldata.Ny+1);          // for middle layer

            t_celldata.NodesNum=nLayer1Nodes*(t_celldata.Nz+1)+nLayer2Nodes*t_celldata.Nz;
            t_celldata.NodesNumPerBulkElmt=20;
            t_celldata.NodesNumPerSurfElmt=8;

            t_celldata.BulkElmtVTKCellType=25;

            t_celldata.LineElmtVTKCellType=4;
            t_celldata.LineElmtMeshType=MeshType::EDGE3;

            t_celldata.SurfElmtVTKCellType=23;
            t_celldata.SurfElmtMeshType=MeshType::QUAD9;
            
            // for the coordinates of each node
            vector<double> nodecoords;
            nodecoords.resize(t_celldata.NodesNum*3,0.0);
            leftnodes.clear();rightnodes.clear();
            bottomnodes.clear();rightnodes.clear();
            backnodes.clear();frontnodes.clear();
            int _Nz,_Ny,_Nx;
            _Nz=t_celldata.Nz;_Ny=t_celldata.Ny;_Nx=t_celldata.Nx;
            for(k=1;k<=_Nz;++k){
                // First for normal layer
                for(j=1;j<=_Ny;++j){
                    // for bottom line of each element
                    for(i=1;i<=2*_Nx+1;i++){
                        kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                    }
                    // for middle line of each element
                    for(i=1;i<=_Nx+1;++i){
                        kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy+dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                    }
                }
                // for top line
                j=_Ny+1;
                for(i=1;i<=2*_Nx+1;i++){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                    nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                    nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
                // Then for middle type layer
                for(j=1;j<=_Ny+1;++j){
                    for(i=1;i<=_Nx+1;++i){
                        kk=(j-1)*(_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz+dz;
                    }
                }
            }
            // for the last top layer
            k=_Nz+1;
            for(j=1;j<=_Ny;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;i++){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                    nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                    nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
                // for middle line of each element
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                    nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy+dy;
                    nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
            }
            // for top line
            j=_Ny+1;
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
            }
            // end-of-node-generation
            t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
            leftconn.resize(t_celldata.Ny*t_celldata.Nz);
            rightconn.resize(t_celldata.Ny*t_celldata.Nz);
            //
            bottomconn.resize(t_celldata.Nx*t_celldata.Nz);
            topconn.resize(t_celldata.Nx*t_celldata.Nz);
            //
            backconn.resize(t_celldata.Nx*t_celldata.Ny);
            frontconn.resize(t_celldata.Nx*t_celldata.Ny);
            
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

                        t_celldata.MeshCell_Total[e-1].Dim=3;
                        t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=20;
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(20);
                        t_celldata.MeshCell_Total[e-1].VTKCellType=t_celldata.BulkElmtVTKCellType;

                        t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i9);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i10);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i11);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i12);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i13);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i14);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i15);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i16);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i17);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i18);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i19);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i20);

                        // for element's nodecoords
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,1)=nodecoords[(21-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,2)=nodecoords[(21-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,3)=nodecoords[(21-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,1)=nodecoords[(i3-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,2)=nodecoords[(i3-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,1)=nodecoords[(i4-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,2)=nodecoords[(i4-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,3)=nodecoords[(i4-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,1)=nodecoords[(i5-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,2)=nodecoords[(i5-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,3)=nodecoords[(i5-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,1)=nodecoords[(i6-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,2)=nodecoords[(i6-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,3)=nodecoords[(i6-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,1)=nodecoords[(i7-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,2)=nodecoords[(i7-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,1)=nodecoords[(i8-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,2)=nodecoords[(i8-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,3)=nodecoords[(i8-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,1)=nodecoords[(i9-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,2)=nodecoords[(i9-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,3)=nodecoords[(i9-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,1)=nodecoords[(i10-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,2)=nodecoords[(i10-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,3)=nodecoords[(i10-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,1)=nodecoords[(i11-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,2)=nodecoords[(i11-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,3)=nodecoords[(i11-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,1)=nodecoords[(i12-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,2)=nodecoords[(i12-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,3)=nodecoords[(i12-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,1)=nodecoords[(i3-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,2)=nodecoords[(i3-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,1)=nodecoords[(i14-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,2)=nodecoords[(i14-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,3)=nodecoords[(i14-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,1)=nodecoords[(i15-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,2)=nodecoords[(i15-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,3)=nodecoords[(i15-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,1)=nodecoords[(i16-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,2)=nodecoords[(i16-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,3)=nodecoords[(i16-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,1)=nodecoords[(i17-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,2)=nodecoords[(i17-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,3)=nodecoords[(i17-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,1)=nodecoords[(i18-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,2)=nodecoords[(i18-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,3)=nodecoords[(i18-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,1)=nodecoords[(i19-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,2)=nodecoords[(i19-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,3)=nodecoords[(i19-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,1)=nodecoords[(i20-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,2)=nodecoords[(i20-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,3)=nodecoords[(i20-1)*3+3-1];

                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;
                        
                        if(i==1){
                            // for left bc elements
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i1);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i5);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i8);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i4);

                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i17);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i16);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i20);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i12);

                            leftnodes.push_back(i1);
                            leftnodes.push_back(i5);
                            leftnodes.push_back(i8);
                            leftnodes.push_back(i4);
                            leftnodes.push_back(i17);
                            leftnodes.push_back(i16);
                            leftnodes.push_back(i20);
                            leftnodes.push_back(i12);
                        }
                        if(i==t_celldata.Nx){
                            // for right bc elements
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i2);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i3);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i7);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i6);

                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i10);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i19);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i14);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i18);

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
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i1);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i2);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i6);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i5);

                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i9 );
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i18);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i13);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i17);

                            bottomnodes.push_back(i1);
                            bottomnodes.push_back(i2);
                            bottomnodes.push_back(i6);
                            bottomnodes.push_back(i5);
                            bottomnodes.push_back(i9 );
                            bottomnodes.push_back(i18);
                            bottomnodes.push_back(i13);
                            bottomnodes.push_back(i17);
                        }
                        if(j==t_celldata.Ny){
                            // for top bc elements
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i4);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i8);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i7);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i3);

                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i20);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i15);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i19);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i11);

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
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i1);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i4);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i3);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i2);

                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i12);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i11);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i10);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i9 );

                            backnodes.push_back(i1);
                            backnodes.push_back(i4);
                            backnodes.push_back(i3);
                            backnodes.push_back(i2);
                            backnodes.push_back(i12);
                            backnodes.push_back(i11);
                            backnodes.push_back(i10);
                            backnodes.push_back(i9 );
                        }
                        if(k==t_celldata.Nz){
                            // for front bc elements
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i5);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i6);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i7);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i8);

                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i13);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i14);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i15);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i16);

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
            } // end-of-element-generation-loop
        } // end-of-hex20-generation
        else if(t_meshtype==MeshType::HEX27){
            t_celldata.MeshOrder=2;
            t_celldata.BulkMeshTypeName="hex27";

            dx=(t_celldata.Xmax-t_celldata.Xmin)/(2.0*t_celldata.Nx);
            dy=(t_celldata.Ymax-t_celldata.Ymin)/(2.0*t_celldata.Ny);
            dz=(t_celldata.Zmax-t_celldata.Zmin)/(2.0*t_celldata.Nz);

            t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
            t_celldata.LineElmtsNum=0;
            t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                      +t_celldata.Nx*t_celldata.Nz
                                      +t_celldata.Ny*t_celldata.Nz);
            
            t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                               +t_celldata.SurfElmtsNum
                               +t_celldata.LineElmtsNum;
            int nLayerNodes=(2*t_celldata.Nx+1)*(2*t_celldata.Ny+1);// for norm layer
            
            t_celldata.NodesNum=(2*t_celldata.Nz+1)*nLayerNodes;
            t_celldata.NodesNumPerBulkElmt=27;
            t_celldata.BulkElmtVTKCellType=29;
            
            t_celldata.NodesNumPerSurfElmt=9;
            t_celldata.SurfElmtVTKCellType=23;
            t_celldata.SurfElmtMeshType=MeshType::QUAD9;
            
            t_celldata.LineElmtVTKCellType=4;
            t_celldata.LineElmtMeshType=MeshType::EDGE3;
            
            // for the coordinates of each node
            vector<double> nodecoords;
            nodecoords.resize(t_celldata.NodesNum*3,0.0);
            
            leftnodes.clear();rightnodes.clear();
            bottomnodes.clear();rightnodes.clear();
            backnodes.clear();frontnodes.clear();
            int _Nz,_Ny,_Nx;
            _Nz=t_celldata.Nz;_Ny=t_celldata.Ny;_Nx=t_celldata.Nx;
            for(k=1;k<=_Nz;++k){
                // For first layer
                for(j=1;j<=2*_Ny+1;++j){
                    // for bottom line of each element
                    for(i=1;i<=2*_Nx+1;++i){
                        kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                    }
                }
                // Then for second layer
                for(j=1;j<=2*_Ny+1;++j){
                    for(i=1;i<=2*_Nx+1;++i){
                        kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes+nLayerNodes;
                        nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                        nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*dy;
                        nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz+dz;
                    }
                }
            }
            // for the last top layer
            k=_Nz+1;
            for(j=1;j<=2*_Ny+1;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                    nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                    nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*dy;
                    nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
            }
            // for the connectivity information of bulk elements
            t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
            leftconn.resize(t_celldata.Ny*t_celldata.Nz);
            rightconn.resize(t_celldata.Ny*t_celldata.Nz);
            //
            bottomconn.resize(t_celldata.Nx*t_celldata.Nz);
            topconn.resize(t_celldata.Nx*t_celldata.Nz);
            //
            backconn.resize(t_celldata.Nx*t_celldata.Ny);
            frontconn.resize(t_celldata.Nx*t_celldata.Ny);
            
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

                        t_celldata.MeshCell_Total[e-1].Dim=3;
                        t_celldata.MeshCell_Total[e-1].VTKCellType=t_celldata.BulkElmtVTKCellType;
                        t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=27;
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(27);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i9);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i10);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i11);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i12);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i13);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i14);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i15);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i16);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i17);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i18);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i19);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i20);

                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i21);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i22);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i23);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i24);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i25);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i26);
                        t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i27);

                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,1)=nodecoords[(i2-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,2)=nodecoords[(i2-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,3)=nodecoords[(i2-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,1)=nodecoords[(i3-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,2)=nodecoords[(i3-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,1)=nodecoords[(i4-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,2)=nodecoords[(i4-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,3)=nodecoords[(i4-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,1)=nodecoords[(i5-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,2)=nodecoords[(i5-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,3)=nodecoords[(i5-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,1)=nodecoords[(i6-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,2)=nodecoords[(i6-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,3)=nodecoords[(i6-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,1)=nodecoords[(i7-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,2)=nodecoords[(i7-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,1)=nodecoords[(i8-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,2)=nodecoords[(i8-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,3)=nodecoords[(i8-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,1)=nodecoords[(i9-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,2)=nodecoords[(i9-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,3)=nodecoords[(i9-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,1)=nodecoords[(i10-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,2)=nodecoords[(i10-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,3)=nodecoords[(i10-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,1)=nodecoords[(i11-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,2)=nodecoords[(i11-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,3)=nodecoords[(i11-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,1)=nodecoords[(i12-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,2)=nodecoords[(i12-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,3)=nodecoords[(i12-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,1)=nodecoords[(i13-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,2)=nodecoords[(i13-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,3)=nodecoords[(i13-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,1)=nodecoords[(i14-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,2)=nodecoords[(i14-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,3)=nodecoords[(i14-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,1)=nodecoords[(i15-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,2)=nodecoords[(i15-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,3)=nodecoords[(i15-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,1)=nodecoords[(i16-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,2)=nodecoords[(i16-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,3)=nodecoords[(i16-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,1)=nodecoords[(i17-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,2)=nodecoords[(i17-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,3)=nodecoords[(i17-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,1)=nodecoords[(i18-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,2)=nodecoords[(i18-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,3)=nodecoords[(i18-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,1)=nodecoords[(i19-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,2)=nodecoords[(i19-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,3)=nodecoords[(i19-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,1)=nodecoords[(i20-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,2)=nodecoords[(i20-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,3)=nodecoords[(i20-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(21,1)=nodecoords[(i21-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(21,2)=nodecoords[(i21-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(21,3)=nodecoords[(i21-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(22,1)=nodecoords[(i22-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(22,2)=nodecoords[(i22-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(22,3)=nodecoords[(i22-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(23,1)=nodecoords[(i23-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(23,2)=nodecoords[(i23-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(23,3)=nodecoords[(i23-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(24,1)=nodecoords[(i24-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(24,2)=nodecoords[(i24-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(24,3)=nodecoords[(i24-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(25,1)=nodecoords[(i25-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(25,2)=nodecoords[(i25-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(25,3)=nodecoords[(i25-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(26,1)=nodecoords[(i26-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(26,2)=nodecoords[(i26-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(26,3)=nodecoords[(i26-1)*3+3-1];
                        //
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(27,1)=nodecoords[(i27-1)*3+1-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(27,2)=nodecoords[(i27-1)*3+2-1];
                        t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(27,3)=nodecoords[(i27-1)*3+3-1];

                        
                        if(i==1){
                            // for left bc elements
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i1);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i5);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i8);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i4);

                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i17);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i16);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i20);
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i12);
                            
                            leftconn[(k-1)*t_celldata.Ny+j-1].push_back(i21);

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
                        if(i==t_celldata.Nx){
                            // for right bc elements
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i2);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i3);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i7);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i6);

                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i10);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i19);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i14);
                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i18);

                            rightconn[(k-1)*t_celldata.Ny+j-1].push_back(i22);

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
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i1);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i2);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i6);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i5);

                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i9 );
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i18);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i13);
                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i17);

                            bottomconn[(k-1)*t_celldata.Nx+i-1].push_back(i23);

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
                        if(j==t_celldata.Ny){
                            // for top bc elements
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i4);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i8);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i7);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i3);

                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i20);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i15);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i19);
                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i11);

                            topconn[(k-1)*t_celldata.Nx+i-1].push_back(i24);

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
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i1);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i4);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i3);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i2);

                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i12);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i11);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i10);
                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i9 );

                            backconn[(j-1)*t_celldata.Nx+i-1].push_back(i25);

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
                        if(k==t_celldata.Nz){
                            // for front bc elements
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i5);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i6);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i7);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i8);

                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i13);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i14);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i15);
                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i16);

                            frontconn[(j-1)*t_celldata.Nx+i-1].push_back(i26);

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
            }// end-of-element-generation

        }//end-of-hex27-mesh

    }// end-of-if(rank==0)
}