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
//+++ Date   : 2022.07.22
//+++ Purpose: Implement least square projection from integration
//+++          points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/LeastSquareProjection.h"

LeastSquareProjection::LeastSquareProjection(){

}
//******************************************************
//*** for local projection action
//******************************************************
void LeastSquareProjection::localProjectionAction(const int &nodesnum,
                                                  const vector<int> &elconn,
                                                  const double &detjac,
                                                  const ShapeFun &shp,
                                                  const MaterialsContainer &mate,
                                                  ProjectionData &data){
    
    if(data.m_scalarmate_num){
        projectLocalScalarMate2Global(nodesnum,elconn,detjac,shp,mate,data.m_scalarmate_num,data.m_scalarmate_namelist,data.m_proj_scalarmate_vec);
    }
    if(data.m_vectormate_num){
        projectLocalVectorMate2Global(nodesnum,elconn,detjac,shp,mate,data.m_vectormate_num,data.m_vectormate_namelist,data.m_proj_vectormate_vec);
    }
    if(data.m_rank2mate_num){
        projectLocalRank2Mate2Global(nodesnum,elconn,detjac,shp,mate,data.m_rank2mate_num,data.m_rank2mate_namelist,data.m_proj_rank2mate_vec);
    }
    if(data.m_rank4mate_num){
        projectLocalRank4Mate2Global(nodesnum,elconn,detjac,shp,mate,data.m_rank4mate_num,data.m_rank4mate_namelist,data.m_proj_rank4mate_vec);
    }
}
//******************************************************
void LeastSquareProjection::projectLocalScalarMate2Global(const int &nodesnum,
                                       const vector<int> &elconn,
                                       const double &detjac,
                                       const ShapeFun &shp,
                                       const MaterialsContainer &mate,
                                       const int &nproj,
                                       const vector<string> &scalar_name,
                                       Vector &scalar_projvec){
    if(nproj!=static_cast<int>(scalar_name.size())){
        MessagePrinter::printErrorTxt("the number of projected scalar material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    for(j=1;j<=nodesnum;j++){
        iInd=elconn[j-1];
        jInd=(iInd-1)*(1+nproj)+0+1;
        w=detjac*shp.shape_value(j);
        scalar_projvec.addValue(jInd,w);
        for(k=1;k<=nproj;k++){
            HasName=false;
            for(const auto &it:mate.getScalarMaterialsCopy()){
                if(it.first==scalar_name[k-1]){
                    HasName=true;
                    jInd=(iInd-1)*(nproj+1)+k+1;
                    scalar_projvec.addValue(jInd,it.second*w);
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::printErrorTxt("can\'t find the scalar material ("+scalar_name[k-1]+") for projection, please check either your input file or your material code");
                MessagePrinter::exitAsFem();
            }
        }
    }
    
}
//******************************************************
void LeastSquareProjection::projectLocalVectorMate2Global(const int &nodesnum,
                                       const vector<int> &elconn,
                                       const double &detjac,
                                       const ShapeFun &shp,
                                       const MaterialsContainer &mate,
                                       const int &nproj,
                                       const vector<string> &vector_name,
                                       Vector &vector_projvec){
    if(nproj!=static_cast<int>(vector_name.size())){
        MessagePrinter::printErrorTxt("the number of projected vector material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    for(j=1;j<=nodesnum;j++){
        iInd=elconn[j-1];
        jInd=(iInd-1)*(1+3*nproj)+0+1;
        w=detjac*shp.shape_value(j);
        vector_projvec.addValue(jInd,w);
        for(k=1;k<=nproj;k++){
            HasName=false;
            for(const auto &it:mate.getVectorMaterialsCopy()){
                if(it.first==vector_name[k-1]){
                    HasName=true;
                    jInd=(iInd-1)*(3*nproj+1)+3*(k-1)+1+1;
                    vector_projvec.addValue(jInd,it.second(1)*w);
                    jInd=(iInd-1)*(3*nproj+1)+3*(k-1)+2+1;
                    vector_projvec.addValue(jInd,it.second(2)*w);
                    jInd=(iInd-1)*(3*nproj+1)+3*(k-1)+3+1;
                    vector_projvec.addValue(jInd,it.second(3)*w);
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::printErrorTxt("can\'t find the vector material ("+vector_name[k-1]+") for projection, please check either your input file or your material code");
                MessagePrinter::exitAsFem();
            }
        }
    }
}
//******************************************************
void LeastSquareProjection::projectLocalRank2Mate2Global(const int &nodesnum,
                                      const vector<int> &elconn,
                                      const double &detjac,
                                      const ShapeFun &shp,
                                      const MaterialsContainer &mate,
                                      const int &nproj,
                                      const vector<string> &rank2_name,
                                      Vector &rank2_projvec){
    if(nproj!=static_cast<int>(rank2_name.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-2 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    for(j=1;j<=nodesnum;j++){
        iInd=elconn[j-1];
        jInd=(iInd-1)*(1+9*nproj)+0+1;
        w=detjac*shp.shape_value(j);
        rank2_projvec.addValue(jInd,w);
        for(k=1;k<=nproj;k++){
            HasName=false;
            for(const auto &it:mate.getRank2MaterialsCopy()){
                if(it.first==rank2_name[k-1]){
                    HasName=true;
                    ii=0;
                    for(i1=1;i1<=3;i1++){
                        for(j1=1;j1<=3;j1++){
                            ii+=1;
                            jInd=(iInd-1)*(9*nproj+1)+9*(k-1)+ii+1;
                            rank2_projvec.addValue(jInd,it.second(i1,j1)*w);
                        }
                    }
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::printErrorTxt("can\'t find the rank-2 material ("+rank2_name[k-1]+") for projection, please check either your input file or your material code");
                MessagePrinter::exitAsFem();
            }
        }
    }
}
//******************************************************
void LeastSquareProjection::projectLocalRank4Mate2Global(const int &nodesnum,
                                      const vector<int> &elconn,
                                      const double &detjac,
                                      const ShapeFun &shp,
                                      const MaterialsContainer &mate,
                                      const int &nproj,
                                      const vector<string> &rank4_name,
                                      Vector &rank4_projvec){
    if(nproj!=static_cast<int>(rank4_name.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-4 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    for(j=1;j<=nodesnum;j++){
        iInd=elconn[j-1];
        jInd=(iInd-1)*(1+36*nproj)+0+1;
        w=detjac*shp.shape_value(j);
        rank4_projvec.addValue(jInd,w);
        for(k=1;k<=nproj;k++){
            HasName=false;
            for(const auto &it:mate.getRank4MaterialsCopy()){
                if(it.first==rank4_name[k-1]){
                    HasName=true;
                    ii=0;
                    for(i1=1;i1<=6;i1++){
                        for(j1=1;j1<=6;j1++){
                            ii+=1;
                            jInd=(iInd-1)*(36*nproj+1)+36*(k-1)+ii+1;
                            rank4_projvec.addValue(jInd,it.second.getVoigtComponent(i1,j1)*w);
                        }
                    }
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::printErrorTxt("can\'t find the rank-4 material ("+rank4_name[k-1]+") for projection, please check either your input file or your material code");
                MessagePrinter::exitAsFem();
            }
        }
    }
}
//******************************************************
//*** for global projection action
//******************************************************
void LeastSquareProjection::globalProjectionAction(const Mesh &mesh,
                                                   ProjectionData &data){
    int i,j,k,iInd,ee,nproj;
    double value,weight,newvalue;
    PetscMPIInt rank,size;

    // assemble all the local projection value
    data.m_proj_scalarmate_vec.assemble();
    data.m_proj_vectormate_vec.assemble();
    data.m_proj_rank2mate_vec.assemble();
    data.m_proj_rank4mate_vec.assemble();


    data.m_proj_scalarmate_vec.makeGhostCopy();
    data.m_proj_vectormate_vec.makeGhostCopy();
    data.m_proj_rank2mate_vec.makeGhostCopy();
    data.m_proj_rank4mate_vec.makeGhostCopy();

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);


    int rankne=mesh.getBulkMeshNodesNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=mesh.getBulkMeshNodesNum();

    data.m_proj_scalarmate_vec.setToZero();
    data.m_proj_vectormate_vec.setToZero();
    data.m_proj_rank2mate_vec.setToZero();
    data.m_proj_rank4mate_vec.setToZero();

    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        //****************************************
        // for scalar materials
        //****************************************
        nproj=data.m_scalarmate_num;
        if(nproj){
            iInd=(i-1)*(1+nproj)+0+1;
            weight=data.m_proj_scalarmate_vec.getIthValueFromGhost(iInd);
            for(j=1;j<=nproj;j++){
                iInd=(i-1)*(nproj+1)+j+1;
                value=data.m_proj_scalarmate_vec.getIthValueFromGhost(iInd);
                if(abs(weight)>1.0e-15){
                    newvalue=value/weight;
                }
                else{
                    newvalue=value;
                }
                data.m_proj_scalarmate_vec.addValue(iInd,newvalue);
            }
        }
        //****************************************
        // for vector materials
        //****************************************
        nproj=data.m_vectormate_num;
        if(nproj){
            iInd=(i-1)*(1+3*nproj)+0+1;
            weight=data.m_proj_vectormate_vec.getIthValueFromGhost(iInd);
            for(j=1;j<=nproj;j++){
                for(k=1;k<=3;k++){
                    iInd=(i-1)*(3*nproj+1)+3*(j-1)+k+1;
                    value=data.m_proj_vectormate_vec.getIthValueFromGhost(iInd);
                    if(abs(weight)>1.0e-15){
                        newvalue=value/weight;
                    }
                    else{
                        newvalue=value;
                    }
                    data.m_proj_vectormate_vec.addValue(iInd,newvalue);
                }
            }
        }
        //****************************************
        // for rank-2 materials
        //****************************************
        nproj=data.m_rank2mate_num;
        if(nproj){
            iInd=(i-1)*(1+9*nproj)+0+1;
            weight=data.m_proj_rank2mate_vec.getIthValueFromGhost(iInd);
            for(j=1;j<=nproj;j++){
                for(k=1;k<=9;k++){
                    iInd=(i-1)*(9*nproj+1)+9*(j-1)+k+1;
                    value=data.m_proj_rank2mate_vec.getIthValueFromGhost(iInd);
                    if(abs(weight)>1.0e-15){
                        newvalue=value/weight;
                    }
                    else{
                        newvalue=value;
                    }
                    data.m_proj_rank2mate_vec.addValue(iInd,newvalue);
                }
            }
        }
        //****************************************
        // for rank-4 materials
        //****************************************
        nproj=data.m_rank4mate_num;
        if(nproj){
            iInd=(i-1)*(1+36*nproj)+0+1;
            weight=data.m_proj_rank4mate_vec.getIthValueFromGhost(iInd);
            for(j=1;j<=nproj;j++){
                for(k=1;k<=36;k++){
                    iInd=(i-1)*(36*nproj+1)+36*(j-1)+k+1;
                    value=data.m_proj_rank4mate_vec.getIthValueFromGhost(iInd);
                    if(abs(weight)>1.0e-15){
                        newvalue=value/weight;
                    }
                    else{
                        newvalue=value;
                    }
                    data.m_proj_rank4mate_vec.addValue(iInd,newvalue);
                }
            }
        }
    }

    data.m_proj_scalarmate_vec.assemble();
    data.m_proj_vectormate_vec.assemble();
    data.m_proj_rank2mate_vec.assemble();
    data.m_proj_rank4mate_vec.assemble();

    data.m_proj_scalarmate_vec.destroyGhostCopy();
    data.m_proj_vectormate_vec.destroyGhostCopy();
    data.m_proj_rank2mate_vec.destroyGhostCopy();
    data.m_proj_rank4mate_vec.destroyGhostCopy();

}