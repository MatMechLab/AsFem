//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.25
//+++ Purpose: assemble and project gauss points' quantities to nodal
//+++          point
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalProjectionToGlobal(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                               const map<string,double> &ProjVariables,
                                               const ScalarMateType &ScalarMate,
                                               const VectorMateType &VectorMate,
                                               const Rank2MateType &Rank2Mate,
                                               const Rank4MateType &Rank4Mate,
                                               SolutionSystem &solutionSystem){
    //*** assemble local projected variables
    AssembleLocalProjVariable2Global(nNodes,DetJac,shp,solutionSystem.GetProjNumPerNode(),solutionSystem.GetProjNameVec(),
                                     ProjVariables,solutionSystem._Proj);

    //*** assemble local projected scalar materials to global
    AssembleLocalProjScalarMate2Global(nNodes,DetJac,shp,solutionSystem.GetScalarMateProjNumPerNode(),
                                       solutionSystem.GetScalarMateNameVec(),ScalarMate,solutionSystem._ProjScalarMate);

    //*** assemble local projected vector materials to global
    AssembleLocalProjVectorMate2Global(nNodes,DetJac,shp,solutionSystem.GetVectorMateProjNumPerNode(),
                                       solutionSystem.GetVectorMateNameVec(),VectorMate,solutionSystem._ProjVectorMate);

    //*** assemble local projected rank-2 materials to global
    AssembleLocalProjRank2Mate2Global(nNodes,DetJac,shp,solutionSystem.GetRank2MateProjNumPerNode(),
                                      solutionSystem.GetRank2MateNameVec(),Rank2Mate,solutionSystem._ProjRank2Mate);

    //*** assemble local projected rank-4 materials to global
    AssembleLocalProjRank4Mate2Global(nNodes,DetJac,shp,solutionSystem.GetRank4MateProjNumPerNode(),
                                      solutionSystem.GetRank4MateNameVec(),Rank4Mate,solutionSystem._ProjRank4Mate);
}
//******************************************************
//@fun: here we assemble the local projected variables to global
void FESystem::AssembleLocalProjVariable2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                                const int &nProj,vector<string> ProjNameVec,
                                                const map<string,double> &elProj,Vec &ProjVec){
    double w;
    int j,k,jInd,iInd;
    bool HasName;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(nProj+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=nProj;k++){
            HasName=false;
            for(const auto &it:elProj){
                if(it.first==ProjNameVec[k-1]){
                    HasName=true;
                    jInd=iInd*(nProj+1)+k;
                    VecSetValue(ProjVec,jInd,w*it.second,ADD_VALUES);
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::PrintErrorTxt("can not find the projected variable(name="
                                               +ProjNameVec[k-1]
                                               +"), please check the 'name=' option in your [projection] block "
                                                "or the gpProj name in your related UEL");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}
//******************************************************
void FESystem::AssembleLocalProjScalarMate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                                  const int &nProj,vector<string> ProjNameVec,
                                                  const ScalarMateType &ScalarMate,Vec &ProjVec){
    double w;
    int j,k,jInd,iInd;
    bool HasName;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(nProj+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=nProj;k++){
            HasName=false;
            for(const auto &it:ScalarMate){
                if(it.first==ProjNameVec[k-1]){
                    HasName=true;
                    jInd=iInd*(nProj+1)+k;
                    VecSetValue(ProjVec,jInd,w*it.second,ADD_VALUES);
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::PrintErrorTxt("can not find the projected scalar material(name="
                                              +ProjNameVec[k-1]
                                              +"), please check the 'scalarmate=' option in your [projection] block "
                                               "or the '_ScalarMaterials' name in your related UMAT");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}
//******************************************************
void FESystem::AssembleLocalProjVectorMate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                                  const int &nProj,vector<string> ProjNameVec,
                                                  const VectorMateType &VectorMate,Vec &ProjVec){
    double w;
    int j,k,jInd,iInd;
    bool HasName;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(nProj*3+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=nProj;k++){
            HasName=false;
            for(const auto &it:VectorMate){
                if(it.first==ProjNameVec[k-1]){
                    HasName=true;
                    // for first component
                    jInd=iInd*(nProj*3+1)+3*(k-1)+1;
                    VecSetValue(ProjVec,jInd,w*it.second(1),ADD_VALUES);
                    // for second component
                    jInd=iInd*(nProj*3+1)+3*(k-1)+2;
                    VecSetValue(ProjVec,jInd,w*it.second(2),ADD_VALUES);
                    // for third component
                    jInd=iInd*(nProj*3+1)+3*(k-1)+3;
                    VecSetValue(ProjVec,jInd,w*it.second(3),ADD_VALUES);
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::PrintErrorTxt("can not find the projected vector material(name="
                                              +ProjNameVec[k-1]
                                              +"), please check the 'vectormate=' option in your [projection] block "
                                               "or the '_VectorMaterials' name in your related UMAT");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}
//******************************************************
void FESystem::AssembleLocalProjRank2Mate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                                 const int &nProj,vector<string> ProjNameVec,
                                                 const Rank2MateType &Rank2Mate,Vec &ProjVec){
    double w;
    int i1,j1,ii,j,k,jInd,iInd;
    bool HasName;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(nProj*9+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=nProj;k++){
            HasName=false;
            for(const auto &it:Rank2Mate){
                if(it.first==ProjNameVec[k-1]){
                    HasName=true;
                    ii=0;
                    for(i1=1;i1<=3;i1++){
                        for(j1=1;j1<=3;j1++){
                            ii+=1;
                            jInd=iInd*(nProj*9+1)+9*(k-1)+ii;
                            VecSetValue(ProjVec,jInd,w*it.second(i1,j1),ADD_VALUES);
                        }
                    }
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::PrintErrorTxt("can not find the projected rank-2 material(name="
                                              +ProjNameVec[k-1]
                                              +"), please check the 'rank2mate=' option in your [projection] block "
                                               "or the '_Rank2Materials' name in your related UMAT");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}
//******************************************************
void FESystem::AssembleLocalProjRank4Mate2Global(const int &nNodes,const double &DetJac,const ShapeFun &shp,
                                                 const int &nProj,vector<string> ProjNameVec,
                                                 const Rank4MateType &Rank4Mate,Vec &ProjVec){
    double w;
    int j,k,jInd,iInd;
    int i1,j1;
    bool HasName;
    for(j=1;j<=nNodes;j++){
        iInd=_elConn[j-1]-1;
        jInd=iInd*(nProj*36+1)+0;
        w=DetJac*shp.shape_value(j);
        VecSetValue(ProjVec,jInd,w,ADD_VALUES);
        for(k=1;k<=nProj;k++){
            HasName=false;
            for(const auto &it:Rank4Mate){
                if(it.first==ProjNameVec[k-1]){
                    HasName=true;
                    // for first component
                    for(i1=1;i1<=6;i1++){
                        for(j1=1;j1<=6;j1++){
                            jInd=iInd*(nProj*36+1)+36*(k-1)+(i1-1)*6+j1;
                            VecSetValue(ProjVec,jInd,w*it.second.VoigtIJcomponent(i1,j1),ADD_VALUES);
                        }
                    }      
                    break;
                }
            }
            if(!HasName){
                MessagePrinter::PrintErrorTxt("can not find the projected rank-4 material(name="
                                              +ProjNameVec[k-1]
                                              +"), please check the 'rank4mate=' option in your [projection] block "
                                               "or the '_Rank4Materials' name in your related UMAT");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}
//******************************************************
void FESystem::Projection(const int &nTotalNodes,SolutionSystem &solutionSystem){
    int i,j,k,iInd,ee,nProj;
    double value,weight,newvalue;

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);

    int rankne=nTotalNodes/_size;
    int eStart=_rank*rankne;
    int eEnd=(_rank+1)*rankne;
    if(_rank==_size-1) eEnd=nTotalNodes;

    // for projected variables
    VecAssemblyBegin(solutionSystem._Proj);
    VecAssemblyEnd(solutionSystem._Proj);
    VecScatterCreateToAll(solutionSystem._Proj,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecSet(solutionSystem._Proj,0.0);
    nProj=solutionSystem.GetProjNumPerNode();
    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nProj)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        VecSetValue(solutionSystem._Proj,iInd,weight,ADD_VALUES);
        for(j=1;j<=nProj;j++){
            iInd=(i-1)*(nProj+1)+j;
            VecGetValues(_ProjSeq,1,&iInd,&value);
            if(abs(weight)>1.0e-15){
                newvalue=value/weight;
            }
            else{
                newvalue=value;
            }
            VecSetValue(solutionSystem._Proj,iInd,newvalue,ADD_VALUES);
        }
    }
    VecAssemblyBegin(solutionSystem._Proj);
    VecAssemblyEnd(solutionSystem._Proj);
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);

    // for projected scalar materials
    VecAssemblyBegin(solutionSystem._ProjScalarMate);
    VecAssemblyEnd(solutionSystem._ProjScalarMate);
    VecScatterCreateToAll(solutionSystem._ProjScalarMate,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,solutionSystem._ProjScalarMate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,solutionSystem._ProjScalarMate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecSet(solutionSystem._ProjScalarMate,0.0);
    nProj=solutionSystem.GetScalarMateProjNumPerNode();
    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nProj)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        VecSetValue(solutionSystem._ProjScalarMate,iInd,weight,ADD_VALUES);
        for(j=1;j<=nProj;j++){
            iInd=(i-1)*(nProj+1)+j;
            VecGetValues(_ProjSeq,1,&iInd,&value);
            if(abs(weight)>1.0e-15){
                newvalue=value/weight;
            }
            else{
                newvalue=value;
            }
            VecSetValue(solutionSystem._ProjScalarMate,iInd,newvalue,ADD_VALUES);
        }
    }
    VecAssemblyBegin(solutionSystem._ProjScalarMate);
    VecAssemblyEnd(solutionSystem._ProjScalarMate);
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);

    // for projected vector materials
    VecAssemblyBegin(solutionSystem._ProjVectorMate);
    VecAssemblyEnd(solutionSystem._ProjVectorMate);
    VecScatterCreateToAll(solutionSystem._ProjVectorMate,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,solutionSystem._ProjVectorMate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,solutionSystem._ProjVectorMate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecSet(solutionSystem._ProjVectorMate,0.0);
    nProj=solutionSystem.GetVectorMateProjNumPerNode();
    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nProj*3)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        VecSetValue(solutionSystem._ProjVectorMate,iInd,weight,ADD_VALUES);
        for(j=1;j<=nProj;j++){
            for(k=1;k<=3;k++){
                iInd=(i-1)*(nProj*3+1)+3*(j-1)+k;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                if(abs(weight)>1.0e-15){
                    newvalue=value/weight;
                }
                else{
                    newvalue=value;
                }
                VecSetValue(solutionSystem._ProjVectorMate,iInd,newvalue,ADD_VALUES);
            }
        }
    }
    VecAssemblyBegin(solutionSystem._ProjVectorMate);
    VecAssemblyEnd(solutionSystem._ProjVectorMate);
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);

    // for projected rank-2 materials
    VecAssemblyBegin(solutionSystem._ProjRank2Mate);
    VecAssemblyEnd(solutionSystem._ProjRank2Mate);
    VecScatterCreateToAll(solutionSystem._ProjRank2Mate,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,solutionSystem._ProjRank2Mate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,solutionSystem._ProjRank2Mate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecSet(solutionSystem._ProjRank2Mate,0.0);
    nProj=solutionSystem.GetRank2MateProjNumPerNode();
    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nProj*9)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        VecSetValue(solutionSystem._ProjRank2Mate,iInd,weight,ADD_VALUES);
        for(j=1;j<=nProj;j++){
            for(k=1;k<=9;k++){
                iInd=(i-1)*(nProj*9+1)+9*(j-1)+k;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                if(abs(weight)>1.0e-15){
                    newvalue=value/weight;
                }
                else{
                    newvalue=value;
                }
                VecSetValue(solutionSystem._ProjRank2Mate,iInd,newvalue,ADD_VALUES);
            }
        }
    }
    VecAssemblyBegin(solutionSystem._ProjRank2Mate);
    VecAssemblyEnd(solutionSystem._ProjRank2Mate);
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);

    // for projected rank-4 materials
    VecAssemblyBegin(solutionSystem._ProjRank4Mate);
    VecAssemblyEnd(solutionSystem._ProjRank4Mate);
    VecScatterCreateToAll(solutionSystem._ProjRank4Mate,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,solutionSystem._ProjRank4Mate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,solutionSystem._ProjRank4Mate,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecSet(solutionSystem._ProjRank4Mate,0.0);
    nProj=solutionSystem.GetRank4MateProjNumPerNode();
    for(ee=eStart;ee<eEnd;++ee){
        i=ee+1;
        iInd=(i-1)*(1+nProj*36)+0;
        VecGetValues(_ProjSeq,1,&iInd,&weight);
        VecSetValue(solutionSystem._ProjRank4Mate,iInd,weight,ADD_VALUES);
        for(j=1;j<=nProj;j++){
            for(k=1;k<=36;k++){
                iInd=(i-1)*(nProj*36+1)+36*(j-1)+k;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                if(abs(weight)>1.0e-15){
                    newvalue=value/weight;
                }
                else{
                    newvalue=value;
                }
                VecSetValue(solutionSystem._ProjRank4Mate,iInd,newvalue,ADD_VALUES);
            }
        }
    }
    VecAssemblyBegin(solutionSystem._ProjRank4Mate);
    VecAssemblyEnd(solutionSystem._ProjRank4Mate);
    VecScatterDestroy(&_scatterproj);
    VecDestroy(&_ProjSeq);
}
