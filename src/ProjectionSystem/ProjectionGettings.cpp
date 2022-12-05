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
//+++ Date   : 2022.09.29
//+++ Purpose: Implement the general access to different projected
//+++          materials via its nodal id and mate name
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/ProjectionSystem.h"

double ProjectionSystem::getIthNodeScalarMateViaMateName(const int &nodeid,const string &matename){
    double matevalue;
    bool hasValue;
    int iInd;
    matevalue=0.0;
    hasValue=false;
    for(int i=1;i<=getScalarMaterialNum();i++){
        if(m_data.m_scalarmate_namelist[i-1]==matename){
            iInd=(nodeid-1)*(1+getScalarMaterialNum())+i+1;
            // m_data.m_proj_scalarmate_vec.makeGhostCopy();
            matevalue=m_data.m_proj_scalarmate_vec.getIthValueFromGhost(iInd);
            // m_data.m_proj_scalarmate_vec.destroyGhostCopy();
            hasValue=true;
            break;
        }
    }
    if(!hasValue){
        MessagePrinter::printErrorTxt("can\'t find projected scalar material(="+matename+"),"
                                      "please check your code or your input file");
        MessagePrinter::exitAsFem();
    }
    return matevalue;
}
//*****************************************************************
Vector3d ProjectionSystem::getIthNodeVectorMateViaMateName(const int &nodeid,const string &matename){
    Vector3d matevalue;
    bool hasValue;
    int iInd;
    matevalue=0.0;
    hasValue=false;
    for(int i=1;i<=getVectorMaterialNum();i++){
        if(m_data.m_vectormate_namelist[i-1]==matename){
            // m_data.m_proj_vectormate_vec.makeGhostCopy();
            iInd=(nodeid-1)*(1+getVectorMaterialNum()*3)+3*(i-1)+1+1;
            matevalue(1)=m_data.m_proj_vectormate_vec.getIthValueFromGhost(iInd);

            iInd=(nodeid-1)*(1+getVectorMaterialNum()*3)+3*(i-1)+2+1;
            matevalue(2)=m_data.m_proj_vectormate_vec.getIthValueFromGhost(iInd);

            iInd=(nodeid-1)*(1+getVectorMaterialNum()*3)+3*(i-1)+3+1;
            matevalue(3)=m_data.m_proj_vectormate_vec.getIthValueFromGhost(iInd);

            // m_data.m_proj_vectormate_vec.destroyGhostCopy();
            hasValue=true;
            break;
        }
    }
    if(!hasValue){
        MessagePrinter::printErrorTxt("can\'t find projected vector material(="+matename+"),"
                                      "please check your code or your input file");
        MessagePrinter::exitAsFem();
    }
    return matevalue;
}
//*******************************************************
Rank2Tensor ProjectionSystem::getIthNodeRank2MateViaMateName(const int &nodeid,const string &matename){
    Rank2Tensor matevalue;
    bool hasValue;
    int iInd,k;
    matevalue=0.0;
    hasValue=false;
    for(int i=1;i<=getRank2MaterialNum();i++){
        if(m_data.m_rank2mate_namelist[i-1]==matename){
            // m_data.m_proj_rank2mate_vec.makeGhostCopy();
            k=0;
            for(int ii=1;ii<=3;ii++){
                for(int jj=1;jj<=3;jj++){
                    k+=1;
                    iInd=(nodeid-1)*(1+getRank2MaterialNum()*9)+9*(i-1)+k+1;
                    matevalue(ii,jj)=m_data.m_proj_rank2mate_vec.getIthValueFromGhost(iInd);
                }
            }
            // m_data.m_proj_rank2mate_vec.destroyGhostCopy();
            hasValue=true;
            break;
        }
    }
    if(!hasValue){
        MessagePrinter::printErrorTxt("can\'t find projected rank-2 material(="+matename+"),"
                                      "please check your code or your input file");
        MessagePrinter::exitAsFem();
    }
    return matevalue;
}
//**********************************************************
Rank4Tensor ProjectionSystem::getIthNodeRank4MateViaMateName(const int &nodeid,const string &matename){
    Rank4Tensor matevalue;
    bool hasValue;
    int iInd;
    matevalue=0.0;
    hasValue=false;
    for(int i=1;i<=getRank4MaterialNum();i++){
        if(m_data.m_rank4mate_namelist[i-1]==matename){
            // m_data.m_proj_rank4mate_vec.makeGhostCopy();
            for(int ii=1;ii<=6;ii++){
                for(int jj=1;jj<=6;jj++){
                    iInd=(nodeid-1)*(1+getRank4MaterialNum()*36)+36*(i-1)+(ii-1)*6+jj+1;
                    matevalue.voigtComponent(ii,jj)=m_data.m_proj_rank4mate_vec.getIthValueFromGhost(iInd);
                }
            }
            // m_data.m_proj_rank4mate_vec.destroyGhostCopy();
            hasValue=true;
            break;
        }
    }
    if(!hasValue){
        MessagePrinter::printErrorTxt("can\'t find projected rank-4 material(="+matename+"),"
                                      "please check your code or your input file");
        MessagePrinter::exitAsFem();
    }
    return matevalue;
}