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
//+++ Date   : 2022.06.05
//+++ Purpose: implement the general gauss point generation and management
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint.h"

QPoint::QPoint(){
    m_qptype=QPointType::GAUSSLEGENDRE;
    m_meshtype=MeshType::NULLTYPE;
    m_order=0;
    m_ngp=0;
    m_dim=0;
    m_coords.clear();
}

QPoint::QPoint(const QPoint &a){
    m_qptype=a.getQPointType();
    m_meshtype=a.getMeshType();
    m_order=a.getOrder();
    m_dim=a.getDim();
    m_coords=a.m_coords;
}

void QPoint::createQPoints(){
    switch (m_qptype){
        case QPointType::GAUSSLEGENDRE:{
            if(m_dim==1){
                QPoint1DGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
            }
            else if(m_dim==2){
                QPoint2DGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
            }
            else if(m_dim==3){
                QPoint3DGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
            }
            else{
                MessagePrinter::printErrorTxt("order=" + to_string(m_dim) + " is not supported for gauss-legendre rule");
                MessagePrinter::exitAsFem();
                break;
            }
            break;
        }
        case QPointType::GAUSSLOBATTO:{
            if(m_dim==1){
                QPoint1DLobattoGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
                break;
            }
            else if(m_dim==2){
                QPoint2DLobattoGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
                break;
            }
            else if(m_dim==3){
                QPoint3DLobattoGenerator::generateQPoints(m_order, m_meshtype, m_ngp, m_coords);
                break;
            }
            else{
                MessagePrinter::printErrorTxt("order=" + to_string(m_dim) + " is not supported for gauss-lobatto rule");
                MessagePrinter::exitAsFem();
                break;
            }
            break;
        }
        case QPointType::USER1:{
            QPointUser1Generator::generateQPoints(m_dim, m_order, m_meshtype, m_ngp, m_coords);
            break;
        }
        case QPointType::USER2:{
            QPointUser2Generator::generateQPoints(m_dim, m_order, m_meshtype, m_ngp, m_coords);
            break;
        }
        case QPointType::USER3:{
            QPointUser3Generator::generateQPoints(m_dim, m_order, m_meshtype, m_ngp, m_coords);
            break;
        }
        case QPointType::USER4:{
            QPointUser4Generator::generateQPoints(m_dim, m_order, m_meshtype, m_ngp, m_coords);
            break;
        }
        case QPointType::USER5:{
            QPointUser5Generator::generateQPoints(m_dim, m_order, m_meshtype, m_ngp, m_coords);
            break;
        }
        default:{
            MessagePrinter::printErrorTxt("unsupported gauss point generation method, error detected in QPoint.cpp(generateQPoints)");
            MessagePrinter::exitAsFem();
            break;
        }
    }
}
//*********************************************************
void QPoint::printQPointsInfo()const{
    char buff[68];
    string str,qtypename;
    MessagePrinter::printNormalTxt("  qpoint information summary");
    if(m_qptype==QPointType::GAUSSLEGENDRE){
        qtypename="gauss-legendre";
    }
    else if(m_qptype==QPointType::GAUSSLOBATTO){
        qtypename="gauss-lobatto";
    }
    else if(m_qptype==QPointType::USER1){
        qtypename="gauss-user1";
    }
    else if(m_qptype==QPointType::USER2){
        qtypename="gauss-user2";
    }
    else if(m_qptype==QPointType::USER3){
        qtypename="gauss-user3";
    }
    else if(m_qptype==QPointType::USER4){
        qtypename="gauss-user4";
    }
    else if(m_qptype==QPointType::USER5){
        qtypename="gauss-user5";
    }
    snprintf(buff,68,"    dim=%2d,order=%2d,npg=%2d,type=%16s",getDim(),getOrder(),getQPointsNum(),qtypename.c_str());
    str=buff;
    MessagePrinter::printNormalTxt(str);
}

void QPoint::printQPointsDetailsInfo()const{
    char buff[68];
    string str,qtypename;
    MessagePrinter::printNormalTxt("  qpoint information summary");
    if(m_qptype==QPointType::GAUSSLEGENDRE){
        qtypename="gauss-legendre";
    }
    else if(m_qptype==QPointType::GAUSSLOBATTO){
        qtypename="gauss-lobatto";
    }
    else if(m_qptype==QPointType::USER1){
        qtypename="gauss-user1";
    }
    else if(m_qptype==QPointType::USER2){
        qtypename="gauss-user2";
    }
    else if(m_qptype==QPointType::USER3){
        qtypename="gauss-user3";
    }
    else if(m_qptype==QPointType::USER4){
        qtypename="gauss-user4";
    }
    else if(m_qptype==QPointType::USER5){
        qtypename="gauss-user5";
    }
    snprintf(buff,68,"   dim=%2d,order=%2d,npg=%2d,type=%16s",getDim(),getOrder(),getQPointsNum(),qtypename.c_str());
    str=buff;
    MessagePrinter::printNormalTxt(str);
    for(int i=1;i<=getQPointsNum();i++){
        if(getDim()==1){
            snprintf(buff,68,"   xi=%12.4e,w=%11.3e",getIthPointJthCoord(i,1),getIthPointJthCoord(i,0));
        }
        else if(getDim()==2){
            snprintf(buff,68,"   xi=%12.4e,eta=%12.4e,w=%11.3e",getIthPointJthCoord(i,1),getIthPointJthCoord(i,2),getIthPointJthCoord(i,0));
        }
        else if(getDim()==3){
            snprintf(buff,68,"   xi=%12.4e,eta=%12.4e,zeta=%12.4e,w=%11.3e",getIthPointJthCoord(i,1),getIthPointJthCoord(i,2),getIthPointJthCoord(i,3),getIthPointJthCoord(i,0));
        }
        str=buff;
        MessagePrinter::printNormalTxt(str);
    }
}

void QPoint::releaseMemory(){
    m_order=0;
    m_ngp=0;
    m_dim=0;
    m_coords.clear();
}