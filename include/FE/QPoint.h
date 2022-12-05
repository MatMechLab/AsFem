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

#pragma once

#include "FE/QPoint1DGenerator.h"
#include "FE/QPoint2DGenerator.h"
#include "FE/QPoint3DGenerator.h"

#include "FE/QPoint1DLobattoGenerator.h"
#include "FE/QPoint2DLobattoGenerator.h"
#include "FE/QPoint3DLobattoGenerator.h"

#include "FE/QPointUser1Generator.h"
#include "FE/QPointUser2Generator.h"
#include "FE/QPointUser3Generator.h"
#include "FE/QPointUser4Generator.h"
#include "FE/QPointUser5Generator.h"

#include "FE/QPointType.h"

/**
 * This class implement the gauss point generation and management
 */
class QPoint:public QPoint1DGenerator,
             public QPoint2DGenerator,
             public QPoint3DGenerator,
             public QPoint1DLobattoGenerator,
             public QPoint2DLobattoGenerator,
             public QPoint3DLobattoGenerator,
             public QPointUser1Generator,
             public QPointUser2Generator,
             public QPointUser3Generator,
             public QPointUser4Generator,
             public QPointUser5Generator{
public:
    /**
     * Constructor
     */
    QPoint();
    /**
     * Constructor
     * @param a right hand side qpoint class
     */
    QPoint(const QPoint &a);

    /**
     * generate qpoints, before calling this, you should setup the meshype, dimension, qpoint order, and qpoin type
     */
    void createQPoints();

    //**********************************************
    //*** general settings
    //**********************************************
    /**
     * set up the order of accuracy
     * @param order integer
     */
    void setOrder(const int &order){m_order=order;}
    /**
     * set up the dimension of current qpoints
     * @param dim integer
     */
    void setDim(const int &dim){m_dim=dim;}
    /**
     * set up the mesh type for current qpoints
     * @param mtype the mesh type enum
     */
    void setMeshType(const MeshType &mtype){m_meshtype=mtype;}
    /**
     * set up the qpoint generation type for current qpoints
     * @param qptype the qpoint type enum
     */
    void setQPointType(const QPointType &qptype){m_qptype=qptype;}
    //**********************************************
    //*** operator overload
    //**********************************************
    /**
     * () access for the qpoint coordinates
     * @param i integer for the index in 1st dim
     * @param j integer for the index in 2nd dim
     */
    inline double operator()(const int &i,const int &j)const{
        if(i<1||i>m_ngp){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(="+to_string(m_ngp)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        if(j<0||j>m_dim){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of range(="+to_string(m_dim)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        return m_coords[(i-1)*(1+m_dim)+j];
    }
    /**
     * () access for the qpoint coordinates
     * @param i integer for the index in 1st dim
     * @param j integer for the index in 2nd dim
     */
    inline double& operator()(const int &i,const int &j){
        if(i<1||i>m_ngp){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(="+to_string(m_ngp)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        if(j<0||j>m_dim){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of range(="+to_string(m_dim)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        return m_coords[(i-1)*(1+m_dim)+j];
    }
    //**********************************************
    //*** general settings
    //**********************************************
    /**
     * get j-th coordinate of i-th qpoints
     * @param i integer for the index in 1st dim
     * @param j integer for the index in 2nd dim
     */
    inline double getIthPointJthCoord(const int &i,const int &j)const{
        if(i<1||i>m_ngp){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(="+to_string(m_ngp)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        if(j<0||j>m_dim){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of range(="+to_string(m_dim)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        return m_coords[(i-1)*(1+m_dim)+j];
    }
    /**
     * get i-th qpoint's weight
     * @param i integer for the index in 1st dim
     */
    inline double getIthWeight(const int &i)const{
        if(i<1||i>m_ngp){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(="+to_string(m_ngp)+") of your qpoints");
            MessagePrinter::exitAsFem();
        }
        return m_coords[(i-1)*(1+m_dim)+0];
    }
    /**
     * get the number of total qpoints
     */
    inline int getQPointsNum()const{return m_ngp;}
    /**
     * get the accuracy order of current qpoints
     */
    inline int getOrder()const{return m_order;}
    /**
     * get the dimension of current qpoints
     */
    inline int getDim()const{return m_dim;}
    /**
     * get the mesh type of current qpoints
     */
    inline MeshType getMeshType()const{return m_meshtype;}
    /**
     * get the qpoint type of current qpoints
     */
    inline QPointType getQPointType()const{return m_qptype;}

    /**
     * print out the qpoints information
     */
    void printQPointsInfo()const;
    /**
     * print out the qpoints information details
     */
    void printQPointsDetailsInfo()const;

    /**
     * release the allocated memory
     */
    void releaseMemory();


private:
    QPointType m_qptype;/**< the qpoint generation type */
    MeshType m_meshtype;/**< the mesh type for setting up qpoints */
    int m_order;/**< accuracy order of gauss point integration */
    int m_ngp;/**< number of total gauss points */
    int m_dim;/**< the dimension of current qpoints, it should be 1, 2, or 3 */
    vector<double> m_coords;/**< the coordinates and weight of all qpoints */
};