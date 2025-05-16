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
//+++ Date   : 2022.05.07
//+++ Purpose: this class stores the coordinates of nodes of 
//+++          a single element
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <vector>

#include "Utils/MessagePrinter.h"

using std::vector;

/**
 * the nodes class, which stores the coordinates for multiple nodes, i.e., the nodes of one single element
 */
class Nodes{
public:
    /**
     * constructor
     */
    Nodes();
    /**
     * constructor
     * @param n integer number for the nodes num
     */
    Nodes(const int &n);
    /**
     * constructor
     * @param nodes another nodes class
     */
    Nodes(const Nodes &nodes);
    /**
     * deconstructor
     */
    ~Nodes();
    //************************************************
    //*** general settings
    //************************************************
    /**
     * set the number of total nodes to store, this will resize the whole vector
     * @param n integer for the number of ndoes
     */
    void setNodesNum(const int &n){m_Size=n;m_Coordinates.resize(n*3,0.0);}
    /**
     * resize the coordinates vector, this will resize the whole vector
     * @param n integer for the number of ndoes
     */
    void resize(const int &n){m_Size=n;m_Coordinates.resize(n*3,0.0);}
    /**
     * set all the coordinates to be zero, but keep its capacity/size
     */
    void setToZeros(){
        fill(m_Coordinates.begin(),m_Coordinates.end(),0.0);
    }
    /**
     * clean the vector and memory
     */
    inline void clear(){
        m_Coordinates.clear();m_Size=0;
    }

    //************************************************
    //*** for operators
    //************************************************
    /**
     * () operator for the I-th node's J-th coordinate
     * @param i i-th node
     * @param j j-th coordinate
     */
    double& operator()(const int &i,const int &j){
        if(i<1||i>m_Size){
            MessagePrinter::printErrorTxt("your i index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>3){
            MessagePrinter::printErrorTxt("your j index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        return m_Coordinates[(i-1)*3+j-1];
    }
    /**
     * const () operator for the I-th node's J-th coordinate
     * @param i i-th node
     * @param j j-th coordinate
     */
    double operator()(const int &i,const int &j)const{
        if(i<1||i>m_Size){
            MessagePrinter::printErrorTxt("your i index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>3){
            MessagePrinter::printErrorTxt("your j index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        return m_Coordinates[(i-1)*3+j-1];
    }
    /**
     * = operator for the assignment with a scalar value
     */
    Nodes& operator=(const double &val){
        fill(m_Coordinates.begin(),m_Coordinates.end(),val);
        return *this;
    }
    /**
     * = operator for the assignment with another nodes class
     */
    Nodes& operator=(const Nodes &a){
        if (m_Size!=a.m_Size) {
            m_Size=a.m_Size;
            m_Coordinates.resize(m_Size*3,0.0);
        }
        for(int i=1;i<=a.m_Size*3;i++) m_Coordinates[i-1]=a.m_Coordinates[i-1];
        return *this;
    }
    //************************************************
    //*** for gettings
    //************************************************
    /**
     * get the size/number of total nodes
     */
    inline int getSize()const{return m_Size;}
    /**
     * get the length of coordinate vector
     */
    inline int getLength()const{return m_Size*3;}
    /**
     * get the i-th node's j-th coordinate
     * @param i i-th node
     * @param j j-th coordinate
     */
    inline double getIthNodeJthCoord(const int &i,const int &j)const{
        if(i<1||i>m_Size){
            MessagePrinter::printErrorTxt("your i index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>3){
            MessagePrinter::printErrorTxt("your j index is out of your nodes class's range");
            MessagePrinter::exitAsFem();
        }
        return m_Coordinates[(i-1)*3+j-1];
    }

    /**
     * get the data of the coordinates vector
    */
    inline double* getData(){return m_Coordinates.data();}
    /**
     * get the reference of coordinate vector
    */
    inline vector<double>& getDataRef(){return m_Coordinates;}
    /**
     * get the copy of coordinate vector
    */
    inline vector<double> getCopy()const{return m_Coordinates;}
    

private:
    int m_Size;/**< the size/number of total nodes */
    vector<double> m_Coordinates;/**< vector for the coordinates of each node */

};