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
//+++ Date   : 2022.05.14
//+++ Purpose: implement the general shape function class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFun1D.h"
#include "FE/ShapeFun2D.h"
#include "FE/ShapeFun3D.h"
#include "FE/ShapeFunUser.h"
#include "FE/ShapeFunType.h"

/**
 * This class implement the calculation and data storage of general shape functions
 */
class ShapeFun:public ShapeFun1D,
               public ShapeFun2D,
               public ShapeFun3D,
               public ShapeFunUser{
public:
    /**
     * constructor
     */
    ShapeFun();

    //***************************************************
    //*** for shape function calculation
    //***************************************************
    /**
     * shape function calc for general case, you can set eta=zeta=0 for 1d case, and zeta=0 for 2d case
     * @param xi f\$\xi f\$ for local coordinates
     * @param eta f\$\eta f\$ for local coordinates
     * @param zeta f\$\zeta f\$ for local coordinates
     * @param t_nodes the nodal coordinates of current mesh
     * @param flag boolean, true for global derivatives, false for local one
     */
    void calc(const double &xi,const double &eta,const double &zeta,const Nodes &t_nodes,const bool &flag=true);

    //***************************************************
    //*** for general settings
    //***************************************************
    /**
     * set the type of shape function, in most cases, it should be 'default' !
     * @param shptype the type of shape function
     */
    void setShapeFunType(const ShapeFunType &shptype){
        m_shp_type=shptype;
        if(shptype!=ShapeFunType::DEFAULT){
            MessagePrinter::printWarningTxt("you are trying to use a user-defined shape function, you must set up the number of shape functions correctly in your code");
        }
    }
    /**
     * set the number of shape functions, it should be used for user-function case!
     * @param nfuns the number of shape functions
     */
    void setShapeFunNums(const int &nfuns){m_funs=nfuns;}
    /**
     * set the type of mesh,
     * @param meshtype the type of mesh
     */
    void setMeshType(const MeshType &meshtype){m_mesh_type=meshtype;}
    /**
     * initialize the shape function class and allocate memory, before you call this function, the dim, type, and funs num 
     * must be set up correctly.
     */
    void init();
    /**
     * get the jacobian determite reference
     */
    inline double& getJacDet(){
        return m_jacdet;
    }
    /**
     * get the i-th shape function value reference
     * @param i integer for i-th shape function
     */
    inline double& shape_value(const int &i){
        if(i<1||i>m_funs){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+"  is out of range in shape_value of ShapeFun.h");
            MessagePrinter::exitAsFem();
        }
        return m_shpvals[i-1];
    }
    /**
     * get the i-th shape function value
     * @param i integer for i-th shape function
     */
    inline double shape_value(const int &i)const{
        if(i<1||i>m_funs){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+"  is out of range in shape_value of ShapeFun.h");
            MessagePrinter::exitAsFem();
        }
        return m_shpvals[i-1];
    }
    /**
     * get the reference of i-th shape function's gradient
     * @param i integer for i-th shape function
     */
    inline Vector3d& shape_grad(const int &i){
        if(i<1||i>m_funs){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+"  is out of range in shape_grad of ShapeFun.h");
            MessagePrinter::exitAsFem();
        }
        return m_shpgrads[i-1];
    }
    /**
     * get the value of i-th shape function's gradient
     * @param i integer for i-th shape function
     */
    inline Vector3d shape_grad(const int &i)const{
        if(i<1||i>m_funs){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+"  is out of range in shape_grad of ShapeFun.h");
            MessagePrinter::exitAsFem();
        }
        return m_shpgrads[i-1];
    }
    /**
     * release the allocated memory
     */
    void releaseMemory();
    

private:
    MeshType m_mesh_type;/**< the type of mesh */
    ShapeFunType m_shp_type;/**< the type of shape function */
    int m_dim;/**< dimension of shape function, i.e., 1d, 2d, and 3d. */
    int m_funs;/**< number of shape functions */
    vector<double> m_shpvals;/**< vector for the shape function values */
    vector<Vector3d> m_shpgrads;/**< vector for the shape function derivatives */
    double m_jacdet;/**< the jacobian's determite */

};