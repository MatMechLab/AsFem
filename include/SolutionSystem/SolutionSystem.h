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
//+++ Date   : 2022.06.13
//+++ Purpose: defines the solution array for FEM analysis, stores
//+++          the necessary results.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <vector>

#include "MathUtils/Vector.h"
#include "MateSystem/MaterialsName.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"


/**
 * This class defines the general solution vectors for FEM analysis, all the results, history variables, and 
 * material properties should be stored here.
 */
class SolutionSystem{
public:
    /**
     * constructor
     */
    SolutionSystem();
    /**
     * initialize the solution system
     * @param t_dofhandler the dof handler system
     * @param t_fe the fe system
     */
    void init(const DofHandler &t_dofhandler,const FE &t_fe);

    /**
     * update solution vectors
     */
    void updateSolution();
    /**
     * update materials vectors
     */
    void updateMaterialsSolution();

    /**
     * release the allocated memory
     */
    void releaseMemory();

    /**
     * get the number of qpoints of each bulk element
     */
    inline int getQPointsNum()const{return m_qpoints_num;}
    /**
     * get the number of bulk elements
     */
    inline int getBulkElmtsNum()const{return m_bulkelmts_num;}
    /**
     * get the number of dofs
     */
    inline int getDofsNum()const{return m_dofs;}
    /**
     * get i-th element's j-th point's scalar material
     */
    inline ScalarMateType getIthElmtJthScalarMaterial(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in scalar materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_qpoints_num){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in scalar material access");
            MessagePrinter::exitAsFem();
        }
        return m_qpoints_scalarmaterials[(i-1)*m_qpoints_num+j-1];
    }
    /**
     * get i-th element's j-th point's vector material
     */
    inline VectorMateType getIthElmtJthVectorMaterial(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in vector materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_qpoints_num){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in vector material access");
            MessagePrinter::exitAsFem();
        }
        return m_qpoints_vectormaterials[(i-1)*m_qpoints_num+j-1];
    }
    /**
     * get i-th element's j-th point's rank-2 material
     */
    inline Rank2MateType getIthElmtJthRank2Material(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in rank-2 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_qpoints_num){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in rank-2 material access");
            MessagePrinter::exitAsFem();
        }
        return m_qpoints_rank2materials[(i-1)*m_qpoints_num+j-1];
    }
    /**
     * get i-th element's j-th point's rank-4 material
     */
    inline Rank4MateType getIthElmtJthRank4Material(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts_num){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in rank-4 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_qpoints_num){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in rank-4 material access");
            MessagePrinter::exitAsFem();
        }
        return m_qpoints_rank4materials[(i-1)*m_qpoints_num+j-1];
    }


public:
    int m_dofs;/**< the size of the vector, which equals to the active dofs number */
    Vector m_u_current;/**< the current solution vector */
    Vector m_u_old;/**< the previous solution vector */
    Vector m_u_older;/**< the pre-previous solution vector */
    Vector m_u_temp;/**< the intermediate/temporary solution vector */
    Vector m_u_copy;/**< this is used for arbitray/intermediate usage, i.e. copy. Its different with u_temp! */
    Vector m_v;/**< the 'velocity' solution vector */
    Vector m_a;/**< the 'acceleration' solution vector */

public:
    //*****************************************************************
    //*** for the different materials of each gauss point
    //*****************************************************************
    vector<ScalarMateType> m_qpoints_scalarmaterials;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_qpoints_vectormaterials;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_qpoints_rank2materials;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_qpoints_rank4materials;/**< for rank-4 materials of each gauss point */
    //*** for the materials of previous solution
    vector<ScalarMateType> m_qpoints_scalarmaterials_old;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_qpoints_vectormaterials_old;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_qpoints_rank2materials_old;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_qpoints_rank4materials_old;/**< for rank-4 materials of each gauss point */

private:
    bool m_allocated;/**< boolean flag for the allocation status */
    int m_bulkelmts_num;/**< the total number of bulk elements */
    int m_qpoints_num;/**< for the number of gauss points of each bulk element */


};