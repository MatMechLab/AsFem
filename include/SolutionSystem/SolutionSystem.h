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
//+++ Date   : 2022.06.13
//+++ Purpose: defines the solution array for FEM analysis, stores
//+++          the necessary results.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

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
    inline int getQPointsNum()const{return m_QpointsNum;}
    /**
     * get the number of bulk elements
     */
    inline int getBulkElmtsNum()const{return m_BulkElmtsNum;}
    /**
     * get the number of bulk elements
     */
    inline int getLocalBulkElmtsNum()const{return m_BulkElmtsNum_Local;}
    /**
     * get the number of dofs
     */
    inline int getDofsNum()const{return m_Dofs;}

    /**
     * get i-th element's j-th point's scalar material of master rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline ScalarMateType getIthElmtJthScalarMaterial(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in scalar materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in scalar material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsScalarMaterials_Total[(i-1)*m_QpointsNum+j-1];
    }
    /**
     * get i-th element's j-th point's scalar material of local rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline ScalarMateType getLocalIthElmtJthScalarMaterial(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum_Local){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for local bulk elements number in scalar materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for local bulk qpoints number in scalar material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsScalarMaterials_Local[(i-1)*m_QpointsNum+j-1];
    }

    /**
     * get i-th element's j-th qpoint's vector material of master rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline VectorMateType getIthElmtJthVectorMaterial(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in vector materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in vector material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsVectorMaterials_Total[(i-1)*m_QpointsNum+j-1];
    }
    /**
     * get i-th element's j-th qpoint's vector material of local rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline VectorMateType getLocalIthElmtJthVectorMaterial(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum_Local){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for local bulk elements number in vector materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for local bulk qpoints number in vector material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsVectorMaterials_Local[(i-1)*m_QpointsNum+j-1];
    }


    /**
     * get i-th element's j-th qpoint's rank-2 material of master rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline Rank2MateType getIthElmtJthRank2Material(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in rank-2 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in rank-2 material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsRank2Materials_Total[(i-1)*m_QpointsNum+j-1];
    }
    /**
     * get i-th element's j-th qpoint's rank-2 material of local rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline Rank2MateType getLocalIthElmtJthRank2Material(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum_Local){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for local bulk elements number in rank-2 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for local bulk qpoints number in rank-2 material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsRank2Materials_Local[(i-1)*m_QpointsNum+j-1];
    }

    /**
     * get i-th element's j-th point's rank-4 material of master rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline Rank4MateType getIthElmtJthRank4Material(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in rank-4 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in rank-4 material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsRank4Materials_Total[(i-1)*m_QpointsNum+j-1];
    }
    /**
     * get i-th element's j-th point's rank-4 material of local rank
     * @param i the i-th element index
     * @param j the j-th qpoint index
     */
    inline Rank4MateType getLocalIthElmtJthRank4Material(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum_Local){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk elements number in rank-4 materials access");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_QpointsNum){
            MessagePrinter::printErrorTxt("j="+to_string(i)+" is out of range for bulk qpoints number in rank-4 material access");
            MessagePrinter::exitAsFem();
        }
        return m_QpointsRank4Materials_Local[(i-1)*m_QpointsNum+j-1];
    }


public:
    int m_Dofs;/**< the size of the vector, which equals to the active dofs number */
    Vector m_Ucurrent;/**< the current solution vector */
    Vector m_Uold;/**< the previous solution vector */
    Vector m_Uolder;/**< the pre-previous solution vector */
    Vector m_Utemp;/**< the intermediate/temporary solution vector */
    Vector m_Ucopy;/**< this is used for arbitray/intermediate usage, i.e. copy. Its different with u_temp! */
    Vector m_V;/**< the 'velocity' solution vector */
    Vector m_A;/**< the 'acceleration' solution vector */

public:
    //*****************************************************************
    //*** for the different materials of each gauss point
    //*****************************************************************
    // for the total one, which is only allocated on master rank
    vector<ScalarMateType> m_QpointsScalarMaterials_Total;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_QpointsVectorMaterials_Total;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_QpointsRank2Materials_Total;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_QpointsRank4Materials_Total;/**< for rank-4 materials of each gauss point */
    // for the local one, owned by each rank
    vector<ScalarMateType> m_QpointsScalarMaterials_Local;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_QpointsVectorMaterials_Local;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_QpointsRank2Materials_Local;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_QpointsRank4Materials_Local;/**< for rank-4 materials of each gauss point */

    //*** for the materials of previous solution
    // for the total one
    vector<ScalarMateType> m_QpointsScalarMaterialsOld_Total;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_QpointsVectorMaterialsOld_Total;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_QpointsRank2MaterialsOld_Total;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_QpointsRank4MaterialsOld_Total;/**< for rank-4 materials of each gauss point */
    // for the local one, owned by each rank
    vector<ScalarMateType> m_QpointsScalarMaterialsOld_Local;/**< for scalar materials of each gauss point */
    vector<VectorMateType> m_QpointsVectorMaterialsOld_Local;/**< for vector materials of each gauss point */
    vector<Rank2MateType>  m_QpointsRank2MaterialsOld_Local;/**< for rank-2 materials of each gauss point */
    vector<Rank4MateType>  m_QpointsRank4MaterialsOld_Local;/**< for rank-4 materials of each gauss point */

private:
    bool m_Allocated;/**< boolean flag for the allocation status */
    int m_BulkElmtsNum;/**< the number of total bulk elements */
    int m_BulkElmtsNum_Local;/**< the number of local bulk elements */
    int m_QpointsNum;/**< for the number of gauss points of each bulk element */

};