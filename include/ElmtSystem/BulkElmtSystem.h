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
//+++ Date   : 2022.05.12
//+++ Purpose: the bulk element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * for AsFem's headers
 */
#include "ElmtSystem/ElmtType.h"
#include "ElmtSystem/ElmtBlock.h"
#include "Utils/MessagePrinter.h"

#include "MathUtils/MatrixXd.h"
#include "MathUtils/VectorXd.h"

#include "FESystem/FECalcType.h"
#include "MateSystem/MaterialsContainer.h"
#include "ElmtSystem/LocalElmtData.h"
#include "FECell/FECell.h"

/**
 * for built-in and user-defined elements
 */
#include "ElmtSystem/PoissonElement.h"
#include "ElmtSystem/DiffusionElement.h"
#include "ElmtSystem/AllenCahnElement.h"
#include "ElmtSystem/MechanicsElement.h"
#include "ElmtSystem/CahnHilliardElement.h"
#include "ElmtSystem/KobayashiElement.h"
#include "ElmtSystem/StressDiffusionElement.h"
#include "ElmtSystem/MieheFractureElement.h"
#include "ElmtSystem/AllenCahnFractureElement.h"
#include "ElmtSystem/StressCahnHilliardElement.h"
#include "ElmtSystem/LaplaceElement.h"
#include "ElmtSystem/ScalarBodySourceElement.h"
#include "ElmtSystem/DiffusionACFractureElement.h"

/**
 * This class implement the elemental calculation for residual, jacobian, projection, etc.
 */
class BulkElmtSystem:public PoissonElement,
                     public DiffusionElement,
                     public AllenCahnElement,
                     public MechanicsElement,
                     public CahnHilliardElement,
                     public KobayashiElement,
                     public StressDiffusionElement,
                     public MieheFractureElement,
                     public AllenCahnFractureElement,
                     public StressCahnHilliardElement,
                     public LaplaceElement,
                     public ScalarBodySourceElement,
                     public DiffusionACFractureElement{
public:
    /**
     * constructor
     */
    BulkElmtSystem();
    /**
     * add element block to the list
     * @param t_elmtblock the element block class read from input file
     */
    void addElmtBlock2List(ElmtBlock t_elmtblock);
    /**
     * initialize the element vector for elmt block info and material block info
     * @param t_fecell the fe cell class
     */
    void init(const FECell &t_fecell);
    /**
     * run the element library for different element(models) calculation
     * here the jacobian and residual will be calculated.
     * @param t_calctype the calculation action type, i.e., form residual or form jacobian
     * @param ctan the 3-d vector for the time derivative coefficients
     * @param t_subelmtid the id of the sub element
     * @param t_materialscontainer_old the material container for pervious step
     * @param t_materialscontainer the material container for current step
     * @param t_elmtinfo the local element info structure
     * @param t_elmtsoln the local element solution
     * @param t_shp the local shape function structure
     * @param K the local K matrix
     * @param R the local residual vector
     */
    void runBulkElmtLibs(const FECalcType &t_calctype,const double (&ctan)[3],
                         const int &t_subelmtid,
                         const MaterialsContainer &t_materialscontainer_old,
                         const MaterialsContainer &t_materialscontainer,
                         const LocalElmtInfo &t_elmtinfo,
                         const LocalElmtSolution &t_elmtsoln,
                         const LocalShapeFun &t_shp,
                         MatrixXd &K,
                         VectorXd &R);

    /**
     * get the total number of bulk element blocks
     */
    inline int getBulkElmtBlocksNum()const{return m_ElmtBlockNum;}

    /**
     * get the i-th bulk element block
     * @param i integer for the block index, start from 1
     */
    inline ElmtBlock getIthBulkElmtBlock(const int &i)const{
        if(i<1||i>m_ElmtBlockNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your bulk element block list");
            MessagePrinter::exitAsFem();
        }
        return m_ElmtBlockList[i-1];
    }
    /**
     * get the bulk element block vector
     */
    inline vector<ElmtBlock> getBulkElmtBlockList()const{
        return m_ElmtBlockList;
    }
    /**
     * get the sub element/modules size of i-th bulk element
     * @param i integer for the sub element index
     */
    inline int getLocalIthBulkElmtSubElmtsNum(const int &i)const{
        if(i<1||i>m_LocalBulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your local bulk elements("+to_string(m_LocalBulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        return static_cast<int>(m_LocalElemental_ElmtBlockID[i-1].size());
    }
    /**
     * get the j-th sub element index id of i-th bulk element
     * @param i integer for the bulk element index
     * @param j integer for the sub element index
     */
    inline int getLocalIthBulkElmtJthSubElmtID(const int &i,const int &j)const{
        if(i<1||i>m_LocalBulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your bulk elements("+to_string(m_LocalBulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>static_cast<int>(m_LocalElemental_ElmtBlockID[i-1].size())){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of range for your sub element("+to_string(m_LocalElemental_ElmtBlockID[i-1].size())+")");
            MessagePrinter::exitAsFem();
        }
        return m_LocalElemental_ElmtBlockID[i-1][j-1];
    }

    //**************************************************************
    //*** for information printing
    //**************************************************************
    /**
     * print out the information summary of element system
     */
    void printBulkElmtSystemInfo()const;

    /**
     * release the allocated memory
     */
    void releaseMemory();

protected:
    int m_ElmtBlockNum;/**< the total number of element block */
    vector<ElmtBlock> m_ElmtBlockList;/**< vector for all the elmt blocks read from input file */

    //****************************
    //*** for elemental info
    //****************************
    int m_GlobalBulkElmtsNum;/**< the total number of bulk elements */
    int m_LocalBulkElmtsNum;/**< the number of local bulk elements */
    vector<vector<int>> m_GlobalElemental_ElmtBlockID;/**< this vector stores the elmt block id for each single fe cell in master rank */
    vector<vector<int>> m_LocalElemental_ElmtBlockID;/**< this vector stores the elmt block id for each single element of local rank */

};