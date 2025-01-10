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
//+++ Date   : 2022.05.09
//+++ Purpose: the bulkdof manager in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

/**
 * for AsFem's header files
 */
#include "Utils/MessagePrinter.h"
#include "FECell/FECell.h"
#include "ElmtSystem/ElmtSystem.h"

using std::vector;
using std::string;


/**
 * This class implements the manager for the degree of freedoms (DoFs) of volume mesh
 */
class BulkDofHandler{
public:
    /**
     * constructor
     */
    BulkDofHandler();
    /**
     * deconstructor
     */
    ~BulkDofHandler();

    //***************************************************
    //*** general settings
    //***************************************************
    /**
     * initialize the dofhandler class, resize all the vectors
     */
    void init();
    /**
     * create the bulk elements' dofs map
     * @param t_fecell the fe cell class
     * @param t_elmtSystem the element class
     */
    void createBulkDofsMap(FECell &t_fecell,const ElmtSystem &t_elmtSystem);
    /**
     * add dof name to the namelist, it must be unique and non-duplicated name
     * @param dofname string for the name of one single dof
     */
    void addDofName2List(const string &dofname);

    //***************************************************
    //*** general gettings
    //***************************************************
    /**
     * get max dofs per node
     */
    inline int getMaxDofsPerNode()const{
        return m_MaxDofsPerNode;
    }
    /**
     * get max dofs per element
     */
    inline int getMaxDofsPerElmt()const{
        return m_MaxDofsPerElmt;
    }
    /**
     * get total dofs
     */
    inline int getTotalDofs()const{
        return m_TotalDofs;
    }
    /**
     * get active dofs
     */
    inline int getActiveDofs()const{
        return m_ActiveDofs;
    }
    /**
     * get nodes
     */
    inline int getNodesNum()const{
        return m_NodesNum;
    }
    /**
     * get bulk elements
     */
    inline int getBulkElmtsNum()const{
        return m_BulkElmtsNum;
    }
    /**
     * get the i-th dof id
     * @param i integer for the i-th dof
     */
    inline int getIthDofID(const int &i)const{
        if(i<1||i>m_MaxDofsPerNode){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your dofs id list");
            MessagePrinter::exitAsFem();
        }
        return m_DofIDList[i-1];
    }
    /**
     * get the i-th dof name
     * @param i integer for the i-th dof
     */
    inline string getIthDofName(const int &i)const{
        if(i<1||i>m_MaxDofsPerNode){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your dofs name list");
            MessagePrinter::exitAsFem();
        }
        return m_DofNameList[i-1];
    }
    /**
     * get dof's id via its name
     * @param dofname string for the name of the input dof
     */
    inline int getDofIDViaName(const string &dofname)const{
        for(int i=0;i<m_MaxDofsPerNode;i++){
            if(dofname==m_DofNameList[i]){
                return m_DofIDList[i];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find dof name="+dofname+", please check your dof name");
        MessagePrinter::exitAsFem();
        return 0;
    }
    
    /**
     * get i-th node's j-th dof id
     * @param i integer for i-th node
     * @param j integer for j-th dof id
     */
    inline int getIthNodeJthDofID(const int &i,const int &j)const{
        if(i<1||i>m_NodesNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of nodes' range(="+to_string(m_NodesNum)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_MaxDofsPerNode){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of nodes' max dof num(="+to_string(m_MaxDofsPerNode)+")");
            MessagePrinter::exitAsFem();
        }
        return m_NodalDofIDs_Global[i-1][j-1];
    }
    /**
     * get i-th node's j-th dof id, start from 0
     * @param i integer for i-th node
     * @param j integer for j-th dof id
     */
    inline int getIthNodeJthDofID0(const int &i,const int &j)const{
        if(i<1||i>m_NodesNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of nodes' range(="+to_string(m_NodesNum)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_MaxDofsPerNode){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of nodes' max dof num(="+to_string(m_MaxDofsPerNode)+")");
            MessagePrinter::exitAsFem();
        }
        return m_NodalDofIDs_Global[i-1][j-1]-1;
    }
    /**
     * get i-th elmt's j-th dof id
     * @param i integer for i-th elmt
     * @param j integer for j-th dof id
     */
    inline int getIthBulkElmtJthDofID(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>static_cast<int>(m_ElementalDofIDs_Global[i-1].size())){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of elmts' dof num(="+to_string(m_ElementalDofIDs_Global[i-1].size())+")");
            MessagePrinter::exitAsFem();
        }
        return m_ElementalDofIDs_Global[i-1][j-1];
    }
    /**
     * get i-th elmt's j-th dof id
     * @param i integer for i-th elmt
     * @param j integer for j-th dof id
     */
    inline int getIthBulkElmtJthDofID0(const int &i,const int &j)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>static_cast<int>(m_ElementalDofIDs_Global[i-1].size())){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of elmts' dof num(="+to_string(m_ElementalDofIDs_Global[i-1].size())+")");
            MessagePrinter::exitAsFem();
        }
        return m_ElementalDofIDs_Global[i-1][j-1]-1;
    }
    /**
     * get i-th elmt's dofs
     * @param i integer for i-th elmt
     * @param elmtdofs vector for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs(const int &i,vector<int> &elmtdofs)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        bool AllDofsAreZero;AllDofsAreZero=true;
        for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Global[i-1].size());j++){
            elmtdofs[j]=m_ElementalDofIDs_Global[i-1][j];
            if(elmtdofs[j]>0) AllDofsAreZero=false;
        }
        if(AllDofsAreZero){
            MessagePrinter::printErrorTxt("element-"+to_string(i)+"'s dofs are all zeros, please check your code");
            MessagePrinter::exitAsFem();
        }
    }
    /**
     * get i-th elmt's dofs, start from 0
     * @param i integer for i-th elmt
     * @param elmtdofs vector for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs0(const int &i,vector<int> &elmtdofs)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        bool AllDofsAreZero;AllDofsAreZero=true;
        for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Global[i-1].size());j++){
            elmtdofs[j]=m_ElementalDofIDs_Global[i-1][j]-1;
            if(elmtdofs[j]>0) AllDofsAreZero=false;
        }
        if(AllDofsAreZero){
            MessagePrinter::printErrorTxt("element-"+to_string(i)+"'s dofs are all zeros, please check your code");
            MessagePrinter::exitAsFem();
        }
    }
    /**
     * get i-th elmt's dofs
     * @param i integer for i-th elmt
     * @param elmtdofs integer pointer for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs(const int &i,int *elmtdofs)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Global[i-1].size());j++){
            elmtdofs[j]=m_ElementalDofIDs_Global[i-1][j];
        }
    }
    /**
     * get i-th elmt's dofs, start from 0
     * @param i integer for i-th elmt
     * @param elmtdofs integer pointer for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs0(const int &i,int *elmtdofs)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Global[i-1].size());j++){
            elmtdofs[j]=m_ElementalDofIDs_Global[i-1][j]-1;
        }
    }
    /**
     * get i-th elemental dofs number
     * @param i integer for i-th elmt
     */
    inline int getIthBulkElmtDofsNum(const int &i)const{
        if(i<1||i>m_BulkElmtsNum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_BulkElmtsNum)+")");
            MessagePrinter::exitAsFem();
        }
        return static_cast<int>(m_ElementalDofIDs_Global[i-1].size());
    }

    /**
     * get the local bulk element's number
     */
    inline int getLocalBulkElmtsNum()const{
        return static_cast<int>(m_ElementalDofIDs_Local.size());
    }
    /**
     * get i-th local element's dof ids
     * @param i integer for i-th local elmt
     */
    inline vector<int> getIthLocalBulkElmtDofIDs(const int &i)const{
        if(i<1||i>static_cast<int>(m_ElementalDofIDs_Local.size())){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of local element's range for dof map");
            MessagePrinter::exitAsFem();
        }
        return m_ElementalDofIDs_Local[i-1];
    }
    /**
     * get i-th local element's dof ids
     * @param i integer for i-th local elmt
     * @param ids integer vector that stores the local element dof ids
     */
    inline void getIthLocalBulkElmtDofIDs(const int &i,vector<int> &ids)const{
        if(i<1||i>static_cast<int>(m_ElementalDofIDs_Local.size())){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of local element's range for dof map");
            MessagePrinter::exitAsFem();
        }
        for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Local[i-1].size());j++){
            ids[j]=m_ElementalDofIDs_Local[i-1][j];
        }
    }
    /**
     * check if the input dof name is a valid name
     * @param dofname the input string of dof
     */
    bool isValidDofName(const string &dofname)const{
        for(const auto &it:m_DofNameList){
            if(dofname==it){
                return true;
            }
        }
        return false;
    }
    /**
     * get the max nonzeros of dof maps
     */
    inline int getMaxNNZ()const{
        return m_MaxNNZ;
    }
    /**
     * get the max nonzeros of each row
     */
    inline int getMaxRowNNZ()const{
        return m_MaxRowNNZ;
    }


    /**
     * release the memory
     */
    void releaseMemory();


    //***************************************************
    //*** printing
    //***************************************************
    /**
     * print out the basic info of bulk dofs
     */
    void printBulkDofsInfo()const;
    /**
     * print out the elemental dofs map, be careful, you'll get lots of message
     * @param flag true for file writing, otherwise show message in the terminal
     */
    void printBulkElementalDofsInfo(const bool &flag=false)const;


protected:
    vector<string> m_DofNameList;/**< vector for the name of dofs of each node */
    vector<int> m_DofIDList;/**< vecotr for the related dof id of each dof(name) */

    int m_BulkElmtsNum;/**< for the number of bulk elements */
    int m_NodesNum;/**< for the number of nodes */
    int m_MaxDofsPerNode;/**< for the maximum dofs of one single node */
    int m_MaxDofsPerElmt;/**< for the maximum dofs of one single element */
    int m_TotalDofs;/**< for the total dofs */
    int m_ActiveDofs;/**< for the active dofs */
    vector<vector<int>> m_ElementalDofIDs_Global;/**< this vector stores the dofs id of each element */
    vector<vector<int>> m_ElementalDofIDs_Local;/**< this vector stores the dofs id of each local element */
    vector<vector<int>> m_NodalDofIDs_Global;/**< this vector stores the dof ids of each node */

    int m_MaxRowNNZ;/**< for the maximum nonzeros of row*/
    int m_MaxNNZ;/**< for the maximum nonzeros of the dof map */

};