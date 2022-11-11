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
//+++ Date   : 2022.05.09
//+++ Purpose: the bulkdof manager in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

/**
 * for AsFem's header files
 */
#include "Utils/MessagePrinter.h"
#include "Mesh/Mesh.h"
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
     * @param t_mesh the mesh class
     * @param t_elmtSystem the element class
     */
    void createBulkDofsMap(const Mesh &t_mesh,const ElmtSystem &t_elmtSystem);
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
        return m_maxdofs_pernode;
    }
    /**
     * get max dofs per element
     */
    inline int getMaxDofsPerElmt()const{
        return m_maxdofs_perelmt;
    }
    /**
     * get total dofs
     */
    inline int getTotalDofs()const{
        return m_total_dofs;
    }
    /**
     * get active dofs
     */
    inline int getActiveDofs()const{
        return m_active_dofs;
    }
    /**
     * get nodes
     */
    inline int getNodesNum()const{
        return m_nodes;
    }
    /**
     * get bulk elements
     */
    inline int getBulkElmtsNum()const{
        return m_bulkelmts;
    }
    /**
     * get the i-th dof id
     * @param i integer for the i-th dof
     */
    inline int getIthDofID(const int &i)const{
        if(i<1||i>m_maxdofs_pernode){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your dofs id list");
            MessagePrinter::exitAsFem();
        }
        return m_dof_idlist[i-1];
    }
    /**
     * get the i-th dof name
     * @param i integer for the i-th dof
     */
    inline string getIthDofName(const int &i)const{
        if(i<1||i>m_maxdofs_pernode){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for your dofs name list");
            MessagePrinter::exitAsFem();
        }
        return m_dof_namelist[i-1];
    }
    /**
     * get dof's id via its name
     * @param dofname string for the name of the input dof
     */
    inline int getDofIDViaName(const string &dofname)const{
        for(int i=0;i<m_maxdofs_pernode;i++){
            if(dofname==m_dof_namelist[i]){
                return m_dof_idlist[i];
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
        if(i<1||i>m_nodes){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of nodes' range(="+to_string(m_nodes)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_maxdofs_pernode){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of nodes' max dof num(="+to_string(m_maxdofs_pernode)+")");
            MessagePrinter::exitAsFem();
        }
        return m_nodal_dofids[i-1][j-1];
    }
    /**
     * get i-th node's j-th dof id, start from 0
     * @param i integer for i-th node
     * @param j integer for j-th dof id
     */
    inline int getIthNodeJthDofID0(const int &i,const int &j)const{
        if(i<1||i>m_nodes){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of nodes' range(="+to_string(m_nodes)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_maxdofs_pernode){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of nodes' max dof num(="+to_string(m_maxdofs_pernode)+")");
            MessagePrinter::exitAsFem();
        }
        return m_nodal_dofids[i-1][j-1]-1;
    }
    /**
     * get i-th elmt's j-th dof id
     * @param i integer for i-th elmt
     * @param j integer for j-th dof id
     */
    inline int getIthBulkElmtJthDofID(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>static_cast<int>(m_elmt_dofids[i-1].size())){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of elmts' dof num(="+to_string(m_elmt_dofids[i-1].size())+")");
            MessagePrinter::exitAsFem();
        }
        return m_elmt_dofids[i-1][j-1];
    }
    /**
     * get i-th elmt's j-th dof id
     * @param i integer for i-th elmt
     * @param j integer for j-th dof id
     */
    inline int getIthBulkElmtJthDofID0(const int &i,const int &j)const{
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>static_cast<int>(m_elmt_dofids[i-1].size())){
            MessagePrinter::printErrorTxt("j="+to_string(j)+" is out of elmts' dof num(="+to_string(m_elmt_dofids[i-1].size())+")");
            MessagePrinter::exitAsFem();
        }
        return m_elmt_dofids[i-1][j-1]-1;
    }
    /**
     * get i-th elmt's dofs
     * @param i integer for i-th elmt
     * @param elmtdofs vector for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs(const int &i,vector<int> &elmtdofs)const{
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        bool AllDofsAreZero;AllDofsAreZero=true;
        for(int j=0;j<static_cast<int>(m_elmt_dofids[i-1].size());j++){
            elmtdofs[j]=m_elmt_dofids[i-1][j];
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
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        bool AllDofsAreZero;AllDofsAreZero=true;
        for(int j=0;j<static_cast<int>(m_elmt_dofids[i-1].size());j++){
            elmtdofs[j]=m_elmt_dofids[i-1][j]-1;
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
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=0;j<static_cast<int>(m_elmt_dofids[i-1].size());j++){
            elmtdofs[j]=m_elmt_dofids[i-1][j];
        }
    }
    /**
     * get i-th elmt's dofs, start from 0
     * @param i integer for i-th elmt
     * @param elmtdofs integer pointer for current element's dof ids
     */
    inline void getIthBulkElmtDofIDs0(const int &i,int *elmtdofs)const{
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=0;j<static_cast<int>(m_elmt_dofids[i-1].size());j++){
            elmtdofs[j]=m_elmt_dofids[i-1][j]-1;
        }
    }
    /**
     * get i-th elemental dofs number
     * @param i integer for i-th elmt
     */
    inline int getIthBulkElmtDofsNum(const int &i)const{
        if(i<1||i>m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of elmts' range(="+to_string(m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        return static_cast<int>(m_elmt_dofids[i-1].size());
    }
    /**
     * check if the input dof name is a valid name
     * @param dofname the input string of dof
     */
    bool isValidDofName(const string &dofname)const{
        for(const auto &it:m_dof_namelist){
            if(dofname==it){
                return true;
            }
        }
        return false;
    }
    /**
     * get the max nonzeros of each row
     */
    inline int getMaxNNZ()const{
        return m_maxnnz;
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
    vector<string> m_dof_namelist;/**< vector for the name of dofs of each node */
    vector<int> m_dof_idlist;/**< vecotr for the related dof id of each dof(name) */

    int m_bulkelmts;/**< for the number of bulk elements */
    int m_nodes;/**< for the number of nodes */
    int m_maxdofs_pernode;/**< for the maximum dofs of one single node */
    int m_maxdofs_perelmt;/**< for the maximum dofs of one single element */
    int m_total_dofs;/**< for the total dofs */
    int m_active_dofs;/**< for the active dofs */
    vector<vector<int>> m_elmt_dofids;/**< this vector stores the dofs id of each element */
    vector<vector<int>> m_nodal_dofids;/**< this vector stores the dof ids of each node */

    int m_maxnnz;/**< for the maximum nonzeros */

};