//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2023.12.30
//+++ Function: finite element cell management class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>


#include "Utils/MessagePrinter.h"
#include "FECell/FECellData.h"


using std::vector;
using std::map;
using std::make_pair;
using std::string;

/**
 * This class defines and manages the basic/single finite element cell data structure,
 * users can use this class for either the elemental loop or system vector/matrix assemble.
*/
class FECell{
public:
    /**
     * Constructor
    */
    FECell();
    /**
     * Constructor
    */
    FECell(const FECell &a);

    /**
     * setup mesh geometry info for 1d case
     * @param nx mesh number along x-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype);
    /**
     * setup mesh geometry info for 2d case
     * @param nx mesh number along x-axis
     * @param ny mesh number along y-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param ymin ymin value of the domain
     * @param ymax ymax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const int &ny,
                     const double &xmin,const double &xmax,
                     const double &ymin,const double &ymax,const MeshType &meshtype);
    /**
     * setup mesh geometry info for 3d case
     * @param nx mesh number along x-axis
     * @param ny mesh number along y-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param ymin ymin value of the domain
     * @param ymax ymax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const int &ny,const int &nz,
                     const double &xmin,const double &xmax,
                     const double &ymin,const double &ymax,
                     const double &zmin,const double &zmax,const MeshType &meshtype);
    
    //*****************************************************
    //*** general gettings
    //*****************************************************
    //**************************************************
    //*** for dims
    //**************************************************
    /**
     * Get the max dim of bulk mesh
     */
    inline int getFECellMaxDim()const{return m_CellData.MaxDim;}
    /**
     * Get the min dim of bulk mesh
     */
    inline int getFeCellMinDim()const{return m_CellData.MinDim;}


    /**
     * Get the reference of fecell data
    */
    inline FECellData& getCellDataRef(){return m_CellData;}
    /**
     * Get the copy of fecell data
    */
    inline FECellData getCellDataCopy()const{return m_CellData;}

    /**
     * Get the number of toal fe cell
     */
    inline int getFECellElmtsNum()const{return m_CellData.ElmtsNum;}
    /**
     * Get the number of bulk FE cell
     */
    inline int getFECellBulkElmtsNum()const{return m_CellData.BulkElmtsNum;}
    /**
     * Get the number of surface FE cell
     */
    inline int getFECellSurfElmtsNum()const{return m_CellData.SurfElmtsNum;}
    /**
     * Get the number of line FE cell
     */
    inline int getFECellLineElmtsNum()const{return m_CellData.LineElmtsNum;}

    /**
     * Get the nodes number of FE cell
     */
    inline int getFECellNodesNum()const{return m_CellData.NodesNum;}
    /**
     * Get the nodes number of bulk FE cell
     */
    inline int getFECellNodesNumPerBulkElmt()const{return m_CellData.NodesNumPerBulkElmt;}
    /**
     * Get the nodes number of surface FE cell
     */
    inline int getFECellNodesNumPerSurfElmt()const{return m_CellData.NodesNumPerSurfElmt;}
    /**
     * Get the nodes number of line FE cell
     */
    inline int getFECellNodesNumPerLineElmt()const{return m_CellData.NodesNumPerLineElmt;}

    /**
     * Get the copy of the local bulk mesh cell vector
     */
    inline vector<SingleMeshCell>  getLocalBulkFECellVecCopy()const{return m_CellData.MeshCell_Local;}
    /**
     * Get the reference of the local bulk mesh cell vector
     */
    inline vector<SingleMeshCell>& getLocalBulkFECellVecRef(){return m_CellData.MeshCell_Local;}

    /**
     * Get the copy of local physical id to fe cell vector mapping
     */
    inline map<int,vector<SingleMeshCell>> getLocalPhyID2MeshCellVectorMapCopy()const{
        return m_CellData.PhyID2MeshCellVectorMap_Local;
    }
    /**
     * Get the reference of local physical id to fe cell vector mapping
     */
    inline map<int,vector<SingleMeshCell>>& getLocalPhyID2MeshCellVectorMapRef(){
        return m_CellData.PhyID2MeshCellVectorMap_Local;
    }

    /**
     * Get the copy of global physical id to fe cell vector mapping
     */
    inline map<int,vector<SingleMeshCell>> getGlobalPhyID2MeshCellVectorMapCopy()const{
        return m_CellData.PhyID2MeshCellVectorMap_Global;
    }
    /**
     * Get the reference of local physical id to fe cell vector mapping
     */
    inline map<int,vector<SingleMeshCell>>& getGlobalPhyID2MeshCellVectorMapRef(){
        return m_CellData.PhyID2MeshCellVectorMap_Global;
    }

    /**
     * Get the copy of the local physical name to fe cell vector mapping
     */
    inline map<string,vector<SingleMeshCell>> getLocalPhyName2MeshCellVectorMapCopy()const{
        return m_CellData.PhyName2MeshCellVectorMap_Local;
    }
    /**
     * Get the reference of the local physical name to fe cell vector mapping
     */
    inline map<string,vector<SingleMeshCell>>& getLocalPhyName2MeshCellVectorMapRef(){
        return m_CellData.PhyName2MeshCellVectorMap_Local;
    }
    /**
     * Get the copy of the global physical name to fe cell vector mapping
     */
    inline map<string,vector<SingleMeshCell>> getGlobalPhyName2MeshCellVectorMapCopy()const{
        return m_CellData.PhyName2MeshCellVectorMap_Global;
    }
    /**
     * Get the reference of the global physical name to fe cell vector mapping
     */
    inline map<string,vector<SingleMeshCell>>& getGlobalPhyName2MeshCellVectorMapRef(){
        return m_CellData.PhyName2MeshCellVectorMap_Global;
    }

    /**
     * Get the number of elements via its physical name
     * @param name string for the physical name 
     */
    inline int getFECellElmtsNumViaPhyName(const string &name)const{
        for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
            if(m_CellData.PhyNameVector_Global[i]==name){
                return m_CellData.PhyGroupElmtsNumVector_Global[i];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find elements number for phyname="+name+" in your fe cell class");
        MessagePrinter::exitAsFem();
        return 0;// if the name is not there, then return 0
    }
    /**
     * get the number of bulk elements via its physical name
     * @param name string for the physical name 
     */
    inline int getFECellBulkElmtsNumViaPhyName(const string &name)const{
        for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
            if(m_CellData.PhyNameVector_Global[i]==name &&
               m_CellData.PhyDimVector_Global[i]==m_CellData.MaxDim){
                return m_CellData.PhyGroupElmtsNumVector_Global[i];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find elements number for phyname="+name+" in your fe cell class");
        MessagePrinter::exitAsFem();
        return 0;
    }
    /**
     * Get the mesh order of bulk elements
     */
    inline int getFECellBulkMeshOrder()const{return m_CellData.MeshOrder;}
    /**
     * Get the dim of elements via its physical name
     * @param name string for the physical name 
     */
    inline int getFECellElmtDimViaPhyName(const string &name)const{
        for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
            if(m_CellData.PhyNameVector_Global[i]==name){
                return m_CellData.PhyDimVector_Global[i];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find element dim for phyname="+name+" in fe cell class");
        MessagePrinter::exitAsFem();
        return 0;// if the name is not there, then return 0
    }
    /**
     * Get the mesh type of bulk elements
     */
    inline MeshType getFECellBulkElmtMeshType()const{return m_CellData.BulkElmtMeshType;}
    /**
     * Get the mesh type of line elements
     */
    inline MeshType getFECellLineElmtMeshType()const{return m_CellData.LineElmtMeshType;}
    /**
     * Get the mesh type of surface elements
     */
    inline MeshType getFECellSurfElmtMeshType()const{return m_CellData.SurfElmtMeshType;}





    //********************************************************
    //*** for physical group info
    //********************************************************
    /**
     * Get the number of physical groups
     */
    inline int getPhysicalGroupNum()const{return m_CellData.PhyGroupNum_Global;}
    /**
     * Get the physical group dim vector copy
     */
    inline vector<int> getPhysicalGroupDimVecCopy()const{return m_CellData.PhyDimVector_Global;}
    /**
     * Get the physical group id vector copy
     */
    inline vector<int> getPhysicalGroupIDVecCopy()const{return m_CellData.PhyIDVector_Global;}
    /**
     * Get the physical group name vector copy
     */
    inline vector<string> getPhysicalGroupNameVecCopy()const{return m_CellData.PhyNameVector_Global;}
    /**
     * Get the physical group elements number vector copy
     */
    inline vector<int> getPhysicalGroupElmtsNumVecCopy()const{return m_CellData.PhyGroupElmtsNumVector_Global;}
    /**
     * Get physical id to physical name mapping 
     */
    inline map<int,string> getPhysicalID2NameMapCopy()const{return m_CellData.PhyID2NameMap_Global;}
    /**
     * Get the physical name to physical id mapping
     */
    inline map<string,int> getPhysicalName2IDMapCopy()const{return m_CellData.PhyName2IDMap_Global;}
    /**
     * Check whether the given string name is a valid physical group name for boundary mesh
     * @param phyname the physical name of the boundary mesh
     */
    inline bool isFECellBCElmtPhyNameValid(const string &phyname)const{
        for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
            if(m_CellData.PhyDimVector_Global[i]<m_CellData.MaxDim && 
               m_CellData.PhyNameVector_Global[i]==phyname){
                return true;
            }
        }
        return false;
    }
    /**
     * Check whether the given string name is a valid physical group name for bulk mesh
     * @param phyname the physical name of the bulk mesh
     */
    inline bool isFECellBulkElmtPhyNameValid(const string &phyname)const{
        for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
            if(m_CellData.PhyDimVector_Global[i]==m_CellData.MaxDim && 
               m_CellData.PhyNameVector_Global[i]==phyname){
                return true;
            }
        }
        return false;
    }


    /**
     * Save the FECell mesh to vtu file
    */
    void saveFECell2VTUFile()const;
    /**
     * Print out the summary information of FECell class, i.e., mpi-based mesh distribution, elements num, nodes num, etc.
    */
    void printSummaryInfo()const;

    void releaseMemory();

private:
    FECellData m_CellData;/**< the FE cell data structure of the whole system */


};