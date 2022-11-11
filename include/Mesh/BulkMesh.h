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
//+++ Date   : 2022.05.06
//+++ Purpose: the bulk mesh class of AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

/**
 * for AsFem's headers
 */
#include "Mesh/MeshData.h"
#include "Mesh/Nodes.h"


/**
 * this class implement the mesh manager for the bulk mesh
 */
class BulkMesh{
public:
    /**
     * default constructor
     */
    BulkMesh();
    /**
     * constructor for assign operator
     */
    BulkMesh(const BulkMesh &mesh);

    //*****************************************************
    //*** general settings
    //*****************************************************
    /**
     * set the x-min value of the mesh
     * @param xmin double value for xmin
     */
    void setXmin(const double &xmin){m_meshdata.m_xmin=xmin;}
    /**
     * set the x-max value of the mesh
     * @param xmax double value for xmax
     */
    void setXmax(const double &xmax){m_meshdata.m_xmax=xmax;}
    /**
     * set the y-min value of the mesh
     * @param ymin double value for ymin
     */
    void setYmin(const double &ymin){m_meshdata.m_ymin=ymin;}
    /**
     * set the y-max value of the mesh
     * @param ymax double value for ymax
     */
    void setYmax(const double &ymax){m_meshdata.m_ymax=ymax;}
    /**
     * set the z-min value of the mesh
     * @param zmin double value for zmin
     */
    void setZmin(const double &zmin){m_meshdata.m_zmin=zmin;}
    /**
     * set the z-max value of the mesh
     * @param zmax double value for zmax
     */
    void setZmax(const double &zmax){m_meshdata.m_zmax=zmax;}
    /**
     * set the number of mesh in x-axis
     * @param nx int value for nx
     */
    void setNx(const int &nx){m_meshdata.m_nx=nx;}
    /**
     * set the number of mesh in y-axis
     * @param ny int value for ny
     */
    void setNy(const int &ny){m_meshdata.m_ny=ny;}
    /**
     * set the number of mesh in z-axis
     * @param nz int value for nz
     */
    void setNz(const int &nz){m_meshdata.m_nz=nz;}
    /**
     * set the dimension of the mesh
     * @param dim int value for dimension
     */
    void setDim(const int &dim){m_meshdata.m_maxdim=dim;}
    /**
     * set the mesh type
     * @param type the mesh type of bulk mesh
     */
    void setMeshType(const MeshType &type){m_meshdata.m_bulkelmt_type=type;}

    //*****************************************************
    //*** general gettings
    //*****************************************************
    //**************************************************
    //*** for dims
    //**************************************************
    /**
     * get the max dim of bulk mesh
     */
    inline int getBulkMeshMaxDim()const{return m_meshdata.m_maxdim;}
    /**
     * get the min dim of bulk mesh
     */
    inline int getBulkMeshMinDim()const{return m_meshdata.m_mindim;}
    //**************************************************
    //*** for elements number
    //**************************************************
    /**
     * get the number of total elements(bulk elements+lower dimension elements)
     */
    inline int getTotalElmtsNum()const{return m_meshdata.m_elements;}
    /**
     * get the number of bulk elements
     */
    inline int getBulkMeshBulkElmtsNum()const{return m_meshdata.m_bulkelmts;}
    /**
     * get the number of line elements
     */
    inline int getBulkMeshLineElmtsNum()const{return m_meshdata.m_lineelmts;}
    /**
     * get the number of surface elements
     */
    inline int getBulkMeshSurfaceElmtsNum()const{return m_meshdata.m_surfaceelmts;}
    /**
     * get the number of elements via its physical name
     * @param name string for the physical name 
     */
    inline int getBulkMeshElmtsNumViaPhyName(const string &name)const{
        for(const auto &it:m_meshdata.m_phygroup_name2elmtconnvec){
            if(it.first==name){
                return static_cast<int>(it.second.size());
            }
        }
        return 0;// if the name is not there, then return 0
    }
    /**
     * get the number of bulk elements via its physical name
     * @param name string for the physical name 
     */
    inline int getBulkMeshBulkElmtsNumViaPhyName(const string &name)const{
        for(int i=0;i<m_meshdata.m_phygroups;i++){
            if(m_meshdata.m_phygroup_phynamevec[i]==name && 
               m_meshdata.m_phygroup_dimvec[i]==m_meshdata.m_maxdim){
                return m_meshdata.m_phygroup_elmtnumvec[i];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find elements number for phyname="+name+" in your Mesh class");
        MessagePrinter::exitAsFem();
        return 0;
    }
    /**
     * get i-th bulk elements global id via the physical group name
     * @param name string for the physical name 
     * @param i the local element index number, start from 1
     */
    inline int getBulkMeshIthBulkElmtIDViaPhyName(const string &name,const int &i)const{
        for(const auto &it:m_meshdata.m_phygroup_name2bulkelmtidvec){
            if(name==it.first){
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range in the bulk element set (phyname="+name+") in your Mesh class");
                    MessagePrinter::exitAsFem();
                }
                return it.second[i-1];
            }
        }
        MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range in the bulk element set (phyname="+name+") in your Mesh class");
        MessagePrinter::exitAsFem();
        return -1;
    }
    /**
     * get the vtk cell type of bulk elements
     */
    inline int getBulkMeshBulkElmtVTKCellType()const{return m_meshdata.m_bulkelmt_vtktype;}
    /**
     * get the type name of bulk elements
     */
    inline string getBulkMeshBulkElmtMeshTypeName()const{return m_meshdata.m_bulkelmt_typename;}
    /**
     * get the mesh type of bulk elements
     */
    inline MeshType getBulkMeshBulkElmtMeshType()const{return m_meshdata.m_bulkelmt_type;}
    /**
     * get the mesh type of line elements
     */
    inline MeshType getBulkMeshLineElmtMeshType()const{return m_meshdata.m_lineelmt_type;}
    /**
     * get the mesh type of surface elements
     */
    inline MeshType getBulkMeshSurfaceElmtMeshType()const{return m_meshdata.m_surfaceelmt_type;}
    /**
     * get the mesh order of bulk elements
     */
    inline int getBulkMeshBulkElmtOrder()const{return m_meshdata.m_order;}
    /**
     * get the dim of elements via its physical name
     * @param name string for the physical name 
     */
    inline int getBulkMeshElmtDimViaPhyName(const string &name)const{
        for(const auto &it:m_meshdata.m_phygroup_name2dimvec){
            if(it.first==name){
                return it.second;
            }
        }
        MessagePrinter::printErrorTxt("can\'t find bulk element dim for phyname="+name+" in BulkMesh class");
        MessagePrinter::exitAsFem();
        return 0;// if the name is not there, then return 0
    }
    //**************************************************
    //*** for mesh data
    //**************************************************
    /**
     * get the reference of mesh data
     */
    inline MeshData& getBulkMeshMeshDataRef(){return m_meshdata;}
    /**
     * get the copy of mesh data
     */
    inline MeshData getBulkMeshMeshDataRef()const{return m_meshdata;}
    //**************************************************
    //*** for nodes number
    //**************************************************
    /**
     * get the total nodes number of bulk element
     */
    inline int getBulkMeshNodesNum()const{return m_meshdata.m_nodes;}
    /**
     * get the number of nodes per bulk element
     */
    inline int getBulkMeshNodesNumPerBulkElmt()const{return m_meshdata.m_nodesperbulkelmt;}
    /**
     * get the number of nodes per line element
     */
    inline int getBulkMeshNodesNumPerLineElmt()const{return m_meshdata.m_nodesperlineelmt;}
    /**
     * get the number of nodes per surface element
     */
    inline int getBulkMeshNodesNumPerSurfaceElmt()const{return m_meshdata.m_nodespersurfaceelmt;}
    /**
     * get the nodes number of i-th bulk element
     * @param i i-th bulk element
     */
    inline int getBulkMeshIthBulkElmtNodesNum(const int &i)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i is out of range for your bulk elements");
            MessagePrinter::exitAsFem();
        }
        return static_cast<int>(m_meshdata.m_bulkelmt_connectivity[i-1].size());
    }
    /**
     * get the nodes number of i-th element via its physical name
     * @param i i-th element
     */
    inline int getBulkMeshIthElmtNodesNumViaPhyName(const string &name,const int &i)const{
        for(const auto &it:m_meshdata.m_phygroup_name2elmtconnvec){
            if(it.first==name){
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i is out of range for your elements(phyname="+name+")");
                    MessagePrinter::exitAsFem();
                }
                return static_cast<int>(it.second[i-1].size());
            }
        }
        MessagePrinter::printErrorTxt("can\'t find nodes number for your elements (phyname="+name+")");
        MessagePrinter::exitAsFem();
        return 0;
    }
    /**
     * get the nodes number of i-th element via its physical name
     * @param name physical name
     * @param i i-th element
     * @param j j-th node id
     */
    inline int getBulkMeshIthElmtJthNodeIDViaPhyName(const string &name,const int &i,const int &j)const{
        for(const auto &it:m_meshdata.m_phygroup_name2elmtconnvec){
            if(it.first==name){
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i is out of range for your elements(phyname="+name+")");
                    MessagePrinter::exitAsFem();
                }
                if(j<1||j>static_cast<int>(it.second[i-1].size())){
                    MessagePrinter::printErrorTxt("j is out of range for your elements(phyname="+name+")");
                    MessagePrinter::exitAsFem();
                }
                
                return it.second[i-1][j-1];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find node id for your elements (phyname="+name+")");
        MessagePrinter::exitAsFem();
        return -1;
    }
    /**
     * get the j-th node id of i-th element
     * @param i i-th element
     * @param j j-th node index
     */
    inline int getBulkMeshIthBulkElmtJthNodeID(const int &i,const int &j)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i is out of range for bulk element(n="+to_string(m_meshdata.m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>m_meshdata.m_nodesperbulkelmt){
            MessagePrinter::printErrorTxt("j is out of range for node index(j<1 or j>"+to_string(m_meshdata.m_nodesperbulkelmt)+")");
            MessagePrinter::exitAsFem();
        }
        return m_meshdata.m_bulkelmt_connectivity[i-1][j-1];
    }
    /**
     * get the i-th element's connectivity info, index start from 0
     * @param i the element index, start from 1
     */
    inline void getBulkMeshIthBulkElmtConnectivity0(const int &i,vector<int> &t_conn)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk element(n="+to_string(m_meshdata.m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=1;j<=m_meshdata.m_nodesperbulkelmt;j++){
            t_conn[j-1]=getBulkMeshIthBulkElmtJthNodeID(i,j)-1;
        }
    }
    /**
     * get the i-th element's connectivity info, index start from 1
     * @param i the element index, start from 1
     */
    inline void getBulkMeshIthBulkElmtConnectivity(const int &i,vector<int> &t_conn)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range for bulk element(n="+to_string(m_meshdata.m_bulkelmts)+")");
            MessagePrinter::exitAsFem();
        }
        for(int j=1;j<=m_meshdata.m_nodesperbulkelmt;j++){
            t_conn[j-1]=getBulkMeshIthBulkElmtJthNodeID(i,j);
        }
    }
    /**
     * get the i-th node id of the required (by name) nodeset
     * @param setname the physical name of the nodeset
     * @param i i-th node index
     */
    inline int getBulkMeshIthNodeIDViaNodeSetName(const string setname,const int &i)const{
        for(const auto &it:m_meshdata.m_nodephygroup_name2nodeidvec){
            if(it.first==setname){
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range("+to_string(static_cast<int>(it.second.size()))+") for a nodeset");
                    MessagePrinter::exitAsFem();
                }
                return it.second[i-1];
            }
        }
        MessagePrinter::printErrorTxt("can\'t find node id for your node set (phyname="+setname+")");
        MessagePrinter::exitAsFem();
        return 0;
    }
    /**
     * get the node id of the required (by name) nodeset
     * @param setname the physical name of the nodeset
     */
    inline vector<int> getBulkMeshNodeIDsViaNodeSetName(const string setname)const{
        for(const auto &it:m_meshdata.m_nodephygroup_name2nodeidvec){
            if(it.first==setname){
                return it.second;
            }
        }
        MessagePrinter::printErrorTxt("can\'t find node ids for your node set (phyname="+setname+")");
        MessagePrinter::exitAsFem();
        return vector<int>(0);
    }
    /**
     * get the number of nodes via its nodal physical name
     * @param name string for the nodal physical name
     */
    inline int getBulkMeshNodesNumViaNodeSetName(const string &name)const{
        for(const auto &it:m_meshdata.m_nodephygroup_name2nodeidvec){
            if(it.first==name){
                return static_cast<int>(it.second.size());
            }
        }
        MessagePrinter::printErrorTxt("can\'t find node number for your node set (phyname="+name+")");
        MessagePrinter::exitAsFem();
        return 0;
    }
    //**************************************************
    //*** for nodes coordinates
    //**************************************************
    /**
     * get the i-th node's j-th coordinate, this is the initial coordinate
     * @param i integer for i-th node
     * @param j integer for j-th coordinate
     */
    inline double getBulkMeshIthNodeJthCoord0(const int &i,const int &j)const{
        if(i<1||i>m_meshdata.m_nodes){
            MessagePrinter::printErrorTxt("i is out of range for total nodes(n="+to_string(m_meshdata.m_nodes)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>3){
            MessagePrinter::printErrorTxt("j is out of range for node coordinate(j<1 or j>3)");
            MessagePrinter::exitAsFem();
        }
        return m_meshdata.m_nodecoords0[(i-1)*3+j-1];
    }
    /**
     * get the i-th node's j-th coordinate, this is the currrent coordinate
     * @param i integer for i-th node
     * @param j integer for j-th coordinate
     */
    inline double getBulkMeshIthNodeJthCoord(const int &i,const int &j)const{
        if(i<1||i>m_meshdata.m_nodes){
            MessagePrinter::printErrorTxt("i is out of range for total nodes(n="+to_string(m_meshdata.m_nodes)+")");
            MessagePrinter::exitAsFem();
        }
        if(j<1||j>3){
            MessagePrinter::printErrorTxt("j is out of range for node coordinate(j<1 or j>3)");
            MessagePrinter::exitAsFem();
        }
        return m_meshdata.m_nodecoords[(i-1)*3+j-1];
    }
    /**
     * get the init coordinates of i-th bulk element, the nodes class must initilized before use!
     * @param i integer for the index of i-th element
     * @param nodes nodes class stores their coordinates
     */
    inline void getBulkMeshIthBulkElmtNodeCoords0(const int &i,Nodes &nodes)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i is out of range for bulk elements");
            MessagePrinter::exitAsFem();
        }
        int j,iInd;
        for(j=1;j<=static_cast<int>(m_meshdata.m_bulkelmt_connectivity[i-1].size());j++){
            iInd=m_meshdata.m_bulkelmt_connectivity[i-1][j-1];
            nodes(j,1)=getBulkMeshIthNodeJthCoord0(iInd,1);
            nodes(j,2)=getBulkMeshIthNodeJthCoord0(iInd,2);
            nodes(j,3)=getBulkMeshIthNodeJthCoord0(iInd,3);
        }
    }
    /**
     * get the current coordinates of i-th bulk element, the nodes class must initilized before use!
     * @param i integer for the index of i-th element
     * @param nodes nodes class stores their coordinates
     */
    inline void getBulkMeshIthBulkElmtNodeCoords(const int &i,Nodes &nodes)const{
        if(i<1||i>m_meshdata.m_bulkelmts){
            MessagePrinter::printErrorTxt("i is out of range for bulk elements");
            MessagePrinter::exitAsFem();
        }
        int j,iInd;
        for(j=1;j<=static_cast<int>(m_meshdata.m_bulkelmt_connectivity[i-1].size());j++){
            iInd=m_meshdata.m_bulkelmt_connectivity[i-1][j-1];
            nodes(j,1)=getBulkMeshIthNodeJthCoord(iInd,1);
            nodes(j,2)=getBulkMeshIthNodeJthCoord(iInd,2);
            nodes(j,3)=getBulkMeshIthNodeJthCoord(iInd,3);
        }
    }
    /**
     * get the init coordinates of i-th element via its physical name, the nodes class must initilized before use!
     * @param name string for its physical name
     * @param i integer for the index of i-th element
     * @param nodes nodes class stores their coordinates
     */
    inline void getBulkMeshIthElmtNodeCoords0ViaPhyName(const string &name,const int &i,Nodes &nodes)const{
        bool IsValidName=false;
        int j,iInd;
        for(const auto &it:m_meshdata.m_phygroup_name2elmtconnvec){
            if(it.first==name){
                // now, we found the correct vector
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i is out of range for phyname="+name+" elements");
                    MessagePrinter::exitAsFem();
                }
                for(j=1;j<=static_cast<int>(it.second[i-1].size());j++){
                    iInd=it.second[i-1][j-1];
                    nodes(j,1)=getBulkMeshIthNodeJthCoord0(iInd,1);
                    nodes(j,2)=getBulkMeshIthNodeJthCoord0(iInd,2);
                    nodes(j,3)=getBulkMeshIthNodeJthCoord0(iInd,3);
                }
                IsValidName=true;break;
            }
        }
        if(!IsValidName){
            MessagePrinter::printErrorTxt("can\'t find the element set for phyname="+name);
            MessagePrinter::exitAsFem();
        }
    }
    /**
     * get the current coordinates of i-th element via its physical name, the nodes class must initilized before use!
     * @param name string for its physical name
     * @param i integer for the index of i-th element in current physical group, not the total one!
     * @param nodes nodes class stores their coordinates
     */
    inline void getBulkMeshIthElmtNodeCoordsViaPhyName(const string &name,const int &i,Nodes &nodes)const{
        bool IsValidName=false;
        int j,iInd;
        for(const auto &it:m_meshdata.m_phygroup_name2elmtconnvec){
            if(it.first==name){
                // now, we found the correct vector
                if(i<1||i>static_cast<int>(it.second.size())){
                    MessagePrinter::printErrorTxt("i is out of range for phyname="+name+" elements");
                    MessagePrinter::exitAsFem();
                }
                for(j=1;j<=static_cast<int>(it.second[i-1].size());j++){
                    iInd=it.second[i-1][j-1];
                    nodes(j,1)=getBulkMeshIthNodeJthCoord(iInd,1);
                    nodes(j,2)=getBulkMeshIthNodeJthCoord(iInd,2);
                    nodes(j,3)=getBulkMeshIthNodeJthCoord(iInd,3);
                }
                IsValidName=true;break;
            }
        }
        if(!IsValidName){
            MessagePrinter::printErrorTxt("can\'t find the element set for phyname="+name);
            MessagePrinter::exitAsFem();
        }
    }
    //**************************************************
    //*** for physical group
    //**************************************************
    /**
     * get the number of elemental physical group
     */
    inline int getBulkMeshPhyGroupNum()const{return m_meshdata.m_phygroups;}
    /**
     * get the dim of i-th physical group
     * @param i integer for i-th physical group
     */
    inline int getBulkMeshIthPhyGroupDim(const int &i)const{return m_meshdata.m_phygroup_dimvec[i-1];}
    /**
     * get the nodes number of i-th physical group(element), the nodes per phycial element
     * @param i integer for i-th physical group
     */
    inline int getBulkMeshIthPhyGroupElmtNodesNum(const int &i)const{return m_meshdata.m_phygroup_nodesnumperelmtvec[i-1];}
    /**
     * get the elements number of i-th physical group
     * @param i integer for i-th physical group
     */
    inline int getBulkMeshIthPhyGroupElmtsNum(const int &i)const{return m_meshdata.m_phygroup_elmtnumvec[i-1];}
    /**
     * get the phy id of i-th physical group
     * @param i integer for i-th physical group
     */
    inline int getBulkMeshIthPhyGroupPhyID(const int &i)const{return m_meshdata.m_phygroup_phyidvec[i-1];}
    /**
     * get the phy name of i-th physical group
     * @param i integer for i-th physical group
     */
    inline string getBulkMeshIthPhyGroupPhyName(const int &i)const{return m_meshdata.m_phygroup_phynamevec[i-1];}
    //******************
    //*** for nodal
    //******************
    /**
     * get the number of nodal physical group
     */
    inline int getBulkMeshNodalPhyGroupNum()const{return m_meshdata.m_nodal_phygroups;}
    /**
     * get i-th nodal physical id
     * @param i integer for i-th nodal physical group
     */
    inline int getBulkMeshIthNodalPhyGroupPhyID(const int &i)const{return m_meshdata.m_nodephygroup_phyidvec[i-1];}
    /**
     * get physical name of i-th nodal physical group
     * @param i integer for i-th physical group
     */
    inline string getBulkMeshIthNodalPhyGroupPhyName(const int &i)const{return m_meshdata.m_nodephygroup_phynamevec[i-1];}
    /**
     * get nodes number of i-th nodal physical group
     * @param i integer for i-th physical group
     */
    inline int getBulkMeshIthNodalPhyGroupNodesNum(const int &i)const{return static_cast<int>(m_meshdata.m_nodephygroup_name2nodeidvec[i-1].second.size());}
    /**
     * check whether the given string name is a valid physical group name for boundary mesh
     * @param phyname the physical name of the boundary mesh
     */
    inline bool isBCElmtPhyNameValid(const string &phyname)const{
        for(int i=0;i<m_meshdata.m_phygroups;i++){
            if(m_meshdata.m_phygroup_dimvec[i]<m_meshdata.m_maxdim && 
               m_meshdata.m_phygroup_phynamevec[i]==phyname){
                return true;
            }
        }
        return false;
    }
    /**
     * check whether the given string name is a valid physical group name for bulk mesh
     * @param phyname the physical name of the bulk mesh
     */
    inline bool isBulkElmtPhyNameValid(const string &phyname)const{
        for(int i=0;i<m_meshdata.m_phygroups;i++){
            if(m_meshdata.m_phygroup_dimvec[i]==m_meshdata.m_maxdim && 
               m_meshdata.m_phygroup_phynamevec[i]==phyname){
                return true;
            }
        }
        return false;
    }
   

    //*****************************************************
    //*** printing and saving
    //*****************************************************
    /**
     * print out the bulk mesh information
     */
    void printBulkMeshInfo()const;
    /**
     * print out the bulk mesh information with more details
     */
    void printBulkMeshDetailedInfo()const;
    /**
     * save mesh to a vtu file, the vtu file name is determined by the inputfilename, i.e., the input file
     * is my.json, then the mesh-vtu file is 'my-mesh.vtu'
     * @param inputfilename the string name of your input file
     */
    void saveBulkMesh2VTU(const string &inputfilename)const;

    /**
     * release memory
     */
    void releaseMemory();

protected:
    MeshData m_meshdata;/**< mesh data for nodal coordinates, elemental connecitivty, etc. */
    
};