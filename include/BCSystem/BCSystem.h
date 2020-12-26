//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.10
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) dirichlet bc, i.e. displacement, temperature ...
//+++               2) neuman bc, i.e. flux, force
//+++               3) robin bc as well as user-defined-bc (ubc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>


//******************************************
//*** for AsFem own header
//******************************************
#include "Utils/MessagePrinter.h"

#include "BCSystem/BCBlock.h"
#include "BCSystem/BCType.h"

#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"
// #include "DofHandler/DofHandler.h" // this line must be comment out to get rid of circular include issue from DofHandler class !!!
#include "FE/FE.h"

#include "Utils/Vector3d.h"

using namespace std;

class Mesh;
class DofHandler;

class BCSystem{
public:
    BCSystem();
    void AddBCBlock2List(BCBlock &bcblock);

    void InitBCSystem(Mesh &mesh);

    inline int GetBCBlockNums()const{return _nBCBlocks;}
    inline BCBlock GetIthBCBlock(const int &i)const{return _BCBlockList[i-1];}

    //**************************************************************
    //*** add basic settings
    //**************************************************************
    void SetBCPenaltyFactor(const double &factor){_PenaltyFactor=factor;}


    //**************************************************************
    //*** Apply boundary conditions
    //**************************************************************
    void ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const double &t,const double (&ctan)[2],Vec &U,Mat &A,Vec &RHS);
    void ApplyInitialBC(Mesh &mesh,DofHandler &dofHandler,const double &t,Vec &U);

    void PrintBCSystemInfo()const;

private:
    //**************************************************************
    //*** some basic get functions
    //**************************************************************
    inline BCType GetIthBCBlockBCType(const int &i)const{return _BCBlockList[i-1]._BCType;}
    inline vector<string> GetIthBCBlockNameVec(const int &i)const{return _BCBlockList[i-1]._BoundaryNameList;}
    inline double GetIthBCBlockBCValue(const int &i)const{return _BCBlockList[i-1]._BCValue;}

    //**************************************************************
    //*** for different boundary conditions
    //**************************************************************
    void ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,const int &dofindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &U,Mat &K,Vec &RHS);
    void ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const int &dofindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &RHS);


private:
    int _nBCBlocks;
    vector<BCBlock> _BCBlockList;

private:
    double _PenaltyFactor;
    int _nBCDim,_nBulkDim,_nNodesPerBCElmt;
    PetscMPIInt _rank,_size;

    //*******************************
    //*** for boundary integration
    //*******************************
    double _xi,_eta,_JxW;
    Nodes _elNodes;
    double _xs[3][3],_dist;

    Vector3d _normals;
    Vec _Useq;
    VecScatter _scatteru;

};