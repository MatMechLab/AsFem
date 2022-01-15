//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
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
#include "FESystem/FECalcType.h"

#include "Utils/Vector3d.h"

/**
 * For built-in and user-defined boundary condition sub-classes
 */

#include "BCSystem/DirichletBC.h"
#include "BCSystem/NeumannBC.h"
#include "BCSystem/FluxBC.h"
#include "BCSystem/NodalForceBC.h"
#include "BCSystem/CyclicDirichletBC.h"
/**
 * For user-defined dirichlet bc
 */
#include "BCSystem/User1DirichletBC.h"
#include "BCSystem/User2DirichletBC.h"
#include "BCSystem/User3DirichletBC.h"
#include "BCSystem/User4DirichletBC.h"
#include "BCSystem/User5DirichletBC.h"
/**
 * For user-defined-interated bc
 */
#include "BCSystem/User1BC.h"
#include "BCSystem/User2BC.h"
#include "BCSystem/User3BC.h"
#include "BCSystem/User4BC.h"
#include "BCSystem/User5BC.h"
#include "BCSystem/User6BC.h"
#include "BCSystem/User7BC.h"
#include "BCSystem/User8BC.h"
#include "BCSystem/User9BC.h"
#include "BCSystem/User10BC.h"

using namespace std;

class Mesh;
class DofHandler;

class BCSystem:public DirichletBC,
               public NeumannBC,
               public FluxBC,
               public NodalForceBC,
               public CyclicDirichletBC,
               // for user-defined-dirichlet-type bc
               public User1DirichletBC,
               public User2DirichletBC,
               public User3DirichletBC,
               public User4DirichletBC,
               public User5DirichletBC,
               // for user-defined-integrated-type bc
               public User1BC,
               public User2BC,
               public User3BC,
               public User4BC,
               public User5BC,
               public User6BC,
               public User7BC,
               public User8BC,
               public User9BC,
               public User10BC
{
public:
    BCSystem();
    void AddBCBlock2List(BCBlock &bcblock);

    void InitBCSystem(const Mesh &mesh);

    inline int GetBCBlockNums()const{return _nBCBlocks;}
    inline BCBlock GetIthBCBlock(const int &i)const{return _BCBlockList[i-1];}

    bool CheckAppliedBCNameIsValid(const Mesh &mesh);
    //**************************************************************
    //*** add basic settings
    //**************************************************************
    void SetBCPenaltyFactor(const double &factor){_PenaltyFactor=factor;}


    //**************************************************************
    //*** Apply boundary conditions
    //**************************************************************
    void ApplyBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const FECalcType &calctype,const double &t,const double (&ctan)[3],Vec &U,Vec &V,Mat &AMATRIX,Vec &RHS);
    
    void ApplyInitialBC(const Mesh &mesh,const DofHandler &dofHandler,const double &t,Vec &U);

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
    void ApplyDirichletBC(const FECalcType &calctype,const BCType &bctype,const vector<string> bcnamelist,const vector<int> &dofindex,const double &bcvalue,const vector<double> &params,const Mesh &mesh,const DofHandler &dofHandler,Vec &U,Mat &K,Vec &RHS);
    
    void ApplyNodalDirichletBC(const FECalcType &calctype,const BCType &bctype,const vector<string> bcnamelist,const vector<int> &dofindex,const double &bcvalue,const vector<double> &params,const Mesh &mesh,const DofHandler &dofHandler,Vec &U,Mat &K,Vec &RHS);

    //**************************************************************
    //*** for nodal type boundary conditions
    //**************************************************************
    void ApplyNodalNeumannBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const vector<int> &dofsindex,const double &bcvalue,const vector<string> &bcnamelist,Vec &RHS);

    //**************************************************************
    //*** for other general boundary conditions
    //**************************************************************
    void RunBCLibs(const FECalcType &calctype,const BCType &bctype,
            const double &bcvalue,const vector<double> &parameters,
            const Vector3d &normals,const double (&ctan)[3],
            const LocalElmtInfo &elmtinfo,
            const LocalElmtSolution &soln,
            const LocalShapeFun &shp,
            VectorXd &localR,
            MatrixXd &localK);


private:
    int _nBCBlocks;
    vector<BCBlock> _BCBlockList;

    /**
     * these three variables are used to store the local information
     */
    LocalElmtInfo _elmtinfo;
    LocalShapeFun _shp;
    LocalElmtSolution _soln;

    vector<int> _dofids;/**< the active dofs id for assemble (start from 0!!!) >*/

private:
    double _PenaltyFactor;
    int _nBCDim,_nDim,_nBulkDim,_nNodesPerBCElmt;
    PetscMPIInt _rank,_size;

    //*******************************
    //*** for boundary integration
    //*******************************
    double _xi,_eta,_JxW;
    Nodes _elNodes;
    double _xs[3][3],_dist;

    Vector3d _normals;
    VectorXd _localR;
    MatrixXd _localK;

    // for PETSc scatter
    Vec _Useq,_Vseq;
    VecScatter _scatteru,_scatterv;

};
