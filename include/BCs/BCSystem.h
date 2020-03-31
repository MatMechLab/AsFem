//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_BCSYSTEM_H
#define ASFEM_BCSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

//***********************************
//*** For AsFem's own header file
//***********************************
#include "MessagePrinter/MessagePrinter.h"
#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"

#include "BCType.h"
#include "BCBlock.h"

class Mesh;
class DofHandler;

class BCSystem{
public:
    BCSystem();
    void InitBCSystem(Mesh &mesh);
    void AddBCBlock(BCBlock &bcblock);

    void SetBCPenaltyFactor(const PetscReal &factor){_PenaltyFactor=factor;}

    void ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const PetscReal &t,const PetscReal (&ctan)[2],Mat &K,Vec &RHS,Vec &U);
    void ApplyInitialBC(Mesh &mesh,DofHandler &dofHandler,const PetscReal &t,Vec &U);
private:
    //***********************************
    //*** apply boundary conditions
    //***********************************
    void ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,Mat &K,Vec &RHS,Vec &U);
    
    void ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
                    const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,
                    Vec &RHS);

    void ApplyPressureBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
                    const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,
                    Vec &RHS);

    void ApplyPresetBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
            const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,const PetscReal (&ctan)[2],
            Mat &K,Vec &RHS,Vec &U);
            
    void ApplyUser1BC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
            const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,const PetscReal (&ctan)[2],
            Mat &K,Vec &RHS,Vec &U);


    //***********************************
    //*** some basic getting functions
    //***********************************
public:
    inline PetscInt GetBCBlocksNum()const{return _nBCBlocks;}
    BCBlock GetIthBCBlock(const PetscInt &i)const{
        return _BCBlockList[i-1];
    }
private:
    inline BCType GetIthBCBlockBCType(const PetscInt &i)const{
        return _BCBlockList[i-1]._BCType;
    }
    inline vector<string> GetIthBCBlockBCName(const PetscInt &i)const{
        return _BCBlockList[i-1]._BoundaryNameList;
    }
    inline PetscReal GetIthBCBlockBCValue(const PetscInt &i)const{
        return _BCBlockList[i-1]._BCValue;
    }


public:
    void PrintBCSystemInfo()const;


private:
    vector<BCBlock> _BCBlockList;
    int _nBCBlocks;
    PetscInt _nDim,_nNodesPerBCElmt;
    PetscReal _PenaltyFactor=1.0;
    PetscMPIInt _rank,_size;

    PetscReal _xi,_eta,_JxW;
    Nodes _elNodes;
    double _xs[3][3],_dist;
    Vector3d _normals;
    Vec _Useq;
    VecScatter _scatteru;
};

#endif // ASFEM_BCSYSTEM_H