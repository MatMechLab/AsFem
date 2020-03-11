//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_ELMTSYSTEM_H
#define ASFEM_ELMTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <vector>


#include "petsc.h"

//***********************************
//*** AsFem's own header file
//***********************************
#include "MessagePrinter/MessagePrinter.h"

#include "ElmtBlock.h"
#include "ElmtType.h"

#include "FE/FE.h"

#include "Utils/MathUtils.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"
#include "Utils/Vector3d.h"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;


class ElmtSystem{
public:
    ElmtSystem();

public:
    void RunElmtLib(const int &isw,const ElmtType &iuel,const int &nDim,const int &nNodes,
                    const double &t,const double &dt,const double (&ctan)[2],
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    const ShapeFun &shp,
                    const vector<double> &ScalarMaterials,
                    const vector<Vector3d> &VectorMaterials,
                    const vector<RankTwoTensor> &Rank2Materials,
                    const vector<RankFourTensor> &Rank4Materials,
                    vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                    MatrixXd &localK,VectorXd &localR);

public:
    //*******************************************
    //*** some setting funs
    //*******************************************
    void AddElmtBlock(ElmtBlock &elmtBlock);

    //*******************************************
    //*** some getting funs
    //*******************************************
    PetscInt GetElmtBlocksNum() const{return _nElmtBlocks;}
    ElmtBlock& GetIthElmtBlock(const PetscInt &i){
        return _ElmtBlockList[i-1];
    }
    ElmtBlock GetIthElmtBlock(const PetscInt &i)const{
        return _ElmtBlockList[i-1];
    }

private:
    //*******************************************************
    //*** For User-Defined-Elements (UEL)
    //*******************************************************
    void Poisson(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //*******************************************************
    //*** For diffusion and thermal conduct problem
    //*******************************************************
    void Diffusion(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    
    //*******************************************************
    //*** For thermal conduct problem
    //*******************************************************
    void Thermal(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //*******************************************************
    //*** For Cahn-Hilliard equation
    //*******************************************************
    void CahnHilliard(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    void TensorCahnHilliard(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    void MechanicalCahnHilliard(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //***************************
    //*** Mechanics
    //***************************
    void Mechanics(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //***************************
    //*** Thermal-Mechanics
    //***************************
    void ThermalMechanics(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //***********************************************
    //*** Miehe's phase field fracture model
    //***********************************************
    void MieheFracture(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //***********************************************
    //*** Yingjie Liu's cohesive phase field fracture model
    //***********************************************
    void CohesivePFFracture(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    //***********************************************
    //*** Miehe's phase field fracture model
    //***********************************************
    void BordenFracture(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);
    
    //***********************************************
    //*** For 2D dendrite model
    //***********************************************
    void DendriteModel(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);

    //*******************************************************
    //*** for User-Defined-Element (UEL) code
    //*******************************************************
    void UserElmt1(const int &isw,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR);

public:
    //********************************************************
    //*** print out the summary information of elmtsystem
    //********************************************************
    void PrintElmtSystemInfo()const;


private:
    vector<ElmtBlock> _ElmtBlockList;
    PetscInt _nElmtBlocks;
};


#endif // ASFEM_ELMTSYSTEM_H