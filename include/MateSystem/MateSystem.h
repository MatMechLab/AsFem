//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MATESYSTEM_H
#define ASFEM_MATESYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>


#include "petsc.h"

//****************************************
//*** For AsFem's own header file
//****************************************
#include "MessagePrinter/MessagePrinter.h"

#include "MateType.h"
#include "MateBlock.h"

#include "Utils/MathUtils.h"
#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"


class MateSystem{
public:
    MateSystem();
    void RunMateLib();
    void InitMateSystem();

    //********************************************
    //*** For User-Defined-Materials (UMAT)
    //********************************************
    void RunMateLib(const MateType &imate,const int &mateindex,const int &nDim,
                    const double &t,const double &dt,
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld);

    //********************************************
    //*** Add MateBlock from input file
    //********************************************
    void AddMateBlock(MateBlock &mateblock);


    //********************************************
    //*** basic getting information functions
    //********************************************
    inline PetscInt GetMateBlocksNum()const{return _nMateBlocks;}
    MateBlock GetIthMateBlock(const PetscInt &i)const{
        return _MateBlockList[i-1];
    }

    void PrintMateSystemInfo()const;

private:
    //*********************************************
    //*** umat implementations
    //*********************************************
    void ConstPoissonMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    void NonlinearPoissonMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    
    //*********************************************
    //*** for general diffusion problem
    //*********************************************
    void ConstDiffusionMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    
    //*********************************************
    //*** for free energy of cahn-hilliard problem
    //*********************************************
    void CahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //*** tensor format
    void TensorCahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    
    //***********************************
    //*** For linear elastic material
    //***********************************
    void LinearElasticMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //***********************************
    //*** For linear elastic-thermal coupled material
    //***********************************
    void LinearThermalElasticMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //***********************************
    //*** For linear elastic-CahnHilliard  material
    //***********************************
    void LinearElasticCahnHilliardMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    
    //***************************************
    //*** For Miehe's model with linear elastic material(small deformation)
    //***************************************
    void MieheLinearElasticMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //***************************************
    //*** For cohesive phase field fracture material
    //*** taken from Y.J. Liu
    //***************************************
    void CohesivePFFractureMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //***************************************
    //*** For anisotropic linear elastic fracture
    //*** taken from Y.J. Liu
    //***************************************
     void AnisoPFFractureMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //***************************************
    //*** For Miehe's model with neo-hookean material(finite deformation)
    //***************************************
    void MieheNeoHookeanMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    
    //***************************************
    //*** For Borden's model with linear elastic material(small deformation)
    //***************************************
    void BordenLinearElasticMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //************************************
    //*** For saint venant material
    //************************************
    void SaintVenantMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //************************************
    //*** For compressive neo-hookean material
    //************************************
    void NeoHookeanMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //******************************************************
    //*** For dendrite material (only in 2d !!!)
    //******************************************************
    void DendriteMaterial(const int &nDim,const double &t,const double &dt,
                        const vector<double> InputParams,
                        const Vector3d &gpCoord,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                        vector<double> &gpHist,const vector<double> &gpHistOld);

    //************************************
    //*** For current induced thermal material
    //************************************
    void CurrentThermalMaterial(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //******************************************************
    //*** For User-Defined-Material (umat) code
    //******************************************************
    void UserMaterial1(const int &nDim,const double &t,const double &dt,
                    const vector<double> InputParams,
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld);
    //----- for 2nd user material
    void UserMaterial2(const int &nDim,const double &t,const double &dt,
                    const vector<double> InputParams,
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld);

                    
public:
    //********************************************
    //*** material paramers to be used in uel
    //********************************************
    vector<double> _ScalarMaterials;
    vector<Vector3d> _VectorMaterials;
    vector<RankTwoTensor> _Rank2Materials;
    vector<RankFourTensor> _Rank4Materials;
private:
    vector<MateBlock> _MateBlockList;
    PetscInt _nMateBlocks=0;
};

using namespace std;




#endif //ASFEM_MATESYSTEM_H