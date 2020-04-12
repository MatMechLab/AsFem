//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::RunMateLib(const MateType &imate,const int &mateindex,const int &nDim,
                const double &t,const double &dt,
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                vector<double> &gpHist,const vector<double> &gpHistOld){
    switch (imate){
    case MateType::NullMate:
        return;
        break;
    case MateType::ConstPoissonMate:
        ConstPoissonMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::NonLinearPoissonMate:
        NonlinearPoissonMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::ConstDiffusionMate:
        ConstDiffusionMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::CahnHilliardMate:
        CahnHilliardMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::TensorCahnHilliardMate:
        TensorCahnHilliardMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::LinearElasticCahnHilliardMate:
        LinearElasticCahnHilliardMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::LinearElasticMate:
        LinearElasticMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::LinearThermalMechanicsMate:
        LinearThermalElasticMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::MieheLinearElasticFracMate:
        MieheLinearElasticMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::CohesivePFFracMate:
        CohesivePFFractureMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::AnisoLinearElasticPhaseFieldFracMate:
        AnisoPFFractureMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::MieheNeoHookeanFracMate:
        MieheNeoHookeanMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::BordenLinearElasticFracMate:
        BordenLinearElasticMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::SaintVenantMate:
        SaintVenantMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::NeoHookeanMate:
        NeoHookeanMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::DendriteMate:
        DendriteMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::CurrentThermalMate:
        CurrentThermalMaterial(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::User1Mate:
        UserMaterial1(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::User2Mate:
        UserMaterial2(nDim,t,dt,_MateBlockList[mateindex-1]._Params,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    default:
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported material type in umate                   !!!   ***\n");
        Msg_AsFem_Exit();
        break;
    }
}